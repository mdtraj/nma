from itertools import combinations_with_replacement

import scipy
from scipy.constants import Boltzmann, Avogadro

import numpy as np

import mdtraj as md

kb = 10**-3 * Boltzmann * Avogadro

__all__ = ['GNMA']


class GNMA(object):
    """A class for Gaussian Network Model Analysis (GNMA) ([IB97]_, [TH97]_)
    using MDTraj.

    See example :ref:`gnm`.

    .. [IB97] Bahar I, Atilgan AR, Erman B. Direct evaluation of thermal
       fluctuations in protein using a single parameter harmonic potential.
       *Folding & Design* **1997** 2:173-181.

    .. [TH97] Haliloglu T, Bahar I, Erman B. Gaussian dynamics of folded
       proteins. *Phys. Rev. Lett.* **1997** 79:3090-3093.
    """

    def __init__(self, mode=0, k=5., nb_cutoff=0.5, steps=10, temp=300.,
                 selection='not element H'):
        self.mode = mode
        self.k = k
        self.nb_cutoff = nb_cutoff
        self.steps = steps
        self.temp = temp
        self.selection = selection
        self.dirty = True

    def _define_interactions(self):
        self.bonded = [(i.index, j.index) for i, j in self.top.bonds]
        self.bonded += [(j.index, i.index) for i, j in self.top.bonds]
        combos = combinations_with_replacement(self.top.atoms, 2)
        self.non_bonded = [(i.index, j.index) for i, j in combos
                           if i != j and (i.index, j.index) not in self.bonded]

    def _solve(self, K):

        vals, vecs = scipy.linalg.eigh(K, turbo=True, eigvals=None)
        ind = np.argsort(vals)[1:]
        self.eigenvalues_, self.eigenvectors_ = vals[ind], vecs[:, ind]

    def _reweight(self):
        self.frequencies_ = np.sqrt(np.abs(self.eigenvalues_)) / (2. * np.pi)
        self.weighted_eigenvectors_ = (self.eigenvectors_ /
                                       np.sqrt(self.masses *
                                               np.ones_like(self.eigenvectors_
                                                            ).T
                                               ).T
                                       )
        self.amplitudes_ = np.sqrt(
            2 * self.temp * kb) / (2 * np.pi * self.frequencies_)
        self.weighted_eigenvectors_ *= self.amplitudes_

    def _animate(self, mode):
        mid = int(self.steps / 2.)
        xyz = np.zeros((self.steps + 1, self.top.n_atoms, 3))
        xyz[mid] = self.xyz[0]

        if isinstance(mode, int):
            mode = [mode]

        grad = np.atleast_2d(self.weighted_eigenvectors_[:, mode].dot(
            self.frequencies_[mode])).T
        grad *= np.ones((self.top.n_atoms, 3))

        for i in range(mid):
            xyz[i + mid + 1] = xyz[i + mid] + grad
            xyz[mid - i - 1] = xyz[mid - i] - grad

        return md.Trajectory(xyz[:-1], self.top)

    def transform(self, X):
        if self.dirty:
            self.fit(X)
            self.dirty = False

        return self._animate(self.mode)

    def _nb_force(self, x):
        return 1E-6 / (x - self.nb_cutoff).clip(1E-3, np.inf)**2.

    def fit(self, X):
        self.ind = X.top.select(self.selection)
        self._structure = X[0].atom_slice(self.ind)
        self.top = self._structure.top
        self.xyz = self._structure.xyz

        # Extract Masses
        self.masses = np.array([i.element.mass for i in self.top.atoms])
        M = np.sqrt(np.outer(self.masses, self.masses))

        # Define Non-Bonded and Bonded Interactions
        self._define_interactions()

        # Initialize Kirchoff Matrix
        self.kirchkoff_ = np.zeros((self.top.n_atoms, self.top.n_atoms))

        # Add Non-Bonded Interactions to Kirchoff Matrix
        d = md.compute_distances(self._structure,
                                 atom_pairs=self.non_bonded).ravel()
        self.kirchkoff_[list(zip(*self.non_bonded))] = - (self.k *
                                                          self._nb_force(d))
        self.kirchkoff_ += self.kirchkoff_.T

        # Add Bonded Interactions to Kirchoff Matrix
        if self.bonded:
            self.kirchkoff_[list(zip(*self.bonded))] = - 1E2 * self.k

        # Set Diagonals
        dind = np.diag_indices_from(self.kirchkoff_)
        self.kirchkoff_[dind] -= self.kirchkoff_.sum(axis=1)

        # Solve Eigenvalue Problem
        self._solve(self.kirchkoff_ / M)

        # Reweight Eigenvalues and Eigenvectors
        self._reweight()

        self.dirty = False

        return self

    def fit_transform(self, X):
        self.dirty = True
        return self.transform(X)
