from itertools import combinations_with_replacement

from scipy import linalg
import numpy as np

import networkx as nx

import mdtraj as md

__all__ = ['ANMA']


def distance_cutoff(structure, non_bonded, nb_cutoff=0.5):
    non_bonded = list(non_bonded)
    distances = md.compute_distances(structure, atom_pairs=non_bonded).ravel()
    filtered =  np.where(distances < nb_cutoff)[0]
    return {non_bonded[i] for i in filtered}


class ANMA(object):
    """Class for Anisotropic Network Model (ANM) analysis of proteins

    Parameters
    ----------
    k_b : float, default=1.
        Bonded spring constant.
    k_nb : float, default=1.
        Non-bonded spring constant.
    mode : int, default=0
        Which mode to animate.
    nb_func : func, default=None
        Function to define a non-bonded interaction which accepts arguments
        structure (mdtraj.Trajectory) and a set of non-bonded contacts.
        If `None`, defaults to a distance cutoff provided by the argument
        `nb_cutoff` (default=0.5) in nanometers.
    n_steps : int, default=50
        Number of frames to animate.
    rmsd=0.15 : float, default=0.15
        RMSD tolerance for animation in nanometers.
    selection : str, default='not element H'
        Atom selection used to compute network model. Non-selected atoms will
        act as rigid bodies.
    turbo : bool, default=True
        Use divide and conquer algorithm for generalized eigenvalue problem
        (faster but expensive in memory).
    **nb_kwargs
        Keyword arguments to pass on to `nb_func`.

    Attributes
    ----------
    eigenvalues_ : array-like, shape (n_atoms,)
        Eigenvalues of the generalized eigenproblem, in decreasing
        order.
    eigenvectors_ : array-like, shape (n_components, n_features)
        Eigenvectors of the generalized eigenproblem. The vectors
        give a set of "directions" through configuration space along
        which the system relaxes. Each eigenvector is associated with
        characteritic normal mode.
    vars_ : array-like, shape (n_atoms,)
        Reciprocal of eigenvalues_.
    trace_ : float
        Sum of eigenvalues_.

    References
    ----------
    .. [1] Doruker P, Atilgan AR, Bahar I. Dynamics of proteins predicted by
       molecular dynamics simulations and analytical approaches: Application to
       a-amylase inhibitor. *Proteins* **2000** 40:512-524.

    .. [2] Atilgan AR, Durrell SR, Jernigan RL, Demirel MC, Keskin O,
       Bahar I. Anisotropy of fluctuation dynamics of proteins with an
       elastic network model. *Biophys. J.* **2001** 80:505-515.
    """
    def __init__(self, k_b=1., k_nb=1., mode=0, n_steps=10, rmsd=0.15,
                 nb_func=None, selection='not element H',
                 turbo=True, **nb_kwargs):
        self.mode = mode
        self.k_b = k_b
        self.k_nb = k_nb
        self.n_steps = n_steps
        self.rmsd = rmsd
        self.selection = selection
        self.turbo = turbo
        self._dirty = True
        self._bonded = None
        self._non_bonded = None
        self.nb_kwargs = nb_kwargs
        self.nb_func = nb_func
        if not self.nb_func:
            self.nb_func = distance_cutoff

    @property
    def bonded(self):
        if not self._bonded:
            self._bonded = {(i.index, j.index) for i, j in self._top.bonds
                            if i.index in self._ind and j.index in self._ind}

        return self._bonded

    @property
    def non_bonded(self):
        if not self._non_bonded:
            combos = combinations_with_replacement(self._ind, 2)
            self._non_bonded = {(i, j) for i, j in combos if i != j}
            self._non_bonded -= self.bonded
            self._non_bonded = self.nb_func(self._structure, self._non_bonded,
                                            **self.nb_kwargs)
        return self._non_bonded

    def _solve(self):
        """Solve the eigenvalue problem for the Hessian"""
        vals, vecs = linalg.eigh(self.hessian_, turbo=self.turbo,
                                 eigvals=(0, 6 + self.mode))

        ind = np.argsort(vals)[6:]
        self.eigenvalues_ = vals[ind]
        self.eigenvectors_ = vecs[:, ind]

        self.vars_ = 1 / self.eigenvalues_
        self.trace_ = self.vars_.sum()

    def _rigid_mode(self, mode):
        """Compute an eigenmode for the system given a rigid-body constraint"""
        # Fill mode with available values from the selected eigenvector of
        # the Hessian
        arr = np.zeros((self._top.n_atoms, 3))
        arr[self._ind, :] = np.reshape(self.eigenvectors_[:, mode],
                                       (self.n_atoms_, 3))

        # Initialize graph of bonded atoms
        G = nx.Graph([(i.index, j.index) for i, j in self._top.bonds])

        # For each non-selected atom search the graph for the nearest bonded
        # atom that has been selected and have the query atom inherit its
        # gradient
        for i in (set(range(self._top.n_atoms)) - set(self._ind)):
            path_lengths = nx.single_source_shortest_path_length(G, i)
            d = np.array(list(path_lengths.values()))[self._ind]
            arr[i, :] = arr[np.argmin(d), :]

        return arr

    def _grad(self, mode):
        """Compute the gradient of a normal mode for the system"""
        step = self.rmsd / self.n_steps
        scale = step * self._top.n_atoms ** 0.5

        arr = self._rigid_mode(mode)
        grad = (arr * scale) / np.sqrt((arr**2).sum())

        return grad

    def _animate(self, mode):
        """Animate a normal mode as an mdtraj.Trajectory"""
        mid = int(self.n_steps / 2.)
        xyz = np.zeros((self.n_steps + 1, self._top.n_atoms, 3))
        xyz[mid] = self._xyz[0]

        grad = self._grad(mode)

        for i in range(mid):
            xyz[i + mid + 1] = xyz[i + mid] + grad
            xyz[mid - i - 1] = xyz[mid - i] - grad

        return md.Trajectory(xyz[:-1], self._top)

    def _set_hessian(self, m, n, g):
        """Set the value of the Hessian for a given pair of atoms (m, n)"""
        i, j = self._ind.index(m), self._ind.index(n)
        i2j = self._xyz[0, m, :] - self._xyz[0, n, :]
        dist2 = np.dot(i2j, i2j)
        super_el = np.outer(i2j, i2j) * (g / dist2)
        res_i3 = i * 3
        res_i33 = res_i3 + 3
        res_j3 = j * 3
        res_j33 = res_j3 + 3
        self.hessian_[res_i3:res_i33, res_j3:res_j33] = super_el
        self.hessian_[res_j3:res_j33, res_i3:res_i33] = super_el
        self.hessian_[res_i3:res_i33, res_i3:res_i33] = \
            self.hessian_[res_i3:res_i33, res_i3:res_i33] - super_el
        self.hessian_[res_j3:res_j33, res_j3:res_j33] = \
            self.hessian_[res_j3:res_j33, res_j3:res_j33] - super_el

    def fit(self, X):
        """Fit the model with X.
        Parameters
        ----------
        X: mdtraj.Trajectory
            A structure with which to perform ANM analysis. By default, only
            the first frame is considered.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._ind = list(X.top.select(self.selection))
        self.n_atoms_ = len(self._ind)
        self._structure = X[0]
        self._top = self._structure.top
        self._xyz = self._structure.xyz

        # Initialize Hessian Matrix
        self.hessian_ = np.zeros((3 * self.n_atoms_,
                                  3 * self.n_atoms_))

        # Add Non-Bonded Interactions to Matrix
        for i, j in self.non_bonded:
            g = - self.k_nb
            self._set_hessian(i, j, g)

        # Add Bonded Interactions to Matrix
        for i, j in self.bonded:
            g = - self.k_b
            self._set_hessian(i, j, g)

        # Solve Eigenvalue Problem
        self._solve()
        self._dirty = False

        return self

    def transform(self, X):
        """Create an animation with X.
        Parameters
        ----------
        X: mdtraj.Trajectory
            A structure with which to perform ANM analysis. By default, only
            the first frame is considered.

        Returns
        -------
        animation : mdtraj.Trajectory (n_steps, )
        """
        if self._dirty:
            self.fit(X)
            self._dirty = False

        return self._animate(self.mode)

    def fit_transform(self, X):
        """Fit the model with X and create an animation with X.
        Parameters
        ----------
        X: mdtraj.Trajectory
            A structure with which to perform ANM analysis. By default, only
            the first frame is considered.

        Returns
        -------
        animation : mdtraj.Trajectory (n_steps, )
        """
        self._dirty = True
        return self.transform(X)
