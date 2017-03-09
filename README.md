GNMA
====

Install
-------

```bash
$ git clone git@github.com:cxhernandez/gnma.git && cd gnma
$ python setup.py Install
```

Usage
-----

```python
# Imports
import mdtraj as md
from nma import ANMA

# Load structure of choice (e.g. GFP)
pdb = md.load_pdb('https://files.rcsb.org/download/1GFL.pdb.gz')

<<<<<<< HEAD
# Initialize GNMA object
anma = ANMA(mode=0, nb_cutoff=1., k=1., n_steps=50, selection='chainid 0 and backbone')
=======
# Initialize ANMA object
anma = ANMA(mode=0, rmsd=0.06, n_steps=50, selection='all')
>>>>>>> 73c0345... hessian

# Transform the PDB into a short trajectory of a given mode
anma_traj = anma.fit_transform(pdb)
```

![](https://raw.githubusercontent.com/cxhernandez/gnma/master/examples/gfp.gif)

Complaints
----------

[Post to the issue tracker.](https://github.com/cxhernandez/gnma/issues)

Shout-outs
----------

+ [ProDy](https://github.com/prody/ProDy) (for loose inspiration)
