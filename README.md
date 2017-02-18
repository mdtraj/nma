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
from gnma import GNMA

# Load structure of choice (e.g. GFP)
pdb = mdtraj.load_pdb('https://files.rcsb.org/download/1EMA.pdb.gz')

# Initialize GNMA object
gnma = GNMA(mode=5, nb_cutoff=1., selection='backbone')

# Transform the PDB into a short trajectory of a given mode
gnma_traj = gnma.fit_transform(pdb)
```

Complaints
----------

[Post to the issue tracker.](https://github.com/cxhernandez/gnma/issues)

Shout-outs
----------

+ https://github.com/prody/ProDy (for loose inspiration)
