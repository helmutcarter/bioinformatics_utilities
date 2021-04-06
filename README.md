# Helmut A Carter's Bioinformatics Utilities and Portfolio

## motif_counter.py
To use: edit the variables at the top with your motif of interest and query sequence. Requires no packages.


## rMATS2HOMER.py
Purpose: Accepts [rMATS](https://github.com/Xinglab/rmats-turbo) output and converts into sequence files that can be used by [HOMER](http://homer.ucsd.edu/homer/motif/), a motif enrichment discovery program

Dependencies:
- [pyfaidx](https://pypi.org/project/pyfaidx/) `pip install pyfaidx`
- [Biopython](https://biopython.org/) `pip install biopython`
