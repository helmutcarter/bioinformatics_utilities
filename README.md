# Helmut A Carter's Bioinformatics Utilities Portfolio

## motif_counter.py
To use: edit the variables at the top with your motif of interest and query sequence.

Dependencies:
- None.

## rMATS2HOMER.py
Purpose: Accepts [rMATS](https://github.com/Xinglab/rmats-turbo) output and converts into sequence files that can be used by [HOMER](http://homer.ucsd.edu/homer/motif/), a motif enrichment discovery program.

Dependencies:
- [pyfaidx](https://pypi.org/project/pyfaidx/) `pip install pyfaidx`
- [Biopython](https://biopython.org/) `pip install biopython`

## rMATS_compare.py
Purpose: Accepts [rMATS](https://github.com/Xinglab/rmats-turbo) output for two different species and compares alternative splicing events between the two species to predict which events are orthologous.

Dependencies:
- [pyfaidx](https://pypi.org/project/pyfaidx/) `pip install pyfaidx`
- [Biopython](https://biopython.org/) `pip install biopython`
