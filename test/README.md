# BISCUIT Testing Framework

## Initial Framework

The initial framework for testing BISCUIT compares the outputs from the HEAD commit of the BISCUIT master branch on
GitHub to the HEAD commit of whatever branch you are working on.

### Setting up BISCUIT and a reference FASTA

To set up BISCUIT, including downloading from GitHub and compiling, run:

```
cd data/
bash get_biscuit.sh
```

Once that is completed, a reference FASTA can be retrieved via:

```
# from data/ still
bash get_fasta.sh
```

This will also create the requisite index files for a BGZIP-compressed FASTA file.
