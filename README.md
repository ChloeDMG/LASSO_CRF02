# Using LASSO to predict HIV-1 CRF02-AG subtype phenotropism on *env* amino acids sequence

To predict the phenotropism of a new HIV-1 CRF02-AG subtype, a valid amino acid sequence of *env* gene must be provided.

The phenotropism prediction is done in 2 steps:
- a pre-pocessing step (python 3 script)
- the LASSO processing (MatLab script)

## Pre-processing the amino acids sequences

In this step a fasta file containing the amino acids sequences of *env* HIV-1 CRF02-AG subtypes are provided as input to the Python 3 script, `preprocessing_lasso.py`. The script requires the [BioPython](https://biopython.org/) library.

Usage is:

```bash
preprocessing_lasso.py [-h] -o OUT input

positional arguments:
  input              path to the fasta file. The file contains unaligned amino
                     acid sequences of HIV-1 env of which phenotropism have to
                     be predicted.

optional arguments:
  -h, --help         show the help message and exit
  -o OUT, --out OUT  path to the output directory.

```
The output directory will contain subdirectories, one by sequence provided in the input fasta file. In each subdirectory (<fasta_seq_id>), 2 CSV files which are the inputs for the MatLab script:
 - <fasta_seq_id>_to_predict.csv
 - <fasta_seq_id>_to_train.csv
 
 ## LASSO phenotropism prediction
