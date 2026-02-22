# BindingDB Tools

This directory contains tools and scripts for working with the BindingDB database.

## Data Source

The full BindingDB database can be downloaded from the [BindingDB Download page](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp).

## Contents

### `2026_02_24_BindingDB_Webinar.pdf`

Slides from a webinar titled "Building Machine Learning Models Using Data from BindingDB" presented February 24, 2026.

### `binding_db_pxr_tutorial.ipynb`

A Jupyter notebook demonstrating how to work with BindingDB data. This notebook can be run on Google Colab without installing any software locally: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/PatWalters/practical_cheminformatics_tutorials/blob/main/bindingdb/binding_db_pxr_tutorial.ipynb)

### `fix_tsv.py`

A script for normalizing BindingDB TSV files and loading them into a DuckDB database. BindingDB TSV files can sometimes have inconsistent column counts across rows, which can cause issues with standard CSV parsers. This script ensures every row matches the header's column count and then imports the data into DuckDB.

#### Requirements

- `duckdb`
- `tqdm`

#### Usage

```bash
python fix_tsv.py <input_tsv> <output_duckdb>
```

Options:
- `--table`: Specify the table name in DuckDB (default: `bindingdb`).
- `--tempfile`: Path for the temporary fixed TSV file.
- `--keep-temp`: Do not delete the temporary file after import.
