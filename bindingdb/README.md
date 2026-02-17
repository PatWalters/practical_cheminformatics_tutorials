# BindingDB Tools

This directory contains tools and scripts for working with the BindingDB database.

## Data Source

The full BindingDB database can be downloaded from the [BindingDB Download page](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp).

## Contents

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
