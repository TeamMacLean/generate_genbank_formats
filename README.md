## Introduction

The script generated Genbank format files from tab separated data files. Data are input from specific columns. Therefore, the input format and the column number should never be changed for the script to work. However, the command to generate genbank file can be re-used.

## Requirements

1) python (v2.7+ or v3.5+)
2) biopython

## Usage:

python3 scripts/generate_genbank_files.py tab-delimited-file.csv       # script generates genbank format files
python3 scritps/get_genbank_record_by_id.py list-of-genban-ids.csv     # script gets geneid for proteins using genbank protein ids
