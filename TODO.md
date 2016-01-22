# Map Comp TODO

## Bugs

## Features

- Tweek mapping params to get higher proportion of mapping markers

- Add utility script to rename LGs
- Add folder for output files
- Rename scripts
- Script to test mapping proportion
  - Script `01_bwa_w_fasta`

## Tutorial
- Provide fake genome
  - Copy it to `/02_raw_data/genome/genome.fasta`
  - Index it with `bwa index /02_raw_data/genome/genome.fasta'
- Provide CSV file with markers from 2 species without totpos
- Prepare the data (add totpos, create fasta)
- Run pipeline

## When paper accepted

- Add Citing info to README.md

## Maybe

- Include reciprocal best hits

- TotPos
  - Correct LGs that do not start at 0 (actually useful)

- Re-order LGs
- Compute metrics
