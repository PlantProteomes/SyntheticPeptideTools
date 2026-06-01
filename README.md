# SyntheticPeptideTools
This repository contains code for analyzing MS1 and MS2 mass spectrometry data from mzML files to support peptide & proteoform analysis in synthetic peptides.

# Usage of the main SyntheticPeptideTools programs
## MS1XICExtractor.py
The MS1XICExtractor tool extracts signal of various precursor ions. The output is a multi-page pdf showing XICs, mass accuracy deviations, total ion current (TIC), and ion injection time trends. 
## MS1XICExtractor Command-line Arguments 

| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--mzml_file` | Path to input mzML file (mass spectrometry raw data) | Required | `"C:\data\sample.mzML"` |
| `--output_file` | Path to output PDF file containing plots | Required | `output.pdf` |
| `--modifications` | Comma-separated list of target analytes in format `Name:mz` | Required | `"TargetPeptide:657.314, Aluminum:669.293"` |
| `--scan_range` | One or more scan ranges in format `start,end;start,end` | None | `"1420,1520;1600,1700"` |
| `--xic_ppm` | Mass tolerance (ppm) for XIC extraction | `4.0` | `5` |
| `--delta_ppm` | Mass tolerance (ppm) for mass delta calculation | `4.0` | `15` |
| `--max_y` | Maximum y-axis value for plots | None | `100000` |
| `--no_normalize` | Disables normalization by ion injection time | False | `--no_normalize` |
