# SyntheticPeptideTools
This repository contains code for analyzing MS1 and MS2 mass spectrometry data from mzML files to support peptide & proteoform analysis in synthetic peptides.

# Main SyntheticPeptideTools programs
## MS2VariantFinder
MS2VariantFinder attempts to identify all MS2 spectra in an MS run based on the assumption that they are a variant of one synthesized target peptide. The variations can be any mass modification known in Unimod or extra/missing residues from the target peptide.
For more information, please consult [MS2VariantFinder.md](https://github.com/PlantProteomes/SyntheticPeptideTools/blob/main/scripts/MS2VariantFinder/MS2VariantFinder.md) in the MS2VariantFinder directory within scripts.
## MS1XICExtractor
MS1XICExtractor creates extracted ion chromatograms (XICs) to analyze the abundance of various precursor ions and, based on either MS1 full range scans or SIM scans, generates corresponding XIC plots, mass accuracy deviations, total ion current, and ion injection time trends.
For more information, please consult [MS1XICExtractor.md](https://github.com/PlantProteomes/SyntheticPeptideTools/blob/main/scripts/MS1XICExtractor.md) within scripts.
