MS2VariantFinder generates a table of all MS2 scans in an mzML file along with relevant data such as retention time, ion injection time, total ion current, etc. It also predictively generates a USI based on fragment ion information.

You can run two scripts in this directory:
- `generate_ms2_table.py`: This is the main MS2VariantFinder script.
- `generate_stdev_plot.py`: This is a supplement script for constraining delta masses of specific modifications.

`generate_ms2_table.py` takes in the following command-line arguments:

| Argument | Description | Default | Example |
| -------- | ----------- | ------- | ------- |
| `--mzml_file` | Filepath of mzML file for analysis. | Required | "C:\data\ms_run.mzML" |
| `--sequence` | Sequence of main peptide. | Required | "AQDSQVLEEER\[Label:13C(6)15N(4)]" |
| `--charge` | Most common charge of peptide. This is required to find the most intense MS1 scan with a peak corresponding to the peptide's precursor m/z. | Required | 2 |
| `--ms2_fragment_tolerance` | Tolerance used to search for MS2 fragment ions. This is used during localization of the peptide. | 20 | 10 |
| `--ms1_precursor_tolerance` | Tolerance used to search for modifications to the peptide and precursor intensities. | 10 | 5 |
| `--run_type` | Type of MS run (DDA or PRM). | Required | "DDA" |
| `--output` | Output file name. The program will output two files: "\[output].csv" and "\[output]\_intensities.csv". | Required | "output.csv" OR "output" |

`generate_ms2_table.py` will output two files: "\[output].csv" and "\[output]\_intensities.csv".
"\[output].csv" is the main output file that summarizes all MS2 spectra in the mzML file. It comprises the following columns of information:
- "scan number": The scan number of the specific MS2 spectrum.
- "retention time": The chromatographic retention time in minutes of the specific MS2 scan.
- "ion injection time": The time in which ions are collected in the ion trap for a specific MS2 scan.
- "total ion current": The total ion current of all peaks in the MS2 spectrum.
- "precursor m/z": The mass-to-charge ratio of the identified precursor ion.
- "precursor charge": The identified charge of the precursor ion.
- "maximum precursor intensity": The intensity of the peak with the closest m/z to the precursor from within the tolerance window of the precursor m/z is calculated for each of the +/- 10 MS1 spectra from the precursor MS1 scan. The maximum of those values is found and returned in this column. Note that this is not equivalent to the intensity calculated via XIC extraction which is a summation of all peak intensities within the m/z tolerance window.
- "relative intensity": The intensity relative to the MS1 scan that has the greatest magnitude base peak within the tolerance window of the main peptide m/z across all MS1 scans in the MS run. Additionally, it should be noted that the first row of this file represents the most intense MS1 spectrum in question.
- "signal to noise ratio": This ratio is calculated crudely for each MS2 spectrum by taking the average of the second and third most intense peaks and dividing it by the average of the least and the second-to-least intense peaks in the spectrum.
- "precursor mass delta": The observed difference in mass between the precursor ion and the main peptide ion.
- "predicted mass delta": The mass delta resulting from the program's predicted modification, if any.
- "mass delta difference": Calculated from observed precursor mass delta - predicted mass delta. 
- "modification": The Unimod name of the predicted modification.
- "modification type": Cation, synthesis error, or blank if neither.
- "localization scores": This code string that logs the scores produced by different localizations to the nearest hundredth. For single-residue modifications, each residue is paired with a score that represents the possibility of the modification localizing there. The possibility of a labile modification is also indicated. The residue abbreviation is capitalized if Unimod lists it as a possible site and decapitalized if otherwise. Because the score is computed through analysis of ion fragmentation, there is no material difference in the scores of N-term and first residue modifications. For example: "Lab-4.91, N-term-5.43, a-5.43, Q-3.25, d-5.97, Q-5.74, r-2.01". For synthesis errors, each possible modified sequence is recorded. For example, if there is a missing Q: "A_DQR-2.86, AQD_R-4.93". The configuration with the highest score is chosen for the final predicted USI.
- "usi": This is the USI representing the most probable modified sequence and charge.
- "confidence": For automatically generated outputs, confidence defaults to "predicted".
- "comments": For manual annotation purposes.

"\[output]\_intensities.csv" is the output file that logs the MS1 intensities.
- "scan number": Scan number of the MS2 spectrum.
- "precursor scan number": Scan number of the precursor MS1 spectrum.
- "precursor m/z": Identified precursor ion m/z.
- "total ion current": Total ion current of the MS2 spectrum.
- "mz\[n]" (where -10<n<10): The closest m/z to the identified precursor m/z within the tolerance window of the MS1 scan n scans away from the precursor MS1 scan.
- "intensity\[n]" (where -10<n<10): The intensity of the peak identified with the closest m/z.
- "total ion intensity": The sum of all the identified intensities +/- 10 MS1 scans away from the precursor MS1 scan.
- "maximum precursor intensity": The maximum of the set of all the identified intensities +/- 10 MS1 scans away from the precursor MS1 scan.

`generate_stdev_plot.py` takes in the following command-line arguments:

| Argument | Description | Default | Example |
| -------- | ----------- | ------- | ------- |
| `--mzml_file` | Filepath of mzML file for analysis. | Required | "C:\data\ms_run.mzML" |
| `--sequence` | Sequence of main peptide. | Required | "AQDSQVLEEER\[Label:13C(6)15N(4)]" |
| `--modification` | Unimod-identified name of the investigated modification. | Required | "Cation:Al\[III]" |
| `--mod_index` | Index of modification on peptide (0-based). If on N-term, use -1. | Required | 2 |
| `--scan_number` | Scan number of MS2 spectrum. | Required | 1484 |
| `--tolerance` | Tolerance in ppm to identify fragment ions. | 6 | 5 |
| `--run_type` | Run type (DDA or PRM). | Required | "DDA" |
| `--output` | Output file name. | Required | "output.csv" or "output" |

`generate_stdev_plot.py` will output a .png image as a graph with the standard deviation score on the y-axis and the deviance from the expected mass delta on the x-axis. Additionally, the 10% increase threshold is plotted as a dotted horizontal line. The modification mass with the lowest score is plotted as a labelled star at the vertex of the graph, and the bounds of the modification uncertainty are marked at the intersection between the horizontal line threshold and the V-shaped curve.

