# How to run MS1XICExtractor
## Setup
Follow instructions in [installation](https://github.com/PlantProteomes/SyntheticPeptideTools/tree/main#installation) to download the repository on your computer and install required packages. 

the repository includes the sample mzML data file in example/example.mzML.gz
## Run the sample data 
Run: 

py MS1XICExtractor.py 
  --mzml_file ..\example\example.mzML.gz 
  --output_file output.pdf 
  --modifications "TargetPeptide:657.314, Aluminum:669.293" 
  --scan_range "1420,1520" 
  --xic_ppm 4 
  --delta_ppm 15

The program will produce output.pdf, which contains analysis results. 
