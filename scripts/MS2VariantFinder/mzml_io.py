from pyteomics import mzml
import re
from models import Scan, MSRun

def read_mzml(filepath, run_type):
    stats = {'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0}
    scans = []
    with mzml.read(filepath) as reader:
        for spectrum in reader:
            stats['counter'] += 1
            if spectrum["ms level"] == 1:
                stats['ms1spectra'] += 1
            if spectrum["ms level"] == 2:
                stats['ms2spectra'] += 1

            scan_number = 1 + int(spectrum['index'])
            filter_string = spectrum['scanList']['scan'][0]['filter string']
            match = re.search(r'NSI (\S+) (\S+).* \[([\d\.]+)\-([\d\.]+)\]', filter_string)

            scan_type = match.group(1) if match else None
            ms_level = int(spectrum['ms level'])
            start = float(match.group(3)) if match else None
            end = float(match.group(4)) if match else None
            isolation_window = (start, end)

            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']
            rt = float(spectrum['scanList']['scan'][0]['scan start time'])
            iit = float(spectrum['scanList']['scan'][0]['ion injection time'])
            tic = float(spectrum['total ion current'])

            precursor_mz = None
            precursor_charge = None
            last_ms1_scan = None
            if ms_level == 2:
                precursor_mz = float(
                    spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                precursor_charge = int(
                    spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                last_ms1_scan = int(
                    re.search(r"scan=(\d+)", spectrum['precursorList']['precursor'][0]['spectrumRef']).group(1))

            scans.append(Scan(
                scan_number=scan_number,
                scan_type=scan_type,
                ms_level=ms_level,
                isolation_window=isolation_window,
                mz_array=mz_array,
                intensity_array=intensity_array,
                rt=rt,
                iit=iit,
                tic=tic,
                precursor_mz=precursor_mz,
                precursor_charge=precursor_charge,
                last_ms1_scan=last_ms1_scan
            ))

    print(f"""Read mzML file {filepath}. 
    Total number of spectra: {stats['counter']}. 
    Number of MS1 spectra: {stats['ms1spectra']}. 
    Number of MS2 spectra: {stats['ms2spectra']}.""")
    if run_type == 'DDA':
        run = MSRun(scans, 'DDA')
        return run
    if run_type == 'PRM':
        run = MSRun(scans, 'PRM')
        return run
    return None


if __name__ == "__main__":
    run = read_mzml("C:\\Users\\carol\\OneDrive - Bellevue School District\\2025-2026\\6-7 Internship\\python_testing\\new peptides\\250402_mEclipse_QC_ncORF-181.mzML", run_type='DDA')
    for scan in run.scans:
        if scan.ms_level == 2:
            print(scan.scan_number, run.get_precursor(scan).scan_number)