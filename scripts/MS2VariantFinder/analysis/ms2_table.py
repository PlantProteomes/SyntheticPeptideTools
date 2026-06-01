from models import MSRun
from models import Peptide
import constants
from usi import generate_usi
from intensity import calculate_precursor_intensity, find_max_ms1
import numpy as np

def generate_ms2_table(run: MSRun, sequence: Peptide, charge: int, tolerance, run_type):         # PRM might need extra param for mod list
    expected_mass = sequence.mass(charge)
    expected_mz = sequence.mz(charge)
    spectra_data = []
    intensity_data = []
    best_ms1_spectrum, max_ms1_intensity, max_ms1_mz, best_mz_row, best_intensity_row, best_sn_ratio = find_max_ms1(run, expected_mz, 10)
    intensity_data.append((best_intensity_row, best_mz_row))
    spectra_data.append({"scan number": best_ms1_spectrum.scan_number,
                         "retention time": best_ms1_spectrum.rt,
                         "ion injection time": best_ms1_spectrum.iit,
                         "total ion current": best_ms1_spectrum.tic,
                         "precursor m/z": max_ms1_mz,
                         "precursor charge": charge,
                         "maximum precursor intensity": max_ms1_intensity,
                         "relative intensity": 1,
                         "signal to noise ratio": best_sn_ratio,
                         "modification": "",
                         "modification type": "",
                         "expected mass delta": 0,
                         "mass delta difference": max_ms1_mz - expected_mz,
                         "localization scores": "",
                         "usi": f"mzspec:PXD{999007}:{run}:{best_ms1_spectrum.scan_number}:{sequence}/{charge}",
                         "confidence": "predicted"})

    for scan in run.ms2_spectra:
        mass_delta = scan.precursor_mz * scan.precursor_charge - expected_mass - constants.PROTON_MASS * (scan.precursor_charge - charge)

        sorted_intensity = np.argsort(scan.intensity_array)
        signal = np.sum(sorted_intensity[-3:-1]) / 2
        noise = np.sum(sorted_intensity[0:2]) / 2
        sn_ratio = signal / noise

        max_precursor_intensity, total_precursor_intensity, max_precursor_mz, mz_row, intensity_row = calculate_precursor_intensity(scan.precursor_mz, run.get_precursor(scan), run, 10)
        relative_intensity = max_precursor_intensity / max_ms1_intensity

        modded_sequence, scores, mod_string, theoretical_delta, mod_type = generate_usi(scan, sequence, mass_delta, tolerance)
        if modded_sequence:
            usi = f"mzspec:PXD{999007}:{run}:{scan.scan_number}:{modded_sequence}/{scan.precursor_charge}" # predict USI
            confidence = "predicted"
        else:
            usi = ""
            confidence = ""

        spectra_data.append({"scan number": scan.scan_number,
                             "retention time": scan.rt,
                             "ion injection time": scan.iit,
                             "total ion current": scan.tic,
                             "precursor m/z": scan.precursor_mz,
                             "precursor charge": scan.precursor_charge,
                             "maximum precursor intensity": max_precursor_intensity,
                             "relative intensity": relative_intensity,
                             "signal to noise ratio": sn_ratio,
                             "modification": mod_string,
                             "modification type": mod_type,
                             "expected mass delta": theoretical_delta,
                             "mass delta difference": mass_delta - theoretical_delta,
                             "localization scores": scores,
                             "usi": usi,
                             "confidence": confidence})
        intensity_data.append((intensity_row, mz_row))

    return spectra_data, intensity_data

