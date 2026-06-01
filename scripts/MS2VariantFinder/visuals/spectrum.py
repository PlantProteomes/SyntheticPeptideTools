from models import Scan
from matplotlib import pyplot as plt

def plot_ms1(spectrum: Scan):
    markers, stems, base = plt.stem(spectrum.mz_array, spectrum.intensity_array, markerfmt=" ")
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.suptitle("Intensity vs. m/z of MS1 scan number " + str(spectrum.scan_number))
    return

def plot_ms2(spectrum):
    markers, stems, base = plt.stem(spectrum.mz_array, spectrum.intensity_array, markerfmt=" ")
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.suptitle("Intensity vs. m/z of MS2 scan number " + str(spectrum.scan_number))
    return