import matplotlib.pyplot as plt

def main():
    runs = ["Average", "No Al", "AlK(SO4)2", "Al(OH)3", "AlCl3", "LM_C"]
    points = [-0.000861, -0.00168, 4E-06, -0.000158, -0.001216, -0.001939]
    neg_xerr = [0.000121, 0.00024, 0.000331, 0.0002, 0.000314, 0.000377]
    pos_xerr = [0.000121, 0.000239, 0.000331, 0.0002, 0.000314, 0.000377]
    xerr = [neg_xerr, pos_xerr]
    plt.figure(figsize=(10,5))
    plt.scatter(points, runs)
    plt.errorbar(points, runs, xerr=xerr, fmt='none', capsize=5)
    plt.axvline(0, linestyle="--")
    plt.xlabel("Minimum standard deviation mass deltas")
    plt.ylabel("Run name")
    plt.title("Minimum standard deviation mass deltas and locations of 10% increase threshold")
    plt.xlim(-0.0030, 0.0010)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()