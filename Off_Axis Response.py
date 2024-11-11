import matplotlib.pyplot as plt

def load_spectrum(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Locate the start of spectrum data (usually after a specific header line)
    spectrum_start = None
    for i, line in enumerate(lines):
        if line.strip().isdigit():  # usually the spectrum data is just numbers
            spectrum_start = i
            break

    # If we found the start, parse the data
    if spectrum_start is not None:
        counts = [int(line.strip()) for line in lines[spectrum_start:] if line.strip().isdigit()]
        return counts
    else:
        print(f"Data not found in file: {filename}")
        return []

# Load both spectra data using filenames (no need for absolute paths in Spyder)
counts_1 = load_spectrum('cobalt_NAI_TI.Spe')
counts_2 = load_spectrum('cobalt_30_NAI_TI.Spe')
counts_3 = load_spectrum('cobalt_60_NAI_TI.Spe')
counts_4 = load_spectrum('cobalt_90_NAI_TI.Spe')

# Plot the spectra
plt.figure(figsize=(10, 6))
plt.plot(counts_1, label='cobalt_NAI_TI.Spe')
plt.plot(counts_2, label='cobalt_30_NAI_TI.Spe')
plt.plot(counts_3, label='cobalt_60_NAI_TI.Spe')
plt.plot(counts_4, label='cobalt_90_NAI_TI.Spe')
plt.xlabel("Channel")
plt.ylabel("Counts")
plt.xlim(0,100)
plt.title("Comparison of Spectra")
plt.legend()
plt.grid(True)
plt.show()
