# Path to your input files
background_path = 'Background.Spe'
cobalt_path = 'cobalt_90_NAI_TI.Spe'
output_path = 'cobalt_minus_background.Spe'

# Function to find the line index where the actual data starts (after $DATA:)
def find_data_start(content, identifier='$DATA:'):
    for idx, line in enumerate(content):
        if line.strip().startswith(identifier):
            return idx + 1  # Data typically starts right after identifier
    return None

# Function to read and extract spectrum data as a list of integers
def extract_spectrum_data(content, start_line):
    return [int(line.strip()) for line in content[start_line:] if line.strip().isdigit()]

# Load and read the content of each file
with open(background_path, 'r') as file:
    background_content = file.readlines()

with open(cobalt_path, 'r') as file:
    cobalt_content = file.readlines()

# Locate data start lines for both files
background_data_start = find_data_start(background_content)
cobalt_data_start = find_data_start(cobalt_content)

# Extract numeric data for both spectra
background_data = extract_spectrum_data(background_content, background_data_start)
cobalt_data = extract_spectrum_data(cobalt_content, cobalt_data_start)

# Ensure both spectra have the same length
if len(background_data) != len(cobalt_data):
    raise ValueError("Error: Spectrum data lengths do not match!")

# Perform background subtraction (ensure no negative values)
subtracted_spectrum = [max(cobalt - background, 0) for cobalt, background in zip(cobalt_data, background_data)]

# Prepare output with metadata from cobalt file, then append subtracted data
output_content = cobalt_content[:cobalt_data_start]  # Copy metadata header
output_content.append('$DATA:\n')  # Data section header
output_content += [f"{count}\n" for count in subtracted_spectrum]  # Append subtracted counts

# Write the result to a new file
with open(output_path, 'w') as file:
    file.writelines(output_content)

print(f"Subtracted spectrum saved to {output_path}")
