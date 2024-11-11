import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import inf as INF
#from datetime import datetime


#Defining Constant Variable
TWO_PI = np.pi*2
ROI = [0, 1000]
xlim = [0, 800]
Gaussian_params = ('mu', 'sig', 'a', 'c0', 'c1', 'c2')
Gaussian_bounds_background = [

    np.array([0, 0, 0, -INF, -INF, -INF]),
    np.array([INF, INF, INF, INF, INF, INF])
]

####################################################### Reading Raw Spectrum ##########################################################


def read_raw_spectrum(filenames):
    # Initialize an empty dictionary to store the spectra data
    spectrum_dict = {}

    for filename in filenames:
        with open(filename, 'r') as file:
            lines = file.readlines()  # Read all lines once for easy processing
            
            # Initialize variable to track where data starts
            data_start = None
            counts = []
            
            # Check for the type of file and where the data section starts
            if any('<<DATA>>' in line for line in lines):  # MCA format
                for i, line in enumerate(lines):
                    if '<<DATA>>' in line:
                        data_start = i + 1  # Data starts right after the <<DATA>> line
                        break
                
                # Read data (assuming 1024 lines of integer counts after <<DATA>>)
                counts = [int(line.strip()) for line in lines[data_start:data_start + 1024]]
            
            elif any('$DATA:' in line for line in lines):  # SPE format
                for i, line in enumerate(lines):
                    if '$DATA:' in line:
                        data_start = i + 2  # Data starts two lines after the $DATA: line
                        break
                
                # Read data (assuming 1024 lines of integer counts)
                counts = [int(line.strip()) for line in lines[data_start:data_start + 1024]]
            
            # Convert counts to numpy array and create corresponding channels
            final_counts = np.asarray(counts)
            channels = np.arange(len(final_counts))
            
            # Store the data in the dictionary with filename as the key
            spectrum_dict[filename] = {"counts": final_counts, "channels": channels}
    
    return spectrum_dict

def plot_the_spectrum(channels, total_counts, filename, popt, pcov):
        
        colours = colourmask(channels, xmin = ROI[0], xmax = ROI[1])
       
        fig, ax = plt.subplots(figsize=(9,6))
        plt.plot(channels, total_counts, label = filename)
        plt.scatter(channels, total_counts, c = colours, marker = '+', label = 'Region of Interest(ROI)' )
        
        #Plot the curve fit with background
        popt_background = popt[3:]
        
        plot_model(ax, gaussian_plus_background, ROI, popt, linestyle = 'dashed', c='g', label='Gaussian-Continuum Fit')
        plot_model(ax, background_slope, ROI, popt_background, linestyle = 'dashed', c = 'b', label = 'Continuum Fit')
        plt.title('Counts VS Channels')
        plt.xlabel('Channels')
        plt.ylabel('Counts')
        plt.xlim(xlim[0], xlim[1])
        plt.legend()
        plt.grid(True)
        plt.show()
        
        
def desired_selected_data(spectrum_dict, filename) :
    
    if filename in spectrum_dict:
        data         = spectrum_dict[filename]
        total_counts = data["counts"]
        channels     = data["channels"]
        return channels, total_counts
    
    else:
        print("File cannot be found in the dictionary")


def gaussian_first_fit(x, mu, sig, a):
    
    return a * np.exp(-0.5 * (x-mu)**2 / sig**2) / np.sqrt(TWO_PI * sig**2)
  

def in_interval(x, xmin = -INF, xmax = INF):

    """
    Boolean Mask with Value True for x in xmax and xmin
    """

    _x = np.asarray(x)
    return np.logical_and(xmin <= _x, _x < xmax)
  

def filter_in_interval(x, y, xmin, xmax):

    """
    Select only elements of x and y where xmin <= _x_wavelength < xmax
    """
    # Check if X and y lists has the same shape
    # Note .shape is only attributed and used with numpy array, not lists
    # Check the type of data before checking their shape 

    x = np.asarray(x)
    y = np.asarray(y)
    
    print("\n\nChecking the dimensions of both X and Y data ......")

    if x.shape != y.shape:
    
        print("X and Y lists has no same shape, check your code\n")

    else:
		
        print("X-data and Y-data have the same shape, you are good to proceed!\n")

    _mask = in_interval(x, xmin, xmax)
            
    return [np.asarray(a)[_mask] for a in (x,y)]



def colourmask(x, xmin = -INF, xmax = INF, cin = 'red', cout = 'gray'):

    """
    This is to colour cin if within the region of interest
    """
    _mask = np.array(in_interval(x, xmin, xmax), dtype = int)

    #Convert to colours 
    colourmap = np.array([cout, cin])
    return colourmap[_mask]


def simple_model_fit(model, channels, total_counts, roi, **kwargs):

    """
    Least Square Estimate of model parameters.
    """

    #Select relevant channels and counts
    _channels, _total_counts = filter_in_interval(channels, total_counts, *roi)

    #Fit the model to the data
    popt, pcov = curve_fit(model, _channels, _total_counts, **kwargs)

    return popt, pcov



def format_result(params, popt, pcov):

    """
    Display parameter best estimates and uncertainties
    """

    #Extract the uncertainties from the covariance matrix
    perr = np.sqrt(np.diag(pcov))

    #Format parameters with best estimates and uncertainties

    _lines = (f"{p} = {o:.4f} ± {e:.4f}" for p,o,e in zip(params, popt, perr))

    return "\n".join(_lines)
   

def plot_model(ax, model, xrange, ps, npoints=1023, **kwargs):

    """
    Plots an 1d model on an Axes smoothly over xrange
    """
    _channels          = np.linspace(*xrange, npoints)
    _total_counts      = model(_channels, *ps)
    
    return ax.plot(_channels, _total_counts, **kwargs)
     

def first_moment(x, y):

    x = np.asarray(x)
    y = np.asarray(y)
    
    return np.sum(x*y) / np.sum(y)


def second_moment(x, y):

    x = np.asarray(x)
    y = np.asarray(y)
    x0 = first_moment(x, y)

    return np.sum(((x-x0)**2) *y) / np.sum(y)

 
def gaussian_initial_estimates(channels, total_counts):

    """
    Estimates of three parameters of the gaussian distribution
    """

    #Initial Guess centroid of the peak wavelength
    mu0  = first_moment(channels, total_counts)

    #Initial Estimates Standard Deviation
    sig0 = np.sqrt(second_moment(channels, total_counts))

    #Initial Estimates Amplitude
    a0 = np.max(total_counts)
    
    c0, c1, c2 = continuum_point(channels, total_counts)

    return (mu0, sig0, a0, c0, c1, c2)


# Subtracting Background Components

def background_slope(x, c0, c1, c2):

    return (c2*(x**2)) + c1*x + c0

  
def continuum_plot(channels, total_counts):

    fig, ax = plt.subplots(figsize=(9,5))
    popt, pcov = curve_fit(background_slope, channels, total_counts)
    
    plot_model(ax, background_slope, (ROI[0], ROI[1]), popt, linestyle = 'dashed', c = 'b', label = 'Continuum Fit')

    plt.legend()


def continuum_point(channels, total_counts):
    
    popt, pcov = curve_fit(background_slope, channels, total_counts)  
    c0, c1, c2 = popt[0], popt[1], popt[2]
    
    return c0, c1, c2


def gaussian_plus_background(x, mu, sig, a, c0, c1, c2):

    """
    This is the definition to illustrate a gaussian on a linear background line.
    """
    return gaussian_first_fit(x, mu, sig, a) + background_slope(x, c0, c1, c2)


def plot_curve_fit(popt, pcov):

    fig, ax = plt.subplots(figsize = (9,5))	

    #Plot the curve fit with background
    popt_background = popt[3:]
    
    plot_model(ax, gaussian_plus_background, ROI, popt, linestyle = 'dashed', c='r', label='Gaussian-Continuum Fit')
    plot_model(ax, background_slope, ROI, popt_background, linestyle = 'dashed', c = 'b', label = 'Continuum Fit')
    
    plt.legend()
    plt.title('Wavelengths VS Fluxes (Gaussian-Continuum Fit)')
    plt.show()


def optimal_value(channels, total_counts):

    c0, c1, c2 = continuum_point(channels, total_counts)

    #Make Initial Estimates or Guesses
    _channels, _total_counts = filter_in_interval(channels, total_counts, *ROI)
    _p0 = gaussian_initial_estimates(channels, total_counts)

    #Show the initial guess
    print("\n> The initial estimates:\n")
    print("\n".join(f"{p} = {o:.4f}" for p, o in zip(Gaussian_params, _p0)))

    #Do the fit 
    popt, pcov = simple_model_fit(gaussian_plus_background, channels, total_counts, ROI, p0 = [_p0[0], _p0[1], _p0[2], c0, c1, c2] , bounds = Gaussian_bounds_background)
    peak_centroid = popt[0]
    
    # Find the corresponding channel number for the centroid
    peak_channel_index = np.argmin(np.abs(channels - peak_centroid))
    peak_channel_number = channels[peak_channel_index]

    return popt, pcov, peak_centroid, peak_channel_number


##########################################################################################################################################
# Define Energy Calibration
# NAITI
peak_energy_channels	= [26, 179, 600, 337]
photopeak_energy        = [59.05, 356.0129, 1173.2, 661.657]


# Define the linear model

def linear_model(x,m,c):
     
     return m*x + c

#Perform Calibration

def calibrate_detector(centroids, known_energies):
    
    # Perform Curve Fitting for Linear Relationship between Peak Energy Level to Centroids
    popt, _ = curve_fit(linear_model, centroids, known_energies)
    m, c = popt
    print(f"Calibration result: slope (m) = {m:.4f}, intercept (c) = {c:.4f}")
    
    return m, c

def convert_channels_to_energy(channels, m, c):
    
    return m * np.array(channels) + c


def plot_channel_energy_relation(energies, total_counts, m, c, peak_energy_level):
    
    plt.subplots(figsize=(10,6))
    plt.plot(energies, total_counts, label='Calibrated Spectrum (Counts VS Photopeak (keV)')
    plt.xlabel("Photopeak Energy (keV)")
    plt.ylabel("Counts")
    plt.title("Calibrated Spectrum (Counts VS Photopeak Energy (keV))")
    plt.axvline(peak_energy_level, color='red', linestyle='--', linewidth=1.5, 
                    label=f'Peak Energy Level: {peak_energy_level:.2f} keV')
    plt.xlim(xlim[0], xlim[1])
    plt.legend()
    plt.grid(True)
    plt.show()
    
#################################################### SOurce Activity Calibration #######################################################
"""
# Conversion factor from µCi to Bq
uCi_to_Bq = 3.7e4

HALF_LIVES = {
    'Am': 432.2,
    'Cs': 30.09,
    'Ba': 10.537,
    'Co': 5.2714
}

INITIAL_ACTIVITIES_uCi = {
    'Am': 11.92,  # Replace with the actual activity for Americium source
    'Cs': 12.41,   # Replace with the actual activity for Cesium source
    'Ba': 10.82,   # Replace with the actual activity for Barium source
    'Co': 11.32    # Replace with the actual activity for Cobalt source
}

# Convert initial activities to Bq (ChatGPT Genrated Code)
INITIAL_ACTIVITIES_Bq = {source: activity * uCi_to_Bq for source, activity in INITIAL_ACTIVITIES_uCi.items()}

# Define the decay constant function
def decay_constant(half_life):
    return np.log(2) / half_life

# Modified function to calculate activity on a specific date
def calculate_activity(initial_activity_Bq, half_life, start_date, current_date):
    
    # Calculate elapsed time in years between start_date and current_date
    elapsed_years = (current_date - start_date).days / 365.25
    
    # Calculate the decay constant
    lambda_val = decay_constant(half_life)
    
    # Calculate the activity on the specified date
    return initial_activity_Bq * np.exp(-lambda_val * elapsed_years)

# Calibration dates
start_dates = {
    'Am': datetime(1979, 2, 1),
    'Cs': datetime(1979, 2, 1),
    'Ba': datetime(1979, 2, 1),
    'Co': datetime(1979, 2, 1)
}

# Define the specific current date for which you want to calculate the activity
specific_date = datetime(2024, 9, 10)  

# Calculate and print the current activity for each source as of the specific date (ChatGPT Generated)
for source, half_life in HALF_LIVES.items():
    activity = calculate_activity(INITIAL_ACTIVITIES_Bq[source], half_life, start_dates[source], specific_date)/uCi_to_Bq
    print(f"\nActivity for {source} as of {specific_date.date()}: {activity:.2f} uCi\n")    
"""
############################################################## Calculate Resolution ####################################################

def calculate_resolution(popt, peak_channel_correspond_energy_level):
    
    fwhm = 2.3548* popt[1]
    energy = peak_channel_correspond_energy_level
    
    
    resolution = fwhm/energy
    
    return resolution

def resolution_model(E, a, b, c):
        
    return a * E**-2 + b * E**-1 + c

# Resolution
NAI_TI_resolution = [0.22, 0.1, 0.04, 0.05] #Clearly Identify Peak Am, Ba, Co, Cs
NAI_TI_energy_level = [59.05, 356.0129, 1173.2, 661.657] #Correspond Collected Resolution

energy_levels = np.array(NAI_TI_energy_level)
energy_resolution = np.array(NAI_TI_resolution)


def r_squared(energy_resolution):
    
    R_squared = energy_resolution**2

    return R_squared    
    
def plot_resolution_model():
    
    R_squared = r_squared(energy_resolution)
    popt, pcov = curve_fit(resolution_model, energy_levels, R_squared)
    
    a, b, c = popt
    
    #Plot the Resolution_energy 
    E_fit = np.linspace(min(energy_levels), max(energy_levels), 1000)
    R_squared_fit = resolution_model(E_fit, *popt)
    
    # Plot the resolution vs energy (log-log plot)
    plt.figure(figsize=(10, 6))
    plt.plot(energy_levels, R_squared, 'ro', label='Resolution for Different Radioactive Source BGO')
    plt.plot(E_fit, R_squared_fit, 'b-', label=f'Fit: a={a:.3f}, b={b:.3f}, c={c:.3f} , aE^-2 + bE^-1 +c')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Resolution Squared (R²)')
    plt.title('Energy Resolution vs Peak Energy')
    plt.legend()
    plt.grid(True)
    plt.show()
    
# Data Collected for Absolute Efficiency Calculation

radius_of_detector = 0.0294 # 2.94cm for BGO and NAL(TI), 0.5cm for CdTe
radius_of_emission = 0.2 #20cm for all detectors
factor_g = (1/4)*((radius_of_detector)/(radius_of_emission))
raw_counts = [32307, 990, 40, 5134] # Am, Ba, Co, Cs NAITI
count_rate = np.array(raw_counts)/300
photon_rate_uci_full = [11.08, 0.54, 0.03, 4.34] # Standard unit in uCi
emission_fraction = [0.3578, 0.6205, 0.9985, 0.8499]
photon_rate_uci = np.array(photon_rate_uci_full)* np.array(emission_fraction)
photon_rate_bp = np.array(photon_rate_uci)*(3.7e4)  
photon_rate_bp_full = np.array(photon_rate_uci_full)*(3.7e4) 

###################################################### Define the efficiency of the detector #############################################

def calculate_full_efficiency(count_rate, photon_rate_bp_full):
    
    full_efficiency = np.array(count_rate)/ np.array(photon_rate_bp_full)
    
    return full_efficiency

def calculate_intrinsic_efficiency(factor_g, full_efficiency):
    
    intrinsic_efficiency = np.array(full_efficiency)/factor_g
    
    return intrinsic_efficiency

# Plotting Efficiency to energy plot

def intrinsic_efficiency_model(energy, a, b, c):
    
    return a + b * np.log(energy) + c * (np.log(energy) ** 2)


def plot_efficiency_energy(energies, efficiencies):
    
    """Fit and plot the intrinsic efficiency model using curve_fit."""
    popt, _ = curve_fit(intrinsic_efficiency_model, energies, efficiencies)
    a, b, c = popt
    
    plt.figure(figsize=(10, 6))
    plt.plot(energies, efficiencies, 'ro', label='Efficiency for Different Radioactive Source BGO')
    plt.plot(energies, intrinsic_efficiency_model(energies, *popt), 'b-', 
             label=f'Fit: a={a:.3f}, b={b:.3f}, c={c:.3f}, a + blnE +c(lnE)^2')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Log Energy (keV)')
    plt.ylabel('Log Efficiency')
    plt.title('Intrinsic Efficiency vs Peak Energy')
    plt.legend()
    plt.grid(True)
    plt.show()        
    
    return a,b,c

def main():
    
    m, c = calibrate_detector(peak_energy_channels, photopeak_energy)
    
    filenames = ['americium_90_NAI_TI.Spe',
                 'cobalt_minus_background_NAITI.Spe',
                 'barium_90_NAI_TI.Spe',
                 'cesium_90_NAI_TI.Spe',
                 'cesium_0_CdTe.mca',
                 'cobalt_0_CdTe.mca',
                 'barium_0_CdTe.mca',
                 'americium_0_CdTe.mca',
                 'barium_90_BGO.Spe',
                 'americium_90_BGO.Spe',
                 'cobalt_90_BGO.Spe',
                 'cesium_90_BGO.Spe',
                 'cobalt_NAI_TI.Spe']
    
    spectrum_data = read_raw_spectrum(filenames)

    while True: 
        
        print("Available files to plot:")
        
        for filename in spectrum_data.keys():  
            print(f" - {filename}")
            
            
        selected_filename = input("\n Enter the filename that you want to plot: ")        
        desired_data = desired_selected_data(spectrum_data, selected_filename)
        channels, total_counts = desired_data           
        
        
        # Evaluating Optimal Value
        popt, pcov, peak_centroid, peak_channel_number = optimal_value(channels, total_counts)
        plot_the_spectrum(channels, total_counts, selected_filename, popt, pcov)
        print("\n> The final fitted estimates:\n")
        print(format_result(Gaussian_params, popt, pcov))
        
        
        # Indexing the channel number
        print(f"\n> The Centroid of the Energy Peak located at Channel:\n{peak_channel_number}")
        
        # Plotting the Energy-Channel Relation Mapping
        energies = convert_channels_to_energy(channels, m, c)
        peak_channel_correspond_energy_level = convert_channels_to_energy(peak_channel_number, m, c)
        plot_channel_energy_relation(energies, total_counts, m, c, peak_energy_level=peak_channel_correspond_energy_level)
        
        print(f"\n> The Correspond Peak Energy:\n {peak_channel_correspond_energy_level: .2f}")
        
        # Calculating Resolution
        resolution = calculate_resolution(popt, peak_channel_correspond_energy_level)
        print(f"\n> The specific resolution for this detector:\n {resolution: .2f}")
        
        # Plotting Resolution Model
        plot_resolution_model()
        
        # Calculating Efficiency
        full_efficiency = calculate_full_efficiency(count_rate, photon_rate_bp_full)
        print(f"\n The absolute efficiency of the sources (Am, Ba, Co, Cs) are: {full_efficiency} respectively")
        
        intrinsic_efficiency = calculate_intrinsic_efficiency(factor_g, full_efficiency)
        print(f"\n The intrinsic efficiency of the sources (Am, Ba, Co, Cs) are: {intrinsic_efficiency} respectively")
                
        # Plotting Intrinsic Efficiency
        energies_list = [59.16, 358.13, 662.43, 1172.84] #Am, Ba, Co, Cs
        efficiencies = [0.00714787, 0.00449429, 0.00289992, 0.00326857]
        energies = np.array(energies_list)
        plot_efficiency_energy(energies, efficiencies)
        
        break
    
    
if __name__ == "__main__":
    
    main()