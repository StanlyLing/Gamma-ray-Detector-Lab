###############################################################################################################
# Calculate Source Activity 
# Half-lives in years (data from manual Table 1)


from datetime import datetime
import numpy as np

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
specific_date = datetime(2024, 9, 23)  

# Calculate and print the current activity for each source as of the specific date (ChatGPT Generated)
for source, half_life in HALF_LIVES.items():
    activity = calculate_activity(INITIAL_ACTIVITIES_Bq[source], half_life, start_dates[source], specific_date)/uCi_to_Bq
    print(f"Activity for {source} as of {specific_date.date()}: {activity:.2f} uCi")
    
    
    
    
    
    
    
    
    
    
    
    