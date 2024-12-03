import numpy as np

def calculate_annealing_time(P0, D, k, target_percentage):
    """
    Calculate the time for a primer to anneal to a DNA template.
    
    Parameters:
    - P0 (float): Initial primer concentration (M)
    - D (float): DNA concentration (M)
    - k (float): Reaction rate constant (M^-1 s^-1)
    - target_percentage (float): Target percentage of annealing (0-100)
    
    Returns:
    - t_anneal (float): Time (seconds) to reach the target annealing percentage
    - time (numpy.ndarray): Time array for plotting (seconds)
    - P (numpy.ndarray): Primer concentration over time (M)
    """
    # Calculate the target primer concentration
    P_target = P0 * (1 - target_percentage / 100)

    # Calculate the annealing time
    t_anneal = -np.log(P_target / P0) / (k * D)

    # Generate time points for plotting
    time = np.linspace(0, t_anneal * 1.5, 1000)  # Simulate up to 1.5 times the annealing time
    P = P0 * np.exp(-k * D * time)  # Primer concentration over time

    return t_anneal, time, P
