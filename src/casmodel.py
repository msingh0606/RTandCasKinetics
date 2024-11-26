import numpy as np 
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
from probability_func import calculate_dNTP_probability

# INPUT REAGENT CONCENTRATIONS
forwardDNA = 'TTTTTTTTTTTTTGATGATGTGAAGGTGTTGTCGTTTATTTATTTATTTATTTATTTCTATCTTTCCTCTTAATTCGACG' # the DNA template sequence
DNAstrand = forwardDNA[::-1] # reverse the DNA sequence
Primer = 'CTATCTTTCCTCTTAATTCGACG' # the primer sequence
casregion = 'ATGATGTGAAGGTGTTGTCG' #Cas12a binding region
Analog = 'T' # the drug analog
dNTPConc_nM = 500 # 500 nM; the dNTP concentration
dNTP_Conc = dNTPConc_nM * 1e-9 # convert nM to M
NRTI_Conc = np.logspace(-11,-5,1000)
templateConc = 5e-9 # 5 nM; the DNA template concentration
#Logspace = Inhibitor_Conc(2)-Inhibitor_Conc(1);

# INPUT RELATIVE INHIBITOR AFFINITY
Kaff = 0.3; #this is the inhibitory affinity factor of the drug analog relative to the dNTPs

# Iterate over the DNA sequence and calculate probabilities
results = []
found_casregion = False

# Start processing from the 24th nucleotide
for position, base in enumerate(DNAstrand):
    # Skip the first 23 nucleotides
    if position < 23:
        continue

    # Index for result tracking (starts at 0 after skipping the first 23 nucleotides)
    result_index = position - 23  

    # Start assigning n = 1 for the function calculation
    function_n = result_index + 1  

    # Calculate the probability
    prob = calculate_dNTP_probability(base, dNTP_Conc, NRTI_Conc[0], Kaff, function_n)

    # Append the result
    results.append({
        "Position": len(DNAstrand) - position,  # Reverse position
        "Base": base,
        "n": result_index,  # Index used for results
        "Function_n": function_n,  # n used for function calculation
        "P_dNTP,n": prob
    })

    # Check for the target sequence (casregion)
    if not found_casregion:
        if DNAstrand[position:position + len(casregion)] == casregion[::-1]:
            found_casregion = True

    # Stop after processing the first 7 nucleotides of casregion
    if found_casregion and DNAstrand[position:position + 7] == casregion[::-1][:7]:
        break  # Exit the loop once the first 7 nucleotides are processed

# Display results
results_df = pd.DataFrame(results)
print(results_df)

import tellurium as te
import matplotlib.pyplot as plt

# Define the model in Antimony format with literature values
model = """
model MichaelisMentenFluorescence
    // Parameters (example literature values)
    kcat = 0.55;      // Turnover number (1/s)
    Km = 663e-9;         // Michaelis constant (M)
    E = 25e-9;       // Enzyme concentration (M)
    S = 500e-9;          // Initial reporter concentration (M)
    P = 0.0;          // Initial fluorescence (product concentration)

    // Reaction (substrate to product)
    v: S -> P; (kcat * E * S) / (Km + S);

    // Species initial concentrations
    S = 1.0; // Initial substrate concentration
    P = 0.0; // Initial fluorescence
end
"""

# Load the model
rr = te.loadAntimonyModel(model)

# Simulate over time
result = rr.simulate(0, 1000000, 2000)  # Simulate from t=0 to t=10 with 100 points

# Extract time and fluorescence (product concentration)
time = result[:, 0]  # Time column
fluorescence = result[:, 2]  # Product (P) column

# Plot fluorescence vs. time
plt.figure(figsize=(8, 5))
plt.plot(time, fluorescence, label="Fluorescence", color="blue")
plt.title("Cas12a Michaelis-Menten Kinetics: Fluorescence vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Fluorescence (RFU)")
plt.legend()
plt.grid(True)
plt.show()