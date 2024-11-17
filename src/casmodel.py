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

    # Start assigning n = 1 from the 24th nucleotide
    n = position - 23 + 1  # Index n starting at 1 for the 24th nucleotide
    prob = calculate_dNTP_probability(base, dNTP_Conc, NRTI_Conc[0], Kaff, n)
    results.append({
        "Position": len(DNAstrand) - position,  # Reverse position
        "Base": base,
        "n": n,
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
