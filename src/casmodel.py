import numpy as np
import tellurium as te
import matplotlib.pyplot as plt

from probability_func import calculate_dNTP_probability  # Ensure this module is correctly defined

# INPUT REAGENT CONCENTRATIONS
forwardDNA = 'TTTTTTTTTTTTTGATGATGTGAAGGTGTTGTCGTTTATTTATTTATTTATTTATTTCTATCTTTCCTCTTAATTCGACG'
DNAstrand = forwardDNA[::-1]
casregion = 'ATGATGTGAAGGTGTTGTCG'  # Cas12a binding region
dNTPConc_nM = 500  # dNTP concentration in nM
dNTP_Conc = dNTPConc_nM * 1e-9  # Convert to M
NRTI_Conc = np.logspace(-10, -5, 5)  # Example: 5 points for simplicity
Kaff = 0.3  # Relative inhibitory affinity

# Define nucleotide addition rate
nucleotide_addition_rate = 18  # Basepairs per second

# Initialize plot
plt.figure(figsize=(10, 6))

# Simulate for each NRTI concentration
for nrti_conc in NRTI_Conc:
    print(f"Simulating for NRTI_Conc = {nrti_conc:.1e} M")
    cumulative_time = 0  # Reset cumulative time for each concentration
    found_casregion = False

    # RT synthesis simulation
    for position, base in enumerate(DNAstrand):
        if position < 23:
            # Skip the first 23 nucleotides
            cumulative_time += 1 / nucleotide_addition_rate
            continue

        # Calculate probability with corrected argument for 'n'
        prob = calculate_dNTP_probability(base, dNTP_Conc, nrti_conc, Kaff, n=position - 23 + 1)

        time_delay = (1 / nucleotide_addition_rate) / prob if prob > 0 else np.inf
        cumulative_time += time_delay

        # Check for Cas12a region
        if not found_casregion and DNAstrand[position:position + len(casregion)] == casregion[::-1]:
            found_casregion = True
            print(f"Cas12a region reached at position {len(DNAstrand) - position} after {cumulative_time:.2f} seconds.")
            break

    # Start Michaelis-Menten kinetics
    if found_casregion:
        initial_S_concentration = dNTP_Conc * prob  # Use probability to scale substrate concentration
        model = f"""
        model MichaelisMentenFluorescence
            kcat = 0.55;      // Turnover number (1/s)
            Km = 663e-9;      // Michaelis constant (M)
            E = 25e-9;        // Enzyme concentration (M)
            S = {initial_S_concentration}; // Adjusted substrate concentration
            P = 0.0;          // Initial fluorescence (product concentration)

            v: S -> P; (kcat * E * S) / (Km + S);
            S = {initial_S_concentration};
            P = 0.0;
        end
        """
        rr = te.loadAntimonyModel(model)
        result = rr.simulate(0, 1000, 1000)  # Simulate Michaelis-Menten kinetics
        time = result[:, 0]
        fluorescence = result[:, 2]

        # Plot fluorescence
        plt.plot(time + cumulative_time, fluorescence, label=f"NRTI_Conc = {nrti_conc:.1e} M")

# Finalize plot
plt.title("Michaelis-Menten Kinetics with RT Incorporation")
plt.xlabel("Time (s)")
plt.ylabel("Fluorescence (RFU)")
plt.legend()
plt.grid(True)
plt.show()
