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
TemplateConc = 5e-9  # Template concentration in M
TemplateSample = int(TemplateConc * 6.022e23 / 1000)  # Template copies in 1 µL sample

# Define nucleotide addition rate
nucleotide_addition_rate = 0.01  # Basepairs per second

# Initialize plot
plt.figure(figsize=(10, 6))

# Simulate for each NRTI concentration
for nrti_conc in NRTI_Conc:
    print(f"Simulating for NRTI_Conc = {nrti_conc:.1e} M")
    cumulative_time = 0  # Reset cumulative time for each concentration
    found_casregion = False

    # RT synthesis simulation with a time limit
    for position, base in enumerate(DNAstrand):
        if position < 23:
            # Skip the first 23 nucleotides
            cumulative_time += 1 / nucleotide_addition_rate
            continue

        # Calculate probability with amplified impact
        prob = calculate_dNTP_probability(base, dNTP_Conc, nrti_conc, Kaff, n=position - 23 + 1)
        prob = max(prob**2, 1e-6)  # Amplify impact of NRTI_Conc on probability

        # Adjust time delay with amplified probability
        time_delay = (1 / nucleotide_addition_rate) / prob
        cumulative_time += time_delay

        # Debugging: Print key values
        print(f"NRTI_Conc: {nrti_conc:.1e}, Position: {position}, Probability: {prob:.6f}, Time Delay: {time_delay:.6f}, Cumulative Time: {cumulative_time:.2f}")

        # Check for Cas12a region
        if not found_casregion and DNAstrand[position:position + len(casregion)] == casregion[::-1]:
            found_casregion = True
            print(f"Cas12a region reached at position {len(DNAstrand) - position} after {cumulative_time:.2f} seconds.")
            break

        # Stop if cumulative time exceeds 120 minutes
        if cumulative_time > 7200:
            print(f"Cas region not reached within 120 minutes for NRTI_Conc = {nrti_conc:.1e}.")
            break

    if found_casregion:
        # Start Michaelis-Menten kinetics
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
        result = rr.simulate(0, 7200, 1000)  # Simulate for 120 minutes (7200 seconds)
        time = result[:, 0]
        fluorescence = result[:, 2]

        # Scale fluorescence by the total number of template copies
        total_fluorescence = fluorescence * TemplateSample
    else:
        # If Cas region is not reached, set total fluorescence to 0
        time = np.linspace(0, 7200, 1000)  # Generate a flat time array
        total_fluorescence = np.zeros_like(time)

    # Plot scaled or flat fluorescence
    plt.plot(time, total_fluorescence, label=f"NRTI_Conc = {nrti_conc:.1e} M")

# Finalize plot with 120-minute bounds
plt.title("Michaelis-Menten Kinetics with RT Incorporation for Scaled Template Copies")
plt.xlabel("Time (s)")
plt.ylabel("Fluorescence (RFU)")
plt.xlim(0, 7200)  # Restrict x-axis to 120 minutes
plt.legend()
plt.grid(True)
plt.show()
