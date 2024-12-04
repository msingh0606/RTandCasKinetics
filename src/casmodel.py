import numpy as np
import tellurium as te
import matplotlib.pyplot as plt
from probability_func import calculate_dNTP_probability

# INPUT PARAMETERS
forwardDNA = 'TTTTTTTTTTTTTGATGATGTGAAGGTGTTGTCGTTTATTTATTTATTTATTTATTTCTATCTTTCCTCTTAATTCGACG'
DNAstrand = forwardDNA[::-1]
casregionforward = 'ATGATGTGAAGGTGTTGTCG'
casregion = casregionforward[::-1]
PrimerConc = 50e-9  # Primer concentration (M)
dNTPConc_nM = 100
dNTP_Conc = dNTPConc_nM * 1e-9  # dNTP concentration (M)
NRTI_Conc = np.logspace(-10, -5, 6)  # varying NRTI concentrations
TemplateConc = 5e-9  # Template Concentration (M)
TemplateCopies = TemplateConc * 6.022e23 * 1000  # copies of template in 1uL reaction
Kaff = 0.3  # K_aff of TFV-DP, the drug we use

# Reaction constants
nucleotide_addition_rate = 5  # Basepairs per second RT incorporates
kcat = 0.55  # Michaelis-Menten turnover number (1/s)
Km = 663e-9  # Michaelis constant (M)
E = 5e-9  # Enzyme concentration (M)

# Simulation parameters
time_total = 7200  # Total simulation time (seconds)

# Estimate primer binding time (diffusion-limited)
k_diff = 1e5  # Diffusion-limited rate constant (M^-1 s^-1)
t_bind = 1 / (k_diff * PrimerConc)  # time of binding
print(f"Estimated primer binding time: {t_bind:.2e} seconds")

# Initialize plot
plt.figure(figsize=(10, 6))

for nrti_conc in NRTI_Conc:
    print(f"Simulating for NRTI_Conc = {nrti_conc:.1e} M")

    cumulative_time = t_bind  # Start with primer binding time
    cumulative_probability = 1.0  # Start with full probability
    reached_casregion = False
    nucleotides_added = 23  # Start synthesis at the 24th nucleotide

    while cumulative_time <= time_total and not reached_casregion:
        if nucleotides_added < len(DNAstrand):
            # Calculate probability for nucleotide addition
            prob = calculate_dNTP_probability(
                base=DNAstrand[nucleotides_added],
                dNTP_Conc=dNTP_Conc,
                NRTI_Conc=nrti_conc,
                Kaff=Kaff,
            )
            prob = max(prob, 1e-6)  # Prevent division by zero or unrealistic values
            cumulative_probability *= prob

            # Calculate nucleotide addition time delay
            nucleotide_delay = (1 / nucleotide_addition_rate) / prob
            cumulative_time += nucleotide_delay

            # Check if Cas region is reached
            if DNAstrand[nucleotides_added:nucleotides_added + len(casregion)] == casregion:
                reached_casregion = True
                print(f"Cas region reached after {cumulative_time:.2f} seconds.")

            nucleotides_added += 1
        else:
            break  # No more nucleotides to process

    # Simulate Michaelis-Menten fluorescence with Tellurium
    if reached_casregion:
        scaled_kcat = kcat * cumulative_probability  # Scale by probability
        model = f"""
        model MichaelisMentenFluorescence
            kcat = {scaled_kcat};      // Turnover number (1/s)
            Km = {Km};         // Michaelis constant (M)
            E = {E};           // Enzyme concentration (M)
            S = 5e-9;  // Initial substrate concentration (M)
            P = 0.0;           // Initial fluorescence (product concentration)

            v: S -> P; (kcat * E * S) / (Km + S);

            S = 5e-9;
            P = 0.0;
        end
        """
        rr = te.loadAntimonyModel(model)
        result = rr.simulate(0, time_total - cumulative_time, 1000)

        # Extract time and fluorescence from Tellurium results
        time = result[:, 0] + cumulative_time  # Adjust time offset by cumulative time
        fluorescence = result[:, 2]  # Product concentration (fluorescence)

        # Scale fluorescence for the total number of templates
        fluorescence *= TemplateCopies

    else:
        # If Cas region is not reached, generate a flat zero fluorescence curve
        time = np.linspace(0, time_total, 1000)
        fluorescence = np.zeros_like(time)
    time_mins = time/60 # convert to minutes
    # Plot fluorescence curve
    plt.plot(time_mins, fluorescence, label=f"NRTI_Conc = {nrti_conc:.1e}")

# Finalize plot
plt.title("Nucleotide Addition and Michaelis-Menten Fluorescence")
plt.xlabel("Time (minutes)")
plt.ylabel("Fluorescence (RFU)")
plt.xlim(0, 120)
plt.legend()
plt.grid(True)
plt.show()
