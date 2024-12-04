def calculate_dNTP_probability(base, dNTP_Conc, NRTI_Conc, Kaff):
    """
    Calculate the probability P_dNTP of incorporating a single dNTP.

    Parameters:
    - base: The nucleotide base ('A', 'T', 'C', 'G')
    - dNTP_Conc: Uniform concentration of dNTP
    - NRTI_Conc: Concentration of NRTI (only affects T)
    - Kaff: Affinity factor of NRTI

    Returns:
    - Probability P_dNTP
    """
    # Adjust fraction if the base is T
    if base == "T":
        fraction = dNTP_Conc / (dNTP_Conc + Kaff * NRTI_Conc)
    else:
        fraction = 1.0  # No competition for non-T bases
    return fraction
