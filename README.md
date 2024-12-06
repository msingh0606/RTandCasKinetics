# RTandCasKinetics
This is a tool that is used for creating a plot that investigates the combined kinetics of HIV reverse transcriptase

## Users
Users do not need to be familiar with python to use this program. The users will use the Graphical User Interface (GUI) to input variables reducing the need to be familiar with python

## Installation
* To install, first install `pdm`, and make sure it's accessible in your $PATH.
* You can install it by going to 'https://pdm-project.org/en/latest/' and following their instructions to download PDM.
Prepare the local python environment by running `pdm install` and follow the instructions provided

## Contributing
To add new packages (e.g. like pip) to your local python environment (like an anaconda environment) run `pdm add matplotlib` for example if you wanted to add matplotlib. 

## Documents
The documents contain three main folders: docs, src, and tests.
* Docs contain the functional and component specifications and a background presentation
* Src contains the functions (casmodel_func.py and probability_func.py) and the GUI interface (casmodel.py)
* Tests contain the two test modules for the functions: casmodel_func.py and probability_func.py.

## Getting Started
* To access the package, go to your terminal and type 
```
pip install RTandCasKinetics
```
* If you cannot find the root directory, enter the following code: `find ~/ -name "RTandCasKinetics"`. It will give many things, but the last will show the root.
* To go to the root directory type in the results from the search and type `cd enter_directory_here`
* To run this, go to the root directory `RTandCasKinetics/` and run the command 
```
pdm run src/rtandcaskinetics/casmodel.py
```

## User Guide:
If you are using a jupyter notebook, use the following code which does not have an interactive interface:
```
import matplotlib.pyplot as plt
from rtandcaskinetics.casmodel_func import compute_fluorescence

# Define the function to simulate fluorescence and plot
def run_simulation(forwardDNA, TemplateConc_nM, PrimerConc_nM, dNTPConc_nM):
    if not forwardDNA or TemplateConc_nM <= 0 or PrimerConc_nM <= 0 or dNTPConc_nM <= 0:
        raise ValueError("All inputs must be valid non-zero values.")

    if not all(base in "ATCG" for base in forwardDNA.upper()):
        raise ValueError("Forward DNA sequence contains invalid characters.")

    # Simulate fluorescence
    results = compute_fluorescence(forwardDNA, TemplateConc_nM, PrimerConc_nM, dNTPConc_nM)

    # NRTI concentrations
    nrti_concs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]

    # Plot the results
    plt.figure(figsize=(8, 6))
    for (time_mins, fluorescence), nrti_conc in zip(results, nrti_concs):
        plt.plot(time_mins, fluorescence, label=f"NRTI_Conc = {nrti_conc:.1e}")
    plt.title("Combined Kinetics of RT Incorporation and Cas12a Cleavage")
    plt.xlabel("Time (minutes)")
    plt.ylabel("Fluorescence (RFU)")
    plt.xlim(0, 120)
    plt.legend()
    plt.grid(True)
    plt.show()

# Example inputs (adapt these to your needs)
forwardDNA = "TTTTTTTTTTTTTGATGATGTGAAGGTGTTGTCGTTTATTTATTTATTTATTTATTTCTATCTTTCCTCTTAATTCGACG"
TemplateConc_nM = 5.0
PrimerConc_nM = 50.0
dNTPConc_nM = 100.0

# Run the simulation
run_simulation(forwardDNA, TemplateConc_nM, PrimerConc_nM, dNTPConc_nM)
```

* If using terminal, follow the steps above