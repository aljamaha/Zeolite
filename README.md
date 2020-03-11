# Zeolite

In Progress ..

## Requirements:
  - molmod
  - ase 
  - Q-Chem

## Step 1: generate_structures: 

### objective: 
Generate unique zeolite structures (1 or 2 Al atoms substituting Si, and compensating cations/protons)
### syntax:
`python main.py`
- functions.py [common useful functions]
- adsorbate.py [scripts specific to adding adsorbates]
### options:
- choice of adsorbate (H, NH3, Pd)
### inputs:
- xyz file of original zeolite structure (placed under original_structures folders)
- index of T-atom (Si atoms) to be replaced with the first Al atom
- number of terminal H atoms in the original structure (needed to eliminate dangling bonds from terminal Si in the cluster)
- adsorbed metals oxidation states and compositions

## outputs:
- generate structures of the zeolites (under structures_saved)
- data file (json) containing information on each generated structure (under data directory)

## Step 2: Create input
### objective: 
Create Q-Chem input files for desired structures
### syntax:
`python create_input.py`
### inputs:
- index of the first structure
- index of the last structure
- calculation type (optimization or single-point)
## outputs:
generate input file for each calculation in a new folder (under calculations directory)

## Step 4: Calculations
Run desired calculations (under calculations directory)

## Step 5: 

