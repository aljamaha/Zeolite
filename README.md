# Zeolite

In Progress ..

## Requirements:
  - molmod
  - ase 
  - Q-Chem

## generate_structures: 
### objective: 
Generate unique zeolite structures (1 or 2 Al atoms substituting Si, and adsorption sites for H/metal)
### syntax:
`python main.py`
### inputs:
- xyz file of original zeolite structure (placed under original_structures folders)
- index of T-atom (Si atoms) to be replaced with the first Al atom
- number of terminal H atoms in the original structure (needed to eliminate dangling bonds from terminal Si in the cluster)
- adsorbed metals oxidation states and compositions

## Missing:
currently only adds adsorbates with oxidation state of +2

## outputs:
- generate structures of the zeolites (under structures_saved)
- data file (json) containing information on each generated structure (under data directory)


## Create input
    
    
    


