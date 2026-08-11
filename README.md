# EpCAM-Trop2-AlphaPulldown
Scripts and supplementary files used for the in silico analysis of potential interaction partners of EpCAM and Trop2.

# In silico analysis of EpCAM and Trop2 interaction partners

This repository contains scripts and supplementary files used in the diploma thesis:

**In silico prediction of potential interaction partners of EpCAM and Trop2**

University of Ljubljana  
Faculty of Chemistry and Chemical Technology  
Biochemistry

## Analysis

Potential protein-protein interactions involving the extracellular domains of
EpCAM and Trop2 were modeled using AlphaPulldown / AlphaFold-Multimer.

Predicted complexes were evaluated using AFM-LIS, with iLIS used as the main
interaction confidence metric.

Additional custom scripts were used for:
- filtering models according to interaction residues,
- analysis of contact residue frequencies,
- filtering against the dimerization interface,
- comparison with known interactors,
- comparison with proteins identified in Žagar's doctoral dissertation.

## Repository structure

- `scripts/` – Python scripts used for data processing and filtering
- `slurm/` – Slurm job scripts used on the VEGA supercomputer
- `supplementary_data/` – input lists and supplementary analysis files

## Software

- AlphaPulldown
- AlphaFold-Multimer
- AFM-LIS
- LocalColabFold
- UCSF ChimeraX

## Author

Simon Kristl
