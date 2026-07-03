# Peptide Snapshot Extraction and Truncation

## Overview

Molecular dynamics (MD) simulations often produce systems that are too large for direct application in density functional theory (DFT) calculations, such as the aqueous ion-peptide systems in the associated manuscript. To address this, we implement a truncation strategy that:

- Preserves the local chemical environment around the metal ion (in this case, rare-earth elements)
- Reduces system size for computational tractability
- Maintains chemical integrity through C-C bond breakage and hydrogen capping at the broken bonds

This approach generates MD-informed, computationally tractable input structures for DFT calculations using representative snapshots of ion-peptide complexes.

This repository contains Python scripts developed to process GROMACS MD trajectories, extract representative snapshots, and apply a truncation strategy to aqueous ion-peptide systems. The workflow produces chemically consistent truncated models for use as DFT input structures.

## Repository Contents

### `extract-snapshots-from-xtc.py`

Selects specified time points from a GROMACS MD trajectory in XTC format and exports the corresponding atomic coordinates as PDB files.

### `peptide-truncation.py`

Applies a radial-cutoff-based truncation strategy that retains the region of interest around a metal ion. The script removes solvent molecules beyond a specified cutoff distance, truncates the peptide beyond a radial cutoff, and caps broken C-C bonds with hydrogen atoms.

These methods are described in detail in Section S9a of the Supporting Information for the associated manuscript.

## License

This project is licensed under the [MIT License](LICENSE).
