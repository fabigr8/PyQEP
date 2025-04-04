
# Python QEPest: Quantitative Estimation of Pesticide-likeness

This module provides functions to calculate quantitative estimates of pesticide-likeness (QEP), including herbicide-likeness (QEH), insecticide-likeness (QEI), and fungicide-likeness (QEF). The implementation is based on the work described in:

Avram, S., Funar-Timofei, S., Borota, A., Chennamaneni, S.R., Manchala, A.K., & Muresan, S. (2014). Quantitative estimation of pesticide-likeness for agrochemical discovery. Journal of Cheminformatics, 6, 42. https://doi.org/10.1186/s13321-014-0042-6

> [!IMPORTANT]  
> We will hopefully soon integrate this into [RDKIT](https://www.rdkit.org/)
  
## Molecular Descriptors
QEP is a derivative of the QED with adapted descriptors (as parameters of a scoring function) to evaluate pesticide likeliness. It implements three sub-metrics: Quantitative Estimates of Herbicide-likeness (QEH), Insecticide-likeness (QEI), and Fungicide-likeness (QEF).

- At its core, it is a scoring function that uses a set of six molecular descriptors:
  - molecular weight (MW),
  - logP (hydrophobicity),
  - number of hydrogen bond acceptors (HBA),
  - number of hydrogen bond donors (HBD),
  - number of rotatable bonds (RB),
  - and number of aromatic rings (arR).

- The specific weights assigned to each molecular descriptor in the scoring function were determined through curve fitting and an optimization process and are derived from 19,459 pesticides from the GVKBio agrochemical patents collection (AgroSAR).  


## Installation

To use this module, you need to have RDKit installed. You can install RDKit using conda:

```sh
pip install rdkit
```

## Usage

### Importing the Module

```python
from QEP import qeh, qei, qef, qep_max, qep_avg, default
from rdkit import Chem
```

### Calculating QEPest Properties

You can calculate the QEH, QEI, and QEF scores for a molecule using the provided functions. Here is an example:

```python
# Create an RDKit molecule object
smiles = "CCO"
mol = Chem.MolFromSmiles(smiles)

# Calculate QEH, QEI, and QEF scores
qeh_score = qeh(mol)
qei_score = qei(mol)
qef_score = qef(mol)

print(f"QEH: {qeh_score}, QEI: {qei_score}, QEF: {qef_score}")

# Calculate the maximum and average QEPest scores
qep_max_score = qep_max(mol)
qep_avg_score = qep_avg(mol)

print(f"QEP Max: {qep_max_score}, QEP Avg: {qep_avg_score}")

# Calculate default QEPest scores
qeh_default, qei_default, qef_default = default(mol)
print(f"Default QEH: {qeh_default}, QEI: {qei_default}, QEF: {qef_default}")
```

## Functions

- `qeh(mol, w=WEIGHT_DEFAULT, qepest_props=None)`: Calculates Quantitative Estimate of Herbicide-likeness (QEH).
- `qei(mol, w=WEIGHT_DEFAULT, qepest_props=None)`: Calculates Quantitative Estimate of Insecticide-likeness (QEI).
- `qef(mol, w=WEIGHT_DEFAULT, qepest_props=None)`: Calculates Quantitative Estimate of Fungicide-likeness (QEF).
- `qep_max(mol, w=WEIGHT_DEFAULT, qepest_props=None)`: Returns the maximum value among QEH, QEI, and QEF.
- `qep_avg(mol, w=WEIGHT_DEFAULT, qepest_props=None)`: Returns the average value of QEH, QEI, and QEF.
- `default(mol)`: Calculates QEH, QEI, and QEF with default weights and returns them as a tuple.


## Unit test
The implementation provides unittest inlcuding tests comparing outputs to original authors data as well as a created testset with well known herbicides. 

```bash
python -m unittest -v UnitTestQEP.py
```
## Cite

As this repository is not part of the original work, please cite the original work:

```bibtex

@article{avram2014quantitative,
  title={Quantitative estimation of pesticide-likeness for agrochemical discovery},
  author={Avram, Sorin and Funar-Timofei, Simona and Borota, Ana and Chennamaneni, Sridhar Rao and Manchala, Anil Kumar and Muresan, Sorel},
  journal={Journal of Cheminformatics},
  volume={6},
  pages={1--11},
  issn = {1758-2946},
  publisher = {BioMed Central},
  doi = {10.1186/s13321-014-0042-6},
  year={2014},
  publisher={Springer}
}

```

## License

This project is licensed under the MIT License.
