"""
QEPest: Quantitative Estimation of Pesticide-likeness

This module provides functions to calculate quantitative estimates and Pesticides (QEP); 
more specific quantitative estimates herbicide-likeness (QEH), insecticide-likeness (QEI), and fungicide-likeness (QEF)
as described in:

Avram, S., Funar-Timofei, S., Borota, A., Chennamaneni, S.R., Manchala, A.K., & Muresan, S. (2014). 
Quantitative estimation of pesticide-likeness for agrochemical discovery. 
Journal of Cheminformatics, 6, 42. https://doi.org/10.1186/s13321-014-0042-6

The QEPest metrics use six molecular descriptors:
- Molecular weight (MW)
- LogP (hydrophobicity)
- Number of hydrogen bond acceptors (HBA)
- Number of hydrogen bond donors (HBD)
- Number of rotatable bonds (RB)
- Number of aromatic rings (ARR)

This implementation is structured to be compatible with RDKit's QED implementation style.
"""

import math
from collections import namedtuple

from rdkit import Chem
from rdkit.Chem import Crippen, Lipinski
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit.Chem.ChemUtils.DescriptorUtilities import setDescriptorVersion

# Define the properties needed for QEPest
QEPestProperties = namedtuple('QEPestProperties', 'MW,ALOGP,HBA,HBD,RB,ARR')

# Define weighting schemes
WEIGHT_DEFAULT = QEPestProperties(1.00, 1.00, 1.00, 1.00, 1.00, 1.00)  # Equal weights
WEIGHT_CUSTOM = QEPestProperties(0.66, 0.75, 0.55, 0.40, 0.60, 0.70)   # Custom weights

# Parameters for herbicide-likeness (QEH)
# Format: [a, b, c, o] for each descriptor
HERB_PARAMS = [
    [70.77, 283.0, 84.97, -1.185],    # MW
    [93.81, 3.077, 1.434, 0.6164],    # ALOGP
    [117.6, 2.409, 1.567, 7.155],     # HBA
    [233.4, 0.4535, -1.48, 4.47],     # HBD
    [84.7, 4.758, -2.423, 5.437],     # RB
    [301.8, 1.101, 0.8869, -22.81]    # ARR
]

# Parameters for insecticide-likeness (QEI)
# Format: [a, b, c, o] for each descriptor
INSECT_PARAMS = [
    [76.38, 298.3, 83.64, 1.912],      # MW
    [74.27, 4.555, -2.193, -2.987],    # ALOGP
    [139.4, 1.363, 1.283, 0.5341],     # HBA
    [670.6, -1.163, 0.7856, 0.7951],   # HBD
    [65.49, 6.219, -2.448, 5.318],     # RB
    [287.5, 0.305, 1.554, -88.64]      # ARR
]

# Parameters for fungicide-likeness (QEF)
# Format: [a, b, c, o] for each descriptor
FUNG_PARAMS = [
    [51.03, 314.2, -56.31, 2.342],     # MW
    [50.73, 3.674, -1.238, 2.067],     # ALOGP
    [73.79, 1.841, 1.326, 0.5158],     # HBA
    [164.7, -0.9762, -2.027, 1.384],   # HBD
    [40.91, 1.822, 2.582, 0.6235],     # RB
    [134.4, 0.8383, 1.347, -31.17]     # ARR
]

# Normalization factors for each class and descriptor
HERB_MAX = [69.5849922, 94.4228257, 120.4572352, 228.1589796, 89.7012502, 276.9634213]
INSECT_MAX = [78.2919965, 71.2829691, 133.9224801, 331.170104, 70.5540709, 193.0023343]
FUNG_MAX = [53.3719946, 52.773116, 73.7976536, 144.9887053, 41.4385926, 102.3024319]


def compute_df(x, a, b, c, o):
    """
    Calculates the desirability function value based on the formula:
    f = a * exp(-exp(-(x-b)/c) - (x-b)/c + 1) + o
    
    Args:
        x: Property value
        a, b, c, o: Parameters for the desirability function
        
    Returns:
        The desirability function value
    """
    return a * math.exp(-1.0 * math.exp(-1.0 * ((x - b) / c)) - (x - b) / c + 1.0) + o


def normalize_herb(d, descriptor_idx):
    """Normalize herbicide desirability value"""
    return d / HERB_MAX[descriptor_idx]


def normalize_insect(d, descriptor_idx):
    """Normalize insecticide desirability value"""
    return d / INSECT_MAX[descriptor_idx]


def normalize_fung(d, descriptor_idx):
    """Normalize fungicide desirability value"""
    return d / FUNG_MAX[descriptor_idx]


def properties(mol):
    """
    Calculates the properties required for QEPest.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        QEPestProperties named tuple with calculated properties
    """
    if mol is None:
        raise ValueError('You need to provide a mol argument.')
    
    mol = Chem.RemoveHs(mol)
    
    qepest_props = QEPestProperties(
        MW=rdmd._CalcMolWt(mol),
        ALOGP=Crippen.MolLogP(mol),
        HBA=rdmd.CalcNumHBA(mol),
        HBD=rdmd.CalcNumHBD(mol),
        RB=rdmd.CalcNumRotatableBonds(mol),
        ARR=Lipinski.NumAromaticRings(mol)
    )
    
    return qepest_props


@setDescriptorVersion(version='1.0.0')
def qeh(mol, w=WEIGHT_DEFAULT, qepest_props=None):
    """
    Calculates Quantitative Estimate of Herbicide-likeness (QEH).
    
    Args:
        mol: RDKit molecule object
        w: Weights for each property
        qepest_props: Pre-calculated properties (optional)
        
    Returns:
        QEH score (0-1)
    """
    if qepest_props is None:
        qepest_props = properties(mol)
    
    prop_values = list(qepest_props)
    qeh_sum = 0.0
    
    for i, value in enumerate(prop_values):
        df = compute_df(value, HERB_PARAMS[i][0], HERB_PARAMS[i][1], 
                        HERB_PARAMS[i][2], HERB_PARAMS[i][3])
        normalized_df = normalize_herb(df, i)
        qeh_sum += w[i] * math.log(max(normalized_df, 1e-10))
    
    return math.exp(qeh_sum / sum(w))


@setDescriptorVersion(version='1.0.0')
def qei(mol, w=WEIGHT_DEFAULT, qepest_props=None):
    """
    Calculates Quantitative Estimate of Insecticide-likeness (QEI).
    
    Args:
        mol: RDKit molecule object
        w: Weights for each property
        qepest_props: Pre-calculated properties (optional)
        
    Returns:
        QEI score (0-1)
    """
    if qepest_props is None:
        qepest_props = properties(mol)
    
    prop_values = list(qepest_props)
    qei_sum = 0.0
    
    for i, value in enumerate(prop_values):
        df = compute_df(value, INSECT_PARAMS[i][0], INSECT_PARAMS[i][1], 
                       INSECT_PARAMS[i][2], INSECT_PARAMS[i][3])
        normalized_df = normalize_insect(df, i)
        qei_sum += w[i] * math.log(max(normalized_df, 1e-10))
    
    return math.exp(qei_sum / sum(w))


@setDescriptorVersion(version='1.0.0')
def qef(mol, w=WEIGHT_DEFAULT, qepest_props=None):
    """
    Calculates Quantitative Estimate of Fungicide-likeness (QEF).
    
    Args:
        mol: RDKit molecule object
        w: Weights for each property
        qepest_props: Pre-calculated properties (optional)
        
    Returns:
        QEF score (0-1)
    """
    if qepest_props is None:
        qepest_props = properties(mol)
    
    prop_values = list(qepest_props)
    qef_sum = 0.0
    
    for i, value in enumerate(prop_values):
        df = compute_df(value, FUNG_PARAMS[i][0], FUNG_PARAMS[i][1], 
                       FUNG_PARAMS[i][2], FUNG_PARAMS[i][3])
        normalized_df = normalize_fung(df, i)
        qef_sum += w[i] * math.log(max(normalized_df, 1e-10))
    
    return math.exp(qef_sum / sum(w))


def qep_max(mol, w=WEIGHT_DEFAULT, qepest_props=None):
    """
    Returns the maximum value among QEH, QEI, and QEF.
    
    Args:
        mol: RDKit molecule object
        w: Weights for each property
        qepest_props: Pre-calculated properties (optional)
        
    Returns:
        Maximum of QEH, QEI, and QEF scores
    """
    if qepest_props is None:
        qepest_props = properties(mol)
    
    qeh_value = qeh(mol, w, qepest_props)
    qei_value = qei(mol, w, qepest_props)
    qef_value = qef(mol, w, qepest_props)
    
    return max(qeh_value, qei_value, qef_value)


def qep_avg(mol, w=WEIGHT_DEFAULT, qepest_props=None):
    """
    Returns the average value of QEH, QEI, and QEF.
    
    Args:
        mol: RDKit molecule object
        w: Weights for each property
        qepest_props: Pre-calculated properties (optional)
        
    Returns:
        Average of QEH, QEI, and QEF scores
    """
    if qepest_props is None:
        qepest_props = properties(mol)
    
    qeh_value = qeh(mol, w, qepest_props)
    qei_value = qei(mol, w, qepest_props)
    qef_value = qef(mol, w, qepest_props)
    
    return (qeh_value + qei_value + qef_value) / 3


def default(mol):
    """
    Calculates QEH, QEI, and QEF with default weights and returns them as a tuple.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Tuple of (QEH, QEI, QEF) scores
    """
    props = properties(mol)
    return (qeh(mol, WEIGHT_DEFAULT, props),
            qei(mol, WEIGHT_DEFAULT, props),
            qef(mol, WEIGHT_DEFAULT, props))
