"""
Unit tests for QEPest (Quantitative Estimation of Pesticide-likeness) module
"""

import doctest
import os.path
import unittest
from collections import namedtuple
import math

from rdkit import Chem
# Import QEPest module - adjust import path as needed
import QEP as qepest

doLong = False
_TestData = namedtuple('_TestData', 'lineNo,smiles,mol,qeh,qei,qef')

# Test data files - create these files with known pesticides and their expected QEPest values
# The paths should be adjusted to your project structure
TEST_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data')
dataHerbicides = os.path.join(TEST_DIR, 'herbicides_qepest.csv')
dataInsecticides = os.path.join(TEST_DIR, 'insecticides_qepest.csv')
dataFungicides = os.path.join(TEST_DIR, 'fungicides_qepest.csv')

# Test data for sample molecules with expected QEPest values
# Sample data properties and expected values from the authors' supplementary material
SAMPLE_DATA = [
    {"name": "mol1", "MW": 240.2127, "LogP": 3.2392, "HBA": 5, "HBD": 1, "RB": 4, "arR": 1, 
     "QEH": 0.8511, "QEI": 0.5339, "QEF": 0.6224},
    {"name": "mol2", "MW": 249.091, "LogP": 3.0273, "HBA": 3, "HBD": 1, "RB": 5, "arR": 1,
     "QEH": 0.9750, "QEI": 0.6913, "QEF": 0.7310},
    {"name": "mol3", "MW": 308.354, "LogP": 2.1086, "HBA": 1, "HBD": 0, "RB": 7, "arR": 1,
     "QEH": 0.7980, "QEI": 0.9018, "QEF": 0.7320},
    {"name": "mol4", "MW": 360.444, "LogP": 4.0137, "HBA": 3, "HBD": 0, "RB": 8, "arR": 0,
     "QEH": 0.5839, "QEI": 0.8382, "QEF": 0.6594},
    {"name": "mol5", "MW": 295.335, "LogP": 4.9335, "HBA": 2, "HBD": 0, "RB": 1, "arR": 1,
     "QEH": 0.8099, "QEI": 0.8118, "QEF": 0.8742}
]

# Sample data of pesticides  SMILES and expected values calculated via the original QEPest Java implementation
SAMPLE_MOLECULES = [
    # name, smiles, qeh, qei, qef
    ("Glyphosate", "C(C(=O)O)NCP(=O)(O)O", 0.136, 0.2136, 0.2092), 
    ("Glyphosate2", "OCC(O)(P(O)(O)=O)NC(=O)C", 0.1443, 0.174, 0.199), # different representation
    ("2,4-D", "C1=CC(=C(C=C1Cl)Cl)OCC(=O)O", 0.8924, 0.6358, 0.7475),
    ("Imidacloprid", "C1CN(C(=N1)N[N+](=O)[O-])CC2=CN=C(C=C2)Cl", 0.5525, 0.4565, 0.5302),
    ("Chlorpyrifos", "CCOP(=S)(OCC)OC1=NC(=C(C=C1Cl)Cl)Cl", 0.7672, 0.7064, 0.6119),
    ("Azoxystrobin", "CO/C=C(/C(=O)OC)c1ccccc1Oc1cc(Oc2ccccc2C#N)ncn1", 0.4402, 0.3398, 0.2722),
    ("Tebuconazole", "CC(C)(C)C(CCC1=CC=C(C=C1)Cl)(CN2C=NC=N2)O", 0.8732, 0.5928, 0.6851)
]


def load_tests(loader, tests, ignore):
    """ Add the Doctests from the module """
    tests.addTests(doctest.DocTestSuite(qepest, optionflags=doctest.ELLIPSIS))
    return tests


class PropertyMolecule:
    """Helper class to simulate a molecule with predefined property values"""
    def __init__(self, name, mw, logp, hba, hbd, rb, arr):
        self.name = name
        self.mw = mw
        self.logp = logp
        self.hba = hba
        self.hbd = hbd
        self.rb = rb
        self.arr = arr
        self.properties = qepest.QEPestProperties(MW=mw, ALOGP=logp, HBA=hba, HBD=hbd, RB=rb, ARR=arr)

    def __str__(self):
        return self.name


class TestCase(unittest.TestCase):
    
    def testVersion(self):
        """Test that the version is as expected"""
        self.assertEqual(qepest.qeh.version, '1.0.0',
                        msg='QEPest version has changed. Update the regression tests if required.')
    
    def testProperties(self):
        """Test the property calculations on Glyphosate"""
        # Glyphosate
        m = Chem.MolFromSmiles("C(C(=O)O)NCP(=O)(O)O")
        p = qepest.properties(m)
        
        self.assertAlmostEqual(p.MW, 169.07, places=2) # may needs to be reduced to 1 hence floating point error
        self.assertAlmostEqual(p.ALOGP, -1.2 , places=2) # may needs to be reduced to 1 hence floating point error
        self.assertEqual(p.HBA, 3)
        self.assertEqual(p.HBD, 4)
        self.assertEqual(p.RB, 4)
        self.assertEqual(p.ARR, 0)
        
        # Check that adding hydrogens will not change the result
        mol_h = Chem.AddHs(m)
        p_h = qepest.properties(mol_h)
        
        self.assertAlmostEqual(p_h.MW, 169.07, places=2) # may needs to be reduced to 1 hence floating point error
        self.assertAlmostEqual(p_h.ALOGP, -1.2, places=2) # may needs to be reduced to 1 hence floating point error
        self.assertEqual(p_h.HBA, 3)
        self.assertEqual(p_h.HBD, 4)
        self.assertEqual(p_h.RB, 4)
        self.assertEqual(p_h.ARR, 0)
    
    def testSampleMolecules(self):
        """Test QEPest on a set of sample molecules with known values"""
        for name, smiles, expected_qeh, expected_qei, expected_qef in SAMPLE_MOLECULES:
            mol = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(mol, f"Failed to create molecule from SMILES: {smiles}")
            
            # Calculate scores
            qeh_score = qepest.qeh(mol)
            qei_score = qepest.qei(mol)
            qef_score = qepest.qef(mol)
            
            # Test with tolerance
            self.assertAlmostEqual(qeh_score, expected_qeh, places=3, 
                                  msg=f"QEH for {name} ({smiles}) does not match expected value")
            self.assertAlmostEqual(qei_score, expected_qei, places=3,
                                  msg=f"QEI for {name} ({smiles}) does not match expected value") 
            self.assertAlmostEqual(qef_score, expected_qef, places=3,
                                  msg=f"QEF for {name} ({smiles}) does not match expected value")
    
    def testHerbicides(self):
        """Test QEH on known herbicides"""
        if not os.path.exists(dataHerbicides):
            self.skipTest(f"Test data file not found: {dataHerbicides}")
            
        for d in readTestData(dataHerbicides):
            self.assertAlmostEqual(qepest.qeh(d.mol), d.qeh, places=3,
                                  msg=f'QEH not equal to expected in line {d.lineNo}')
            
            # Check that adding hydrogens will not change the result
            mol = Chem.AddHs(d.mol)
            self.assertAlmostEqual(qepest.qeh(mol), d.qeh, places=3,
                                  msg=f'QEH (with Hs) not equal to expected in line {d.lineNo}')
    
    def testInsecticides(self):
        """Test QEI on known insecticides"""
        if not os.path.exists(dataInsecticides):
            self.skipTest(f"Test data file not found: {dataInsecticides}")
            
        for d in readTestData(dataInsecticides):
            self.assertAlmostEqual(qepest.qei(d.mol), d.qei, places=3,
                                  msg=f'QEI not equal to expected in line {d.lineNo}')
            
            # Check that adding hydrogens will not change the result
            mol = Chem.AddHs(d.mol)
            self.assertAlmostEqual(qepest.qei(mol), d.qei, places=3,
                                  msg=f'QEI (with Hs) not equal to expected in line {d.lineNo}')
    
    def testFungicides(self):
        """Test QEF on known fungicides"""
        if not os.path.exists(dataFungicides):
            self.skipTest(f"Test data file not found: {dataFungicides}")
            
        for d in readTestData(dataFungicides):
            self.assertAlmostEqual(qepest.qef(d.mol), d.qef, places=3,
                                  msg=f'QEF not equal to expected in line {d.lineNo}')
            
            # Check that adding hydrogens will not change the result
            mol = Chem.AddHs(d.mol)
            self.assertAlmostEqual(qepest.qef(mol), d.qef, places=3,
                                  msg=f'QEF (with Hs) not equal to expected in line {d.lineNo}')
    
    def testCombinedMetrics(self):
        """Test combined QEP metrics (max and avg)"""
        for name, smiles, expected_qeh, expected_qei, expected_qef in SAMPLE_MOLECULES[:2]:  # Just test first two
            mol = Chem.MolFromSmiles(smiles)
            
            # Calculate max and avg
            qep_max_value = qepest.qep_max(mol)
            qep_avg_value = qepest.qep_avg(mol)
            
            # Check max
            expected_max = max(expected_qeh, expected_qei, expected_qef)
            self.assertAlmostEqual(qep_max_value, expected_max, places=3,
                                 msg=f"QEP_max for {name} does not match expected value")
            
            # Check avg
            expected_avg = (expected_qeh + expected_qei + expected_qef) / 3
            self.assertAlmostEqual(qep_avg_value, expected_avg, places=3,
                                 msg=f"QEP_avg for {name} does not match expected value")
    
    def testCustomWeights(self):
        """Test using custom weights"""
        mol = Chem.MolFromSmiles(SAMPLE_MOLECULES[0][1])  # First molecule
        
        # Define custom weights
        custom_weights = qepest.QEPestProperties(1.5, 1.0, 0.5, 0.5, 1.0, 1.0)
        
        # Calculate with custom weights and default weights
        qeh_custom = qepest.qeh(mol, custom_weights)
        qeh_default = qepest.qeh(mol)
        
        # They should be different
        self.assertNotEqual(qeh_custom, qeh_default,
                           msg="Custom weights didn't change the QEH value")
                           
    def testExactCalculation(self):
        """Test exact calculation steps for a specific property"""
        # Create a PropertyMolecule for the first test sample
        mol = PropertyMolecule(
            SAMPLE_DATA[0]["name"],
            SAMPLE_DATA[0]["MW"],
            SAMPLE_DATA[0]["LogP"],
            SAMPLE_DATA[0]["HBA"],
            SAMPLE_DATA[0]["HBD"],
            SAMPLE_DATA[0]["RB"],
            SAMPLE_DATA[0]["arR"]
        )
        
        # Test detailed calculation for MW of the first molecule
        mw = mol.mw
        df_h = qepest.compute_df(mw, qepest.HERB_PARAMS[0][0], qepest.HERB_PARAMS[0][1], 
                                qepest.HERB_PARAMS[0][2], qepest.HERB_PARAMS[0][3])
        norm_h = qepest.normalize_herb(df_h, 0)
        
        # Verify intermediate calculation results
        self.assertAlmostEqual(df_h / qepest.HERB_MAX[0], norm_h, places=6,
                              msg="Normalization calculation is incorrect")
        
        # Calculate QEH manually to verify each step
        qeh_sum = 0.0
        for i, value in enumerate([mol.mw, mol.logp, mol.hba, mol.hbd, mol.rb, mol.arr]):
            df = qepest.compute_df(value, qepest.HERB_PARAMS[i][0], qepest.HERB_PARAMS[i][1], 
                                 qepest.HERB_PARAMS[i][2], qepest.HERB_PARAMS[i][3])
            normalized_df = qepest.normalize_herb(df, i)
            qeh_sum += qepest.WEIGHT_DEFAULT[i] * math.log(max(normalized_df, 1e-10))
        
        manual_qeh = math.exp(qeh_sum / sum(qepest.WEIGHT_DEFAULT))
        
        # Test with the QEH function
        function_qeh = qepest.qeh(None, qepest.WEIGHT_DEFAULT, mol.properties)
        
        self.assertAlmostEqual(manual_qeh, function_qeh, places=6, 
                              msg="Manual calculation and function calculation differ")
    
    def testSupplementaryData(self):
        """Test against supplementary data from the authors"""
        #print("\nTesting against authors' supplementary data:")
        
        # Test each molecule from SAMPLE_DATA
        for mol_data in SAMPLE_DATA:
            # Create a PropertyMolecule
            mol = PropertyMolecule(
                mol_data["name"],
                mol_data["MW"],
                mol_data["LogP"],
                mol_data["HBA"],
                mol_data["HBD"],
                mol_data["RB"],
                mol_data["arR"]
            )
            
            # Calculate scores
            qeh_score = qepest.qeh(None, qepest.WEIGHT_DEFAULT, mol.properties)
            qei_score = qepest.qei(None, qepest.WEIGHT_DEFAULT, mol.properties)
            qef_score = qepest.qef(None, qepest.WEIGHT_DEFAULT, mol.properties)
            
            # Compare with expected values 
            # # allowing small tolerance due difference in FP calcs between implementations
            self.assertAlmostEqual(qeh_score, mol_data["QEH"], places=2,
                                  msg=f"QEH for {mol.name} does not match expected value")
            self.assertAlmostEqual(qei_score, mol_data["QEI"], places=2,
                                  msg=f"QEI for {mol.name} does not match expected value")
            self.assertAlmostEqual(qef_score, mol_data["QEF"], places=2,
                                  msg=f"QEF for {mol.name} does not match expected value")


def readTestData(filename):
    """ Read test data from file """
    with open(filename, 'r') as f:
        for lineNo, line in enumerate(f, 1):
            if line[0] == '#':
                continue
            parts = line.strip().split(',')
            if len(parts) == 4:  # smiles,qeh,qei,qef
                smiles, qeh, qei, qef = parts
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    raise AssertionError(f'Molecule construction failed on line {lineNo}')
                yield _TestData(lineNo, smiles, mol, float(qeh), float(qei), float(qef))


def createTestData():
    """Create test data files with known pesticides"""
    os.makedirs(TEST_DIR, exist_ok=True)
    
    # Known herbicides
    herbicides = [
        r"CCNC1=NC(=NC(=N1)Cl)NC(C)C ",  # Atrazine PubChem
        r"CCNc1nc(Cl)nc(NC(C)C)n1"       # Astracine CHEMBL
        r"C1=CC(=C(C=C1Cl)Cl)OCC(=O)O",  # 2,4-D PubChem
        r"C(C(=O)O)NCP(=O)(O)O",  # Glyphosate  PubChem
        r"C1=CC(=CC=C1[N+](=O)[O-])OC2=C(C=C(C=C2)Cl)Cl",  # Nitrofen PubChem
    ]
    
    # Known insecticides
    insecticides = [
        r"C1CN(C(=N1)N[N+](=O)[O-])CC2=CN=C(C=C2)Cl",  # Imidacloprid PubChem
        r"O=[N+]([O-])/N=C1\NCCN1Cc1ccc(Cl)nc1"        #Imidacloprid CHEMBL
        r"CCOP(=S)(OCC)OC1=NC(=C(C=C1Cl)Cl)Cl",  # Chlorpyrifos PubChem
        r"C1=C(C=C(C(=C1Cl)N2C(=C(C(=N2)C#N)S(=O)C(F)(F)F)N)Cl)C(F)(F)F",  # Fipronil PubChem
        r"COC1=NN(C(=O)S1)CSP(=S)(OC)OC",  # Methidation PubChem
    ]
    
    # Known fungicides
    fungicides = [
        r"CO/C=C(\C1=CC=CC=C1OC2=NC=NC(=C2)OC3=CC=CC=C3C#N)/C(=O)OC ",  # Azoxystrobin PubChem: 3034285
        r"CO/C=C(/C(=O)OC)c1ccccc1Oc1cc(Oc2ccccc2C#N)ncn1"   #Azoxystrobin CHEMBL230001 values should be same as above
        r"CC(C)(C)C(CCC1=CC=C(C=C1)Cl)(CN2C=NC=N2)O",  # Tebuconazole PubChem
        r"C(C(Cl)(Cl)Cl)(NC=O)N1CCN(C(C(Cl)(Cl)Cl)NC=O)CC1",  # Triforine CAS SMILES
        r"CCCC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)Cl)Cl",  # Propiconazole PubChem
    ]
    
    # Create test data files
    with open(dataHerbicides, 'w') as f:
        f.write("# Test data for QEH descriptor\n")
        for smiles in herbicides:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                qeh_val = qepest.qeh(mol)
                qei_val = qepest.qei(mol)
                qef_val = qepest.qef(mol)
                f.write(f"{smiles},{qeh_val:.6f},{qei_val:.6f},{qef_val:.6f}\n")
    
    with open(dataInsecticides, 'w') as f:
        f.write("# Test data for QEI descriptor\n")
        for smiles in insecticides:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                qeh_val = qepest.qeh(mol)
                qei_val = qepest.qei(mol)
                qef_val = qepest.qef(mol)
                f.write(f"{smiles},{qeh_val:.6f},{qei_val:.6f},{qef_val:.6f}\n")
    
    with open(dataFungicides, 'w') as f:
        f.write("# Test data for QEF descriptor\n")
        for smiles in fungicides:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                qeh_val = qepest.qeh(mol)
                qei_val = qepest.qei(mol)
                qef_val = qepest.qef(mol)
                f.write(f"{smiles},{qeh_val:.6f},{qei_val:.6f},{qef_val:.6f}\n")


if __name__ == '__main__':  # pragma: nocover
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', default=False, action='store_true', dest='doLong',
                        help='Run long tests')
    parser.add_argument('-c', default=False, action='store_true', dest='createTestData',
                        help='Create test data files')
    args = parser.parse_args()

    # Handle possible arguments
    doLong = args.doLong
    if args.doLong:
        sys.argv.remove('-l')

    if args.createTestData:
        createTestData()
        sys.argv.remove('-c')

    unittest.main()