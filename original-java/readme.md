# Data Creation for Python Unit Test

We used the author's original Java program (.jar) to create test data.  
To do so, we:  
- Used known herbicide structures SMILES strings.  
- Calculated the needed input values for QEP calculation (MW, LogP, HBA, HBD, RB, arR) with RDKit.  
- Calculated the QEP values with the original Java programm.  
- integrated them into the unit test file (line 42 following).

## Create Your Own Test Data

To create your own test data, you:
1. Need an input file such as `input.csv` (with column headings, e.g., SMILES and SMILES Strings in it). 

2. Use the Python file `create_test_data.py` to calculate the QEPest needed inputs in a file, e.g., `original-java/data.csv`, while specifying the column with SMILES Strings. 

~~~bash
python create_test_data.py testdata/input.csv SMILES original-java/data.csv
~~~

3. Run the original Java code to calculate QEP: 
    - First, navigate into the directory where the jar and the data need to be located: `cd original-java`
    - Then start the calculation: `java -jar QEPest.jar`
    - The Java program creates a `data.txt.out` file with the QEP values.

4. Replace the SMILES and the target QEP values used in the `UnitTestQEP.py` (lines 42 following) with your created values.
