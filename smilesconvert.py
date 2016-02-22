#!/usr/bin/env python
import sys
import pybel

class SmilesConvert:
    '''Convert SMILES string to another molecule file format which openbabel supporting and make a 1500 steps geometry optimisation using MMFF94 forcefield'''
    def __init__(self, smiles, file_format='pdb', filename='main'):
        self.smiles = smiles
        self.format = file_format
        self.filename = filename + '.' + file_format
        mol = pybel.readstring('smi', self.smiles) 
        mol.make3D()
        mol.localopt(forcefield='mmff94', steps=1500)
        output = pybel.Outputfile(self.format, self.filename)
        output.write(mol)
        output.close()
