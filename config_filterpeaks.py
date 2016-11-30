# -*- coding: utf-8 -*-

"""
FilterPeaks Python configuration file.

version  : 1.0.0
author   : Marc van Dijk
status   : release beta1
date     : 1 January 2017
copyright: Copyright (c) 2017 Marc van Dijk

This Python style configuration file re-defines the default
FilterPeaks program configuration and defines the residue 
and atom translation scheme required to convert the often
ambiguous (Sparky) residue and atom naming to the standard
IUPAC PDB naming convention. This enables the matching 
between peak assignments and atoms in the structure.

This is the default configuration file used by filterPeaks.py
as defined in the CONFIG_FILE variable. By default is should
be available in the same directory as the main program unless
specified otherwise in the CONFIG_FILE variable.
You may edit this file or reserve a copy for your custom
configuration and load it using the -i/--config command line
option.

NOTE: 
the Python style configuration file uses a Python dictionary
to store all settings. A dictionary name other than 'settings' 
or a syntax error will not load the file.

Residue and atom name translation:
The FilterPeaks program uses a two step translation proces;
first the residue name is translated followed by the atom name
in the scope of the translated residue name. If the residue 
translation was unsuccessful, the atom translation will fallback
to general, residue non-specific, atom name translations if 
defined in the 'atom_translation' dictionary.
If the translation was unsuccessful the non-translated name 
will be used.
The translation process has support for regular expression
based matching, activated using the -r/--use_regex_matching
command line option.

The translation process supports ambiguity by definition of
a list of translated residue or atom names. The residue name
'DX' for example represents all four standard DNA residue names
defined as ['THY', 'ADE', 'CYT', 'GUA'].
The program duplicates the assignment for 'DX' for every atom
belonging to each of the non-ambiguous residue names.

The residue and atom translation rules are based on atomnames.py
file of the Sparky program.
"""

settings = {
  
  # Chemical shift and atom-atom distance tolerance settings
  'chem_shift_ctol'           : 0.5,       # Carbon chemical shift tolerance in PPM
  'chem_shift_ptol'           : 0.05,      # Proton chemical shift tolerance in PPM
  'struc_dist_tol'            : 6.0,       # Distance cutoff in A for atom-atom distance caluclation
  
  # Output settings
  'loglevel'                  : 'warning', # Logging verbosity level 'critical','error','warning','info','debug'
   
  # Residue and atom translation
  'use_atom_translation'      : False,
  'use_residue_translation'   : False,
  'use_regex_matching'        : False,     # Use regular expression matching in residue and atom translation.
  
  # Residue name translation dictionary
  'residue_translation': {
    'A':   ['ALA'], 'R':   ['ARG'], 'N':   ['ASN'], 'D':   ['ASP'], 'B':   ['ASX'],
    'C':   ['CYS'], 'E':   ['GLU'], 'Q':   ['GLN'], 'Z':   ['GLX'], 'G':   ['GLY'],
    'H':   ['HIS'], 'I':   ['ILE'], 'L':   ['LEU'], 'K':   ['LYS'], 'M':   ['MET'],
    'F':   ['PHE'], 'P':   ['PRO'], 'S':   ['SER'], 'T':   ['THR'], 'W':   ['TRP'],
    'Y':   ['TYR'], 'V':   ['VAL'],
    'ALA': ['ALA'], 'ARG': ['ARG'], 'ASN': ['ASN'], 'ASP': ['ASP'], 'ASX': ['ASX'],
    'CYS': ['CYS'], 'GLU': ['GLU'], 'GLN': ['GLN'], 'GLX': ['GLX'], 'GLY': ['GLY'],
    'HIS': ['HIS'], 'ILE': ['ILE'], 'LEU': ['LEU'], 'LYS': ['LYS'], 'MET': ['MET'],
    'PHE': ['PHE'], 'PRO': ['PRO'], 'SER': ['SER'], 'THR': ['THR'], 'TRP': ['TRP'],
    'TYR': ['TYR'], 'VAL': ['VAL'],
    'DT':  ['THY'], 'DA':  ['ADE'], 'DC':  ['CYT'], 'DG':  ['GUA'],
    'THY': ['THY'], 'ADE': ['ADE'], 'CYT': ['CYT'], 'GUA': ['GUA'],
  },
  
  # Protein residues
  'protein_residues': ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS',
                       'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
  
  # Atom name translation dictionary
  'atom_translation': {
    'ADE': {"C1'": ["C1'"], 'C2': ['C2'], "C2'": ["C2'"], "C3'": ["C3'"], 'C4': ['C4'], "C4'": ["C4'"], 'C5': ['C5'], 
            "C5'": ["C5'"], 'C6': ['C6'], 'C8': ['C8'], "H1'": ["H1'"], 'H2': ['H2'], "H2'2": ["H2'2"], "H2'1": ["H2'1"], 
            "H2'": ["H2'1"], "H3'": ["H3'"], "H4'": ["H4'"], "H5'2": ["H5'2"], "H5'1": ["H5'1"], 'H61': ['H61'],
            "H5'": ["H5'1"], 'H5"': ["H5'2"], "H5''": ["H5'2"], "H2''": ["H2'2"], 'H2"': ["H2'2"],
            'H62': ['H62'], 'H8': ['H8'], 'N1': ['N1'], 'N3': ['N3'], 'N6': ['N6'], 'N7': ['N7'], 'N9': ['N9'], 
            'O1P': ['O1P'], "O2'": ["O2'"], 'O2P': ['O2P'], "O3'": ["O3'"], "O4'": ["O4'"], "O5'": ["O5'"], 'P': ['P']},
    'ALA': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'H': ['H'], 'HA': ['HA'], 'HB1': ['HB1'], 'HB2': ['HB2'], 'HB3': ['HB3'], 
            'N': ['N'], 'MB':['HB1', 'HB2', 'HB3'], 'H1': ['H']},
    'ARG': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD': ['CD'], 'CG': ['CG'], 'CZ': ['CZ'], 'H': ['H'], 'HA': ['HA'], 
            'HB2': ['HB2'], 'HB3': ['HB3'], 'HD2': ['HD2'], 'HD3': ['HD3'], 'HE': ['HE'], 'HG2': ['HG2'], 'HG3': ['HG3'], 
            'HH11': ['HH11'], 'HH12': ['HH12'], 'HH21': ['HH21'], 'HH22': ['HH22'], 'N': ['N'], 'NE': ['NE'], 'NH1': ['NH1'],
            'NH2': ['NH2'], 'NQH': ['NH1', 'NH2'], 'QB': ['HB2', 'HB3'], 'QD': ['HD1', 'HD2'], 'QG': ['HG2', 'HG3'], 
            'QH1': ['HH11', 'HH12'], 'QH2': ['HH21', 'HH22'], 'QQH': ['HH11', 'HH12', 'HH21', 'HH22']},
    'ASN': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CG': ['CG'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'], 
            'HD21': ['HD21'], 'HD22': ['HD22'], 'N': ['N'], 'ND2': ['ND2'], 'QB': ['HB2', 'HB3'], 'QD2': ['HD21', 'HD22'], 
            'H1': ['H']},
    'ASP': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CG': ['CG'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'],
            'HD2': ['HD2'], 'N': ['N'], 'QB': ['HB2', 'HB3'], 'H1': ['H']},
    'CYS': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'], 'HG': ['HG'], 
            'N': ['N'], 'QB':['HB2', 'HB3'], 'H1': ['H']},
    'CYT': {"C1'": ["C1'"], 'C2': ['C2'], "C2'": ["C2'"], "C3'": ["C3'"], 'C4': ['C4'], "C4'": ["C4'"], 'C5': ['C5'], 
            "C5'": ["C5'"], 'C6': ['C6'], "H1'": ["H1'"], "H2'2": ["H2'2"], "H2'1": ["H2'1"], "H2'": ["H2'1"], "H3'": ["H3'"], 
            "H4'": ["H4'"], 'H41': ['H41'], 'H42': ['H42'], 'H5': ['H5'], "H5'2": ["H5'2"], "H5'1": ["H5'1"], 'H6': ['H6'],
            "H5'": ["H5'1"], 'H5"': ["H5'2"], "H5''": ["H5'2"], "H2''": ["H2'2"], 'H2"': ["H2'2"],
            'N1': ['N1'], 'N3': ['N3'], 'N4': ['N4'], 'O1P': ['O1P'], 'O2': ['O2'], "O2'": ["O2'"], 'O2P': ['O2P'], 
            "O3'": ["O3'"], "O4'": ["O4'"], "O5'": ["O5'"], 'P': ['P']},
    'GLN': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD': ['CD'], 'CG': ['CG'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'],
            'HB3': ['HB3'], 'HE21': ['HE21'], 'HE22': ['HE22'], 'HG2': ['HG2'], 'HG3': ['HG3'], 'N': ['N'], 'NE2': ['NE2'],
            'QB': ['HB2', 'HB3'], 'QE2': ['HE21', 'HE22'], 'QG': ['HG2', 'HG3'], 'H1': ['H']},
    'GLU': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD': ['CD'], 'CG': ['CG'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 
            'HB3': ['HB3'], 'HE2': ['HE2'], 'HG2': ['HG2'], 'HG3': ['HG3'], 'N': ['N'], 'QB': ['HB2', 'HB3'], 'QG': ['HG2', 'HG3'],
            'H1': ['H']},
    'GLY': {'C': ['C'], 'CA': ['CA'], 'H': ['H'], 'HA2': ['HA2'], 'HA3': ['HA3'], 'N': ['N'], 'QA': ['HA2', 'HA3'], 'H1': ['H']},
    'GUA': {"C1'": ["C1'"], 'C2': ['C2'], "C2'": ["C2'"], "C3'": ["C3'"], 'C4': ['C4'], "C4'": ["C4'"], 'C5': ['C5'], 
            "C5'": ["C5'"], 'C6': ['C6'], 'C8': ['C8'], 'H1': ['H1'], "H1'": ["H1'"], "H2'2": ["H2'2"], "H2'1": ["H2'1"], 
            "H2'": ["H2'1"], 'H21': ['H21'], 'H22': ['H22'], "H3'": ["H3'"], "H4'": ["H4'"], "H5'2": ["H5'2"], "H5'1": ["H5'1"],
            "H5'": ["H5'1"], 'H5"': ["H5'2"], "H5''": ["H5'2"], "H2''": ["H2'2"], 'H2"': ["H2'2"],
            'H8': ['H8'], 'N1': ['N1'], 'N2': ['N2'], 'N3': ['N3'], 'N7': ['N7'], 'N9': ['N9'], 'O1P': ['O1P'], "O2'": ["O2'"],
            'O2P': ['O2P'], "O3'": ["O3'"], "O4'": ["O4'"], "O5'": ["O5'"], 'O6': ['O6'], 'P': ['P']},
    'HIS': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD2': ['CD2'], 'CE1': ['CE1'], 'CG': ['CG'], 'H': ['H'], 'HA': ['HA'], 
            'HB2': ['HB2'], 'HB3': ['HB3'], 'HD1': ['HD1'], 'HD2': ['HD2'], 'HE1': ['HE1'], 'HE2': ['HE2'], 'N': ['N'], 
            'ND1': ['ND1'], 'NE2': ['NE2'], 'QB': ['HB2', 'HB3'], 'H1': ['H']},
    'ILE': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD1': ['CD1'], 'CG1': ['CG1'], 'CG2': ['CG2'], 'CQG': ['CG1', 'CG2'], 
            'H': ['H'], 'HA': ['HA'], 'HB': ['HB'], 'HD11': ['HD11'], 'HD12': ['HD12'], 'HD13': ['HD13'], 'HG12': ['HG12'],
            'HG13': ['HG13'], 'HG21': ['HG21'], 'HG22': ['HG22'], 'HG23': ['HG23'], 'MD1': ['HD11', 'HD12', 'HD13'], 'MG2': ['HG21', 'HG22', 'HG23'], 'N': ['N'],
            'QG1': ['HG12', 'HG13'], 'H1': ['H']},
    'LEU': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD1': ['CD1'], 'CD2': ['CD2'], 'CG': ['CG'], 'CQD': ['CD1', 'CD2'], 'H': ['H'], 
            'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'], 'HD11': ['HD11'], 'HD12': ['HD12'], 'HD13': ['HD13'], 'HD21': ['HD21'],
            'HD22': ['HD22'], 'HD23': ['HD23'], 'HG': ['HG'], 'MD1': ['HD11', 'HD12', 'HD13'], 'MD2': ['HD21', 'HD22', 'HD23'], 'N': ['N'], 'QB': ['HB2', 'HB3'], 
            'QMD': ['HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'], 'H1': ['H']},
    'LYS': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD': ['CD'], 'CE': ['CE'], 'CG': ['CG'], 'H': ['H'], 'HA': ['HA'], 
            'HB2': ['HB2'], 'HB3': ['HB3'], 'HD2': ['HD2'], 'HD3': ['HD3'], 'HE2': ['HE2'], 'HE3': ['HE3'], 'HG2': ['HG2'], 
            'HG3': ['HG3'], 'HZ1': ['HZ1'], 'HZ2': ['HZ2'], 'HZ3': ['HZ3'], 'MZ': ['HZ1', 'HZ2', 'HZ3'], 'N': ['N'], 'NZ': ['NZ'], 
            'QB': ['HB2', 'HB3'], 'QD': ['HD1', 'HD2'], 'QE': ['HE1', 'HE2'], 'QG': ['HG2', 'HG3'], 'H1': ['H']},
    'MET': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CE': ['CE'], 'CG': ['CG'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 
            'HB3': ['HB3'], 'HE1': ['HE1'], 'HE2': ['HE2'], 'HE3': ['HE3'], 'HG2': ['HG2'], 'HG3': ['HG3'], 'ME': ['HE1', 'HE2', 'HE3'], 
            'N': ['N'], 'QB': ['HB2', 'HB3'], 'QG': ['HG2', 'HG3'], 'H1': ['H']},
    'PHE': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD1': ['CD1'], 'CD2': ['CD2'], 'CE1': ['CE1'], 'CE2': ['CE2'], 
            'CG': ['CG'], 'CQD': ['CD1', 'CD2'], 'CQE': ['CE1', 'CE2'], 'CZ': ['CZ'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'], 
            'HD1': ['HD1'], 'HD2': ['HD2'], 'HE1': ['HE1'], 'HE2': ['HE2'], 'HZ': ['HZ'], 'N': ['N'], 'QB': ['HB2', 'HB3'], 'QD': ['HD1', 'HD2'], 
            'QE': ['HE1', 'HE2'], 'H1': ['H']},
    'PRO': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD': ['CD'], 'CG': ['CG'], 'H2': ['H2'], 'H3': ['H3'], 'HA': ['HA'], 
            'HB2': ['HB2'], 'HB3': ['HB3'], 'HD2': ['HD2'], 'HD3': ['HD3'], 'HG2': ['HG2'], 'HG3': ['HG3'], 'N': ['N'], 
            'QB': ['HB2', 'HB3'], 'QD': ['HD1', 'HD2'], 'QG': ['HG2', 'HG3'], 'H1': ['H']},
    'SER': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'], 'HG': ['HG'], 
            'N': ['N'], 'QB': ['HB2', 'HB3'], 'H1': ['H']},
    'THR': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CG2': ['CG2'], 'H': ['H'], 'HA': ['HA'], 'HB': ['HB'], 'HG1': ['HG1'], 
            'HG21': ['HG21'], 'HG22': ['HG22'], 'HG23': ['HG23'], 'MG2': ['HG21', 'HG22', 'HG23'], 'N': ['N'], 'H1': ['H']},
    'THY': {"C1'": ["C1'"], 'C2': ['C2'], "C2'": ["C2'"], "C3'": ["C3'"], 'C4': ['C4'], "C4'": ["C4'"], 'C5': ['C5'], 
            "C5'": ["C5'"], 'C6': ['C6'], 'C7': ['C7'], "H1'": ["H1'"], "H2'2": ["H2'2"], "H2'1": ["H2'1"], "H2'": ["H2'1"], 
            'H3': ['H3'], "H3'": ["H3'"], "H4'": ["H4'"], "H5'2": ["H5'2"], "H5'1": ["H5'1"], 'H6': ['H6'], 'H71': ['H71'], 
            'H72': ['H72'], "H5''": ["H5'2"], 'H5"': ["H5'2"], "H2''": ["H2'2"], 'H2"': ["H2'2"], "H5'": ["H5'1"],
            'H73': ['H73'], 'M7': ['H71', 'H72', 'H73'], 'N1': ['N1'], 'N3': ['N3'], 'O1P': ['O1P'], 'O2': ['O2'], 'O2P': ['O2P'], 
            "O3'": ["O3'"], 'O4': ['O4'], "O4'": ["O4'"], "O5'": ["O5'"], 'P': ['P'], 'H1': ['H'], 'M5': ["H5'1", "H5'2"]},
    'TRP': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD1': ['CD1'], 'CD2': ['CD2'], 'CE2': ['CE2'], 'CE3': ['CE3'], 'CG': ['CG'],
            'CH2': ['CH2'], 'CZ2': ['CZ2'], 'CZ3': ['CZ3'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'], 
            'HD1': ['HD1'], 'HE1': ['HE1'], 'HE3': ['HE3'], 'HH2': ['HH2'], 'HZ2': ['HZ2'], 'HZ3': ['HZ3'], 'N': ['N'], 
            'NE1': ['NE1'], 'QB': ['HB2', 'HB3'], 'H1': ['H']},
    'TYR': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CD1': ['CD1'], 'CD2': ['CD2'], 'CE1': ['CE1'], 'CE2': ['CE2'], 'CG': ['CG'], 
            'CQD': ['CD1', 'CD2'], 'CQE': ['CE1', 'CE2'], 'CZ': ['CZ'], 'H': ['H'], 'HA': ['HA'], 'HB2': ['HB2'], 'HB3': ['HB3'], 'HD1': ['HD1'], 
            'HD2': ['HD2'], 'HE1': ['HE1'], 'HE2': ['HE2'], 'HH': ['HH'], 'N': ['N'], 'QB': ['HB2', 'HB3'], 'QD': ['HD1', 'HD2'], 
            'QE': ['HE1', 'HE2'], 'H1': ['H']},
    'URI': {"C1'": ["C1'"], 'C2': ['C2'], "C2'": ["C2'"], "C3'": ["C3'"], 'C4': ['C4'], "C4'": ["C4'"], 'C5': ['C5'], 
            "C5'": ["C5'"], 'C6': ['C6'], "H1'": ["H1'"], "H2'1": ["H2'1"], "H2'": ["H2'1"], 'H3': ['H3'], "H3'": ["H3'"], "H4'": ["H4'"], 
            'H5': ['H5'], "H5'2": ["H5'2"], "H5'1": ["H5'1"], 'H6': ['H6'], 'N1': ['N1'], 'N3': ['N3'], 'O1P': ['O1P'],
            "H5'": ["H5'1"], 'H5"': ["H5'2"], "H5''": ["H5'2"], "H2'2": ["H2'2"], "H2''": ["H2'2"], 'H2"': ["H2'2"],
            "O2'": ["O2'"], 'O2P': ['O2P'], "O3'": ["O3'"], "O4'": ["O4'"], "O5'": ["O5'"], 'P': ['P']},
    'VAL': {'C': ['C'], 'CA': ['CA'], 'CB': ['CB'], 'CG1': ['CG1'], 'CG2': ['CG2'], 'CQG': ['CG1', 'CG2'], 'H': ['H'], 'HA': ['HA'],
            'HB': ['HB'], 'HG11': ['HG11'], 'HG12': ['HG12'], 'HG13': ['HG13'], 'HG21': ['HG21'], 'HG22': ['HG22'], 
            'HG23': ['HG23'], 'MG1': ['HG11', 'HG12', 'HG13'], 'MG2': ['HG21', 'HG22', 'HG23'], 'N': ['N'], 
            'QMG': ['HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23'], 'H1': ['H']},
    
    # Proton 123 prefixes to suffixes
    '1HD2': ['HD21'],
    "1H2'": ["H2'1"],
    '2HH1': ['HH12'],
    '2HD1': ['HD12'],
    '1H2': ['H21'],
    "2H2'": ["H2'2"],
    '1H6': ['H61'],
    '1H7': ['H71'],
    '1H4': ['H41'],
    '3HG1': ['HG13'],
    '3HG2': ['HG23'],
    '2H2': ['H22'],
    '2HD2': ['HD22'],
    '2H4': ['H42'],
    '2H7': ['H72'],
    '2H6': ['H62'],
    '3HA': ['HA3'],
    '3HB': ['HB3'],
    '1HG2': ['HG21'],
    '3HE': ['HE3'],
    '3HG': ['HG3'],
    '3HD': ['HD3'],
    '3HD1': ['HD13'],
    '1HD1': ['HD11'],
    '2HH2': ['HH22'],
    '2HZ': ['HZ2'],
    '3H7': ['H73'],
    '1HZ': ['HZ1'],
    '1HE2': ['HE21'],
    '1HG1': ['HG11'],
    "1H5'": ["H5'1"],
    '1HB': ['HB1'],
    '2HH': ['HH2'],
    '1HA': ['HA1'],
    '1HG': ['HG1'],
    '1HD': ['HD1'],
    '1HE': ['HE1'],
    '2HA': ['HA2'],
    '2HB': ['HB2'],
    '2HE': ['HE2'],
    '2HD': ['HD2'],
    '2HG': ['HG2'],
    '2HG2': ['HG22'],
    '2HG1': ['HG12'],
    '2HE2': ['HE22'],
    "2H5'": ["H5'2"],
    '1HH2': ['HH21'],
    '3HD2': ['HD23'],
    '3HZ': ['HZ3'],
    '1HH1': ['HH11'],
    
    # Custom
    'HX': ["H1'"],
    'X1': ["H1'"],
    'X2': ["H1'"],
  },
}