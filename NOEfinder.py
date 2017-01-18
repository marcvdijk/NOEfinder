#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import copy
import re
import imp
import logging

from   math import sqrt

__docformat__ = 'restructuredtext'
__version__   = '{major:d}.{minor:d}.{micro:d}'.format(major=1, minor=0, micro=0)
__author__    = 'Marc van Dijk'
__status__    = 'release beta1'
__date__      = '1 January 2017'
__license__   = 'Apache 2.0'
__copyright__ = 'Copyright (c) 2016 Marc van Dijk'
__rootpath__  = os.path.dirname(__file__)

USAGE = """
NOEfinder: Making better use of unused peaks

NOEfinder is designed to resolve peaks that have not been used for the structure
calculation process of a typical protein. These peaks may contain additional
protein assignments but also protein-ligand NOE peaks. The latter are often not
a target of structure calculation programs because they focus on resolving the 
protein structure but they contain valuable information on the orientation of the
ligand.

Input:
- peakfile: H-C-H noesy peakfile in Xeasy format
- Sparky project: a Sparky .proj file
- PDB (obtional): PDB structure file
  
  python NOEfinder.py -p <peakfile> -y <Sparky project> -s <PDB file>

Output:
The NOEfinder program reports results as tabular output written to standard out.
The table reports for each assignment the peak ID and the associated carbon and
two protons from the w2, w3 and w1 dimensions of the 3D experiment respectively.
For each dimension the associated residue name and number, atom name and number
and chemical shift are reported. Each assignment is labeled with a four
character capitalized string (Cat. column). The first two characters indicate:

- IR: Inter Residue assignment, any C-H-H assignment between two different
      residues.
- TB: Through Bond, assignment includes atoms covalently linked (calculated 
      distance < 1.09 Å)
- TR: Through Residue, assignment includes atoms part of the same residue 
      (name and number)
- NX: No carbon atom found in structure.

The second two characters indicate:
- PP: Assignment includes protein residues only (standard amino-acids).
- PO: Assignment includes a protein residue and a non-protein residue judged by
      residue name.
- XX: Not assigned, in case of NC. NC=? = NX?

The label allows for easy sorting and filtering of the output table using
standard UNIX command line tools. The last column of the table holds the
calculated distances in Å between the carbon atom of the identified 
carbon-proton pair and the proton in the w1 dimension in case the PDB structure
file of the associated protein or complex is provided. If such structure file
is not provided or the distance could not be resolved the distance value is 
set to 999.000

Settings:
Basic program settings are defined using the command line interface. 
Advanced options such as residue and atom renaming conventions (see below) and
persistent storage of settings are enabled via the python style configuration
file (config_noefinder.py). Please read the documentation on how to use the
configuration file in the header of that file.

Residue and atom naming:
Identification of residues and atoms in structure files requires their naming
to match those of the assignments which is very often not the case. 
The NOEfinder program offers advanced residue and atom translation services for
this purpose. Please consult the config_noefinder.py file for more information
on how to use it.

Use the -h/--help option for more information on program options
"""

EPILOG = "Version: {version},   Author: {author},   License: {license}".format(
    author=__author__, version=__version__, license=__license__)

# Primary program configuration
STRUCTURE_DICT  = {'id':0, 'atom':None, 'resi':None, 'resn':0, 'coor':(0.00, 0.00, 0.00)}
PEAKS_DICT      = {'id':0, 'w1':0.000, 'w2':0.000, 'w3':0.000, 'spec_type':'U', 'vol':0.000, 'vol_err':0.000, 'intm':'e', 'unused':0, 'ass1':0, 'ass2':0, 'ass3':0}
ASSIGNMENT_DICT = {'id':0, 'shift':0.000, 'nuc':None, 'resi':None, 'resn':0, 'atom':None}
SHIFTCTOL       = 0.5
SHIFTPTOL       = 0.05
DISTTOL         = 6.0
USE_ATOM_TRANS  = False
USE_RESI_TRANS  = False
USE_REGEX_MATCH = False
LOGLEVEL        = 'warning'
CONFIG_FILE     = 'config_noefinder'
NONASSIGNEDID   = set(['NC.', 'NC=?', '= NX?'])

def import_settings(pyfile):
    """
    Import NOEfinder settings from Python style settings file.
    The importer expects to find a dictionary named 'settings' in
    the Python settings file. Only that dictionary is imported.
    
    :param pyfile: name or path to NOEfinder settings file.
    :type pyfile:  str
    :return:       settings
    :rtype:        dict
    """
    
    settings = {}
    
    try:
        fp, pathname, description = imp.find_module(pyfile)
    except ImportError:
        logging.warning('Filter peaks settings file {0} not found'.format(pyfile))
        return settings
    
    try:
        config_filterpeaks = imp.load_module('config_filterpeaks', fp, pathname, description)
        if hasattr(config_filterpeaks, 'settings') and type(config_filterpeaks.settings) == dict:
            settings = config_filterpeaks.settings
    except:
        logging.error('Unable to import Filter peaks settings from {0}'.format(pyfile))
    
    return settings

def translate_residue_names(residue):
    """
    Translate residue name to a IUPAC PDB compatible naming
    convention. By default this means translating one and
    two-letter codes to three letter codes.
    
    Translation is done using the 'residue_translation'
    dictionary in the settings. This dictionary supports
    ambiguity and string lookup using regular expressions.
    
    The function always returns a list one one or multiple
    (ambiguity) translated residues. In case of an unknown
    residue, the original upper cased residue is returned
    as list.
    
    :param residue: residue name to translate
    :type residue:  string
    :return:        translated residue(s)
    :rtype:         list
    """
    
    if not options.use_residue_translation:
        return [residue]
    
    resi = None
    if len(residue):
        resi = options.residue_translation.get(residue.upper(), None)
    
    # Check for regular expression based match if direct lookup failed
    if not resi and options.use_regex_matching:
        for key in options.residue_translation:
            match = re.search(key,residue)
            if match:
                if match.group() == residue:
                    resi = options.residue_translation[key]
                    break
    
    if not resi:
        logging.warning('Unknown residue name {0}'.format(residue))
        return [residue]
    
    if resi != [residue]:
        logging.debug('Translate residue {0} to {1}'.format(residue, ' '.join(resi)))
    
    if type(resi) != list:
        return [resi]
    return resi

def translate_atom_names(atom, residue=None):
    """
    Translate atom names to a IUPAC PDB compatible naming
    convention.

    Translation is done using the 'atom_translation'
    dictionary in the settings. This dictionary supports
    ambiguity and string lookup using regular expressions.

    The function always returns a list one one or multiple
    (ambiguity) translated atoms. In case of an unknown
    atom, the original upper cased atom is returned
    as list.

    :param residue: atom name to translate
    :type residue:  string
    :return:        translated atom(s)
    :rtype:         list
    """
  
    if not options.use_atom_translation:
        return [atom]
  
    atom_name = None
    if len(atom):
        atom = atom.upper()
   
        if residue and residue in options.atom_translation:
            if atom in options.atom_translation[residue]:
                atom_name = options.atom_translation[residue][atom]
            if not atom_name:
                for key in options.atom_translation[residue]:
                    if re.match(key,atom):
                        atom_name = options.atom_translation[residue][key]
                        break
   
        if not atom_name and atom in options.atom_translation:
            atom_name = options.atom_translation[atom]
   
        if not atom_name and options.use_regex_matching:
            for key in options.atom_translation:
                match = re.search(key,atom)
                if match:
                    if match.group() == atom:
                        atom_name = options.atom_translation[key]
                        break
   
        if not atom_name:
            logging.warning('Unknown atom name {0} of residue {1}'.format(atom, residue))
            atom_name = atom
  
    if atom != [atom_name]:
        logging.debug('Translate atom {0} to {1}'.format(atom, ' '.join(atom_name)))
  
    if type(atom_name) != list:
        return [atom_name]
    return atom_name

class _Common(object):
    """
    Base class with methods common to the Assignments,
    Peaks and Structure storage classes.
    """
    
    def __repr__(self):
        """
        Implements the class __repr__ magic method.
        """
        
        obj = type(self).__name__
        return '<{0} object with {1} {2}>'.format(obj.title(),len(self._store),obj)
    
    def __len__(self):
        """
        Implements the class __len__ magic method.
        """
        
        return len(self._store)
    
    def __iter__(self):
        """
        Implements the class __iter__ magic method.
        """
        
        for element in self._store:
            yield element
    
    def __getitem__(self, key):
        """
        Implements the class __getitem__ magic method.
        """
        
        return self._store[key]
    
    def add(self, element):
        """
        Add a new record to the dictionary database.
        Default record fields are defined in the self._base dictionary wich is
        specific to the class that inherits from the _Common base class.
        
        :param element: record to add
        """
        
        newelement = copy.copy(self._base)
        newelement['id'] = element.get('id',self._id)
        newelement.update(element)
        self._store[self._id] = newelement
        self._id += 1
    
    def get(self, key, default=None):
        """
        Emulated the default dictionary get method.
        """
        
        return self._store.get(key, default)
    
    def iter_values(self):
        """
        Iterate over the dictionary values.
        """
        
        for key in sorted(self._store.keys()):
            yield self._store[key]
    
    def iter_items(self):
        """
        Iterate over the dictionary items.
        """
        
        for key in sorted(self._store.keys()):
            yield key, self._store[key]

class Structure(_Common):
    """
    Class that controls a dictionary style storage of PDB style structures.
    Inherits from the _Common baseclass.
    """
    
    def __init__(self):
        
        self._id    = 1
        self._store = {}
        self._base  = STRUCTURE_DICT

class Assignments(_Common):
    """
    Class that controls a dictionary style storage of NMR assignments.
    Inherits from the _Common baseclass.
    """
    
    def __init__(self):
        
        self._id    = 1
        self._store = {}
        self._base  = ASSIGNMENT_DICT

class Peaks(_Common):
    """
    Class that controls a dictionary style storage of NMR peaks.
    Inherits from the _Common baseclass.
    """
    
    def __init__(self):
        
        self._id    = 1
        self._store = {}
        self._base  = PEAKS_DICT
    
    def iter_assigned(self):
        """
        Iterate over all assigned peaks
        
        :return: Assigned peak
        :rtype: dict
        """
        
        for peak_id, peak in self._store.items():
            if self.is_assigned(peak_id):
                yield peak_id, peak
    
    def iter_nonassigned(self):
        """
        Iterate over all non-assigned peaks
        
        :return: None assigned peak
        :rtype: dict
        """
        
        for peak_id, peak in self._store.items():
            if not self.is_assigned(peak_id):
                yield peak_id, peak
    
    def is_assigned(self, peak_id):
        """
        Return True or False if the peak is assigned.
        e.a if all dimensions have an assignment associated.
        
        :param peak_id: peak identifier
        :type peak_id:  int
        """
        
        assert peak_id in self._store, logging.error('No peak with ID: {0}'.format(peak_id))
        return all([self._store[peak_id][ass] != 0 for ass in ('ass1','ass2','ass3')])

def calculatedistance(vector1, vector2):
    """
    Calculate the distance between two vectors (atoms)
    Returns the distance squared or None if vectors not aligned.
    
    :param vector1: First atom coordinate (x,y,z)
    :type vector1:  list
    :param vector2: Second atom coordinate (x,y,z)
    :type vector2:  list
    :return:        squared distance between two atoms
    """
    
    d= float(0)
    if len(vector1) == len(vector2):
        for i in range(len(vector1)):
            d=d+(vector1[i]-vector2[i])**2
        return d
    else:
        return None

def parse_sparky_proj(proj_file):
    """
    Parse Sparky project file (.proj) to Assignments object

    :param proj_file: Sparky project file path
    :type proj_file:  string
    """
  
    assert os.path.exists(proj_file), 'Sparky project file {0} does not exist.'.format(proj_file)
    logging.debug('Parsing Sparky project file {0}'.format(proj_file))
  
    read = False
    assignments = Assignments()
    with open(proj_file) as sparkyproject:
        for line in sparkyproject.readlines():
            if line.startswith('<resonances>'):
                read = True
                continue
            elif line.startswith('<end resonances>'):
                read = False
        
            if read:
                line = line.split()
                
                if not len(line) == 3:
                    logging.warning("Malformed assignment in %s, line: %s" % (proj_file, line))
                    continue
                
                atomid = [a for a in line[0].split('|') if len(a)]
                
                resi = [n for n in re.split('\d', atomid[0]) if len(n)]
                resn = [int(n) for n in re.split('\D', atomid[0]) if len(n)]
                
                if len(resn):
                    resn = resn[0]
                else:
                    resn = None
                
                # Translate residue names
                resi = translate_residue_names("".join(resi))
                if not resi:
                    logging.warning('No valid residue name in line: {0}'.format(line))
                
                for r in resi:
                    
                    # Translate atom names, returns multiple atoms in case of ambiguity
                    atom = translate_atom_names(atomid[1], residue=r)
                    if not len(atom):
                        logging.warning('No valid atom name in line: {0}'.format(line))
                    
                    for a in atom:
                        assignment = {
                          'resi': r,
                          'resn': resn,
                          'atom': a,
                          'shift': float(line[1]),
                          'nuc': line[2],
                        }
                        
                        assignments.add(assignment)
                
    return assignments

def parse_xeasy_peaks(peak_file):
    """
    Parse Xeasy3D peakfile to Peaks object
    
    Xeasy file format stores a column labled 'unused'
    to indicate rather the peak has been used in a
    structure calculation procedure (0 or 1). This
    column is, however, not assigned automatically
    and may not be set at all by the user.
    
    :param peak_file: Xeasy3D peak file path
    :type peak_file:  string
    """
    
    assert os.path.exists(peak_file), 'Xeasy3D peakfile {0} does not exist.'.format(peak_file)
    logging.debug('Parsing Xeasy3D peakfile {0}'.format(peak_file))
    
    peaks = Peaks()
    with open(peak_file) as peakfile:
        for line in peakfile.readlines():
            if not line.startswith('#') and len(line):
                line = line.strip().split()
                if len(line) > 10:
                    
                    peak = {
                        'id':int(line[0]),
                        'w1':float(line[1]),
                        'w2':float(line[2]),
                        'w3':float(line[3]),
                        'spec_type':line[5],
                        'vol':float(line[6]),
                        'vol_err':float(line[7]),
                        'intm':line[8],
                        'unused':int(line[9]),
                        'ass1':int(line[10]),
                        'ass2':int(line[11]),
                        'ass3':int(line[12])
                    }
                    
                    peaks.add(peak)
    
    return peaks

def parse_pdb(pdb_file):
    """
    Parse RCSB PDB file to Structure object
    
    :param pdb_file: PDB structure file path
    :type pdb_file:  str
    """
    
    assert os.path.exists(pdb_file), 'RCSB PDB file {0} does not exist.'.format(pdb_file)
    logging.debug('Parsing RCSB PDB file {0}'.format(pdb_file))
    
    structure = Structure()
    with open(pdb_file) as pdb:
        for line in pdb.readlines():
            if line.startswith('ATOM') and len(line):
                line = line.strip().split()
                
                atom = {
                  'id': int(line[1]),
                  'atom': line[2],
                  'resi': translate_residue_names(line[3]),
                  'resn': int(line[4]),
                  'coor':(float(line[5]), float(line[6]), float(line[7]))
                }
                
                structure.add(atom)
            
            elif line.startswith('END'):
                break
    
    return structure

def get_chh_couple(peaks, assignments, chem_shift_ctol=0.5, chem_shift_ptol=0.05):
    """
    Find carbond-proton assignment pairs for unassigned peaks matching carbon
    chemical shift +/- tolerance and proton chemical shift +/- tolerance in 
    w2 and w3 dimensions respectivly. In addition look for proton assignments 
    in the w1 dimension for previous identified C-H couple.
    
    :param peaks:           peaks
    :type peaks:            peak object
    :param assignments:     assignments
    :type assignments:      assignemnt object
    :param chem_shift_ctol: carbon chemical shift tolerance range in ppm
    :type chem_shift_ctol:  float
    :param chem_shift_ctol: Proton chemical shift tolerance range in ppm
    :type chem_shift_ctol:  float
    :return:                enumerated dictionary with tuples containing
                            in sequence; peak ID and assignment ID's for 
                            C (w2), H (w3) and H (w1).
    :rtype:                 dict
    """
    
    logging.info('Looking for Carbon-Proton-Proton pairs. C chemical shift tolerance: {0:.3f}, H chemical shift tolerance: {1:.3f}'.format(chem_shift_ctol, chem_shift_ptol))
     
    # For all unassigned peaks, search assignments for residues
    # that match the carbon shift +/- chem_shift_ctol
    chh_couple = {}
    enum = 1
    for pid,peak in peaks.iter_nonassigned():
        
        for id2,ass2 in assignments.iter_items():
            if not (peak['w2'] > (ass2['shift'] - chem_shift_ctol) and
                    peak['w2'] < (ass2['shift'] + chem_shift_ctol) and
                    ass2['nuc'] == '13C'):
                continue
            
            # Carbon found, continue to look for proton
            for id3,ass3 in assignments.iter_items():
                if (peak['w3'] > (ass3['shift'] - chem_shift_ptol) and
                    peak['w3'] < (ass3['shift'] + chem_shift_ptol) and
                    ass2['resn'] == ass3['resn'] and
                    ass3['nuc'] == '1H'):
                    
                    chh_couple[enum] = (peak['id'], id2, id3)
                    logging.debug('Peak {0}: carbon-proton assignment pair found {1}-{2}'.format(peak['id'],ass2['id'],ass3['id']))
                    
                    # Continue to look for unassigned peaks in w1 dimension
                    w1_found = False
                    for id1,ass1 in assignments.iter_items():
                        if (peak['w1'] > (ass1['shift'] - chem_shift_ptol) and
                            peak['w1'] < (ass1['shift'] + chem_shift_ptol) and
                            ass1['nuc'] == '1H'):
                    
                            chh_couple[enum] = (peak['id'], id2, id3, id1)
                            logging.debug('Peak {0}: found proton assignment in w1: {1} for C-H couple: {2}-{3}'.format(peak['id'],ass1['id'],ass2['id'],ass3['id']))
                            enum += 1
                            w1_found = True
                    
                    if not w1_found:
                        enum += 1
                            
            if not peak['id'] in chh_couple:
                logging.info('No C-H couple found for peak {0}'.format(peak['id']))
    
    logging.info('Identified {0} C-H-(H) assignment couples for unassigned peaks'.format(len(chh_couple)))
    return chh_couple

def find_unused(chh_couple, peaks, assignments, structure=None, dist_cutoff=6.0, chem_shift_ctol=0.5, chem_shift_ptol=0.05):
    """
    For each previously identified Carbon-Proton-(Proton) (w2-w3-(w1)) pair do:
    - Check if the carbon coordinates can be retrieved from the structure if
      available.
    - If carbon coordinates, check if the attached proton (w3) can be found in 
      structure and validate the covalent bond length.
    - If w3 proton coordinates, check if the assignment in the w1 dimension 
      (the unassigned peaks) can be found in the structure.
    - If coordinates where found, check if the distance is below dist_cutoff
      then store the distance.
    - Catagorize the carbon-proton-proton pair using the 'catagorize_assignment'
      and print the results to standard out.
    
    :param chh_couple:      C-H-(H) (W2-W3-(W1)) couples as identified by 
                            get_chh_couple function.
    :type chh_couple:       dict
    :param peaks:           peaks object
    :param assignments:     assignments object
    :param structure:       structure object
    :param dist_cutoff:     distance cutoff in Angstrom
    :type dist_cutoff:      float
    :param chem_shift_ctol: carbon chemical shift tolerance range in ppm
    :type chem_shift_ctol:  float
    :param chem_shift_ctol: Proton chemical shift tolerance range in ppm
    :type chem_shift_ctol:  float
    """
    
    logging.info('Looking for Carbon-Proton-Proton pairs. C chemical shift tolerance: {0:.3f}, H chemical shift tolerance: {1:.3f}, Distance cutoff: {2:.3f}'.format(chem_shift_ctol, chem_shift_ptol, dist_cutoff))
    
    dist_cutoff = dist_cutoff**2
    print("      <Carbon dimension (w2)>         <Proton dimension (w3)           <Proton dimension (w1)>")
    print("Peak  ID resi resn atom shift         ID resi resn atom shift          ID resi resn atom shift           Cat.  Distance (A)")
    
    # Using peaks for which CH(H) couple was identified, validate distances
    for pid in sorted(chh_couple.keys()):
        
        peak       = peaks[chh_couple[pid][0]]
        assign_cw2 = assignments[chh_couple[pid][1]]
        assign_hw3 = assignments[chh_couple[pid][2]]
        assign_hw1 = assignments[chh_couple[pid][3]] or None
        
        # Look for carbon atom in structure and store coordinates.
        dist = 999
        if structure:
            carbon = None
            for aid_cw2,atom_cw2 in structure.iter_items():
                if atom_cw2['resn'] == assign_cw2['resn'] and atom_cw2['atom'] == assign_cw2['atom']:
                    carbon = atom_cw2['coor']
                    break
        
            # If carbon found, continue
            if carbon and assign_hw1:
                
                # Find proton W3 and calculate distance to carbon to ensure it is
                # covalently attached.
                w3proton = None
                for aid_w3,atom_w3 in structure.iter_items():
                    if atom_w3['resn'] == assign_hw3['resn'] and atom_w3['atom'] == assign_hw3['atom']:
                        cw2_w3_dist = calculatedistance(carbon, atom_w3['coor'])
                        if 1.16 < cw2_w3_dist < 1.21:
                            w3proton = atom_w3['coor']
                
                if not w3proton:
                    logging.debug('Proton {0} {1}-{2} of carbon-proton pair {3}-{4} (peak {5}) not found in structure'.format(assign_hw3['atom'],
                        assign_hw3['resi'], assign_hw3['resn'], assign_hw3['id'], assign_hw3['id'], pid))
                    continue
                
                # Calculate distance between proton W3 and proton W1.
                # There could be ambiquity in atom numbering/naming resulting in multiple 
                # distances calculated (but unlikely). Report the smallest distance
                found_distance = []
                for aid_w1,atom_w1 in structure.iter_items():
                    if atom_w1['resn'] == assign_hw1['resn'] and atom_w1['atom'] == assign_hw1['atom']:
                        dist = calculatedistance(w3proton, atom_w1['coor'])
                        if dist <= dist_cutoff:
                            found_distance.append(dist)
                
                dist = 999
                if found_distance:
                    if len(found_distance) > 1:
                        logging.warning('Multiple atomic distances calculated for CH{0}-H{1}. Report smallest distance'.format(atom_cw2['atom'],atom_w1['atom']))
                    dist = sqrt(min(found_distance))
            
                cat = catagorize_assignment(assign_cw2, assign_hw3, assign_hw1, dist)
                print_chh_info(assign_cw2, assign_hw3, assign_hw1, pid, cat, dist)
                
                continue
                   
            logging.warning('Carbon atom {0} {1}-{2} of carbon-proton pair {3}-{4} (peak {5}) not found in structure'.format(assign_cw2['atom'],
                assign_cw2['resi'], assign_cw2['resn'], assign_cw2['id'], assign_hw3['id'], pid))
        
        # If no structure or no carbon print C-H-(H) couple
        cat = catagorize_assignment(assign_cw2, assign_hw3, assign_hw1, dist)
        print_chh_info(assign_cw2, assign_hw3, assign_hw1, pid, cat, dist)
        
def catagorize_assignment(c, ch, h, dist):
    """
    Catagorize the identified C-H-(H) assigned couple by labeling them with a
    four character capitalized string (Cat. column). The first two characters
    indicate:

    - IR: Inter Residue assignment, any C-H-H assignment between two different
          residues.
    - TB: Through Bond, assignment includes atoms covalently linked (calculated 
          distance < 1.09 Å)
    - TR: Through Residue, assignment includes atoms part of the same residue 
          (name and number)
    - NX: No carbon atom found in structure, no distance calculated

    The second two characters indicate:
    - PP: Assignment includes protein residues only (standard amino-acids).
    - PO: Assignment includes a protein residue and a non-protein residue judged by
          residue name.
    - XX: Not assigned, in case of NC. NC=? = NX? (conatined in NONASSIGNEDID set)
    
    :param c:    carbon W2 assignment
    :type c:     assignments object
    :param ch:   proton W3 assignment
    :type ch:    assignments object
    :param h:    proton W1 assignment
    :type h:     assignments object
    :param dist: calculated distance or 999
    :type dist:  int
    :return:     assignment catagory
    :rtype:      str
    """
    
    resi = [c['resi'],ch['resi'],h['resi']]
    resn = [c['resn'],ch['resn'],h['resn']]
    
    cat = 'IR'
    if dist < 1.09:
        cat = 'TB'
    elif dist == 999:
        cat = 'NX'
    elif len(set(resi)) == 1 and len(set(resn)) == 1:
        cat = 'TR'
    
    if all([r in options.protein_residues for r in resi]):
        cat += 'PP'
    elif len(NONASSIGNEDID.intersection(set(resi))):
        cat += 'XX'
    else:
        cat += 'PO'
    
    return cat

def print_chh_info(assign_cw2, assign_hw3, assign_hw1, pid, cat, dist):
    """
    Report C-H-(H) couple assignment information in tabular form
    """
    
    if assign_hw1:
        print("%-4i C: %4i %3s %-3i %4s %-10.3f C-H: %4i %3s %-3i %4s %-10.3f H: %4i %3s %-4i %5s %3.3f %9s %5.3f" % (pid,
            assign_cw2['id'], assign_cw2['resi'], assign_cw2['resn'], assign_cw2['atom'], assign_cw2['shift'],
            assign_hw3['id'], assign_hw3['resi'], assign_hw3['resn'], assign_hw3['atom'], assign_hw3['shift'], 
            assign_hw1['id'], assign_hw1['resi'], assign_hw1['resn'], assign_hw1['atom'], assign_hw1['shift'], cat, dist))
    else:
        print("%-4i C: %4i %3s %-3i %4s %-10.3f C-H: %4i %3s %-3i %4s %-10.3f H:  not found in structure   %9s %5.3f" % (pid,
            assign_cw2['id'], assign_cw2['resi'], assign_cw2['resn'], assign_cw2['atom'], assign_cw2['shift'],
            assign_hw3['id'], assign_hw3['resi'], assign_hw3['resn'], assign_hw3['atom'], assign_hw3['shift'], cat, dist))
    
def main(options):
    """
    Main program.
    - Parse peaks, assignments and structure files
    - Look for carbon-proton pairs in the assignments
    - Look for missing proton and check distance in structure
    - Print the results
    
    :param options: argparse command line options
    :type options:  argparse argument object
    """
    
    # Parse Xeasy peakfile into peaks class
    peaks = parse_xeasy_peaks(options.peak_file)
    
    # Parse Sparky project file
    assignments = parse_sparky_proj(options.sparky_file)
    
    # Parse PDB structure file into structure class
    structure = None
    if options.pdb_file:
        structure = parse_pdb(options.pdb_file)
    
    # Look for carbon-proton pairs
    chh_couple = get_chh_couple(peaks, assignments,
                            chem_shift_ctol=options.chem_shift_ctol,
                            chem_shift_ptol=options.chem_shift_ptol)
    
    # Look for unused peaks in w1 dimension that could be a missing proton in the ch pair
    find_unused(chh_couple, peaks, assignments,
              structure=structure,
              dist_cutoff=options.struc_dist_tol,
              chem_shift_ctol=options.chem_shift_ctol,
              chem_shift_ptol=options.chem_shift_ptol)

if __name__ == '__main__':
    
    # Import standard configuration file from current directory if any
    default_options = {}
    default_options.update(import_settings(CONFIG_FILE))
    
    # Running from command line, parse argument
    import argparse
    
    parser = argparse.ArgumentParser(
        description='filterPeaks command line options', usage=USAGE, epilog=EPILOG,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( "-p", "--peaks", action="store", dest="peak_file", type=str, help="Xeasy peak file", required=True)
    parser.add_argument( "-a", "--assignemnts", action="store", dest="assignment_file", type=str, help="Xeasy assignement file")
    parser.add_argument( "-y", "--sparky", action="store", dest="sparky_file", type=str, help="Sparky .proj project file", required=True)
    parser.add_argument( "-s", "--pdb", action="store", dest="pdb_file", type=str, help="PDB structure file")
    parser.add_argument( "-c", "--shiftctol", action="store", dest="chem_shift_ctol", type=float, default=default_options.get('chem_shift_ctol',SHIFTCTOL), help="Carbon chemical shift tolerance in PPM")
    parser.add_argument( "-b", "--shiftptol", action="store", dest="chem_shift_ptol", type=float, default=default_options.get('chem_shift_ptol',SHIFTPTOL), help="Proton chemical shift tolerance in PPM")
    parser.add_argument( "-d", "--disttol", action="store", dest="struc_dist_tol", type=float, default=default_options.get('struc_dist_tol',DISTTOL), help="Distance cutoff in A")
    parser.add_argument( "-l", "--loglevel", action="store", dest="loglevel", choices=['critical','error','warning','info','debug'], type=str, default=default_options.get('loglevel',LOGLEVEL), help="Verbosity level")
    parser.add_argument( "-i", "--config", action="store", dest="config_file", type=str, help="FilterPeaks Python configuration file")
    parser.add_argument( "-t", "--use_atom_translation", action="store_true", dest="use_atom_translation", default=default_options.get('use_atom_translation',USE_ATOM_TRANS), help="Translate atom names according to translation rules in config_filterpeaks.")
    parser.add_argument( "-e", "--use_residue_translation", action="store_true", dest="use_residue_translation", default=default_options.get('use_residue_translation',USE_RESI_TRANS), help="Translate residue names according to translation rules in config_filterpeaks.")
    parser.add_argument( "-r", "--use_regex_matching", action="store_true", dest="use_regex_matching", default=default_options.get('use_regex_matching',USE_REGEX_MATCH), help="Use regular expression matching in residue and atom translation.")
    
    # Add default_options to options and make it global
    global options
    options = parser.parse_args()
    for k,v in default_options.items():
        if k not in options:
            options.__setattr__(k,v)
    
    # If config file defined. Parse and update options
    if options.config_file:
        config_filterpeaks = import_settings(options.config_file)
        for k,v in config_filterpeaks.items():
            options.__setattr__(k,v)
    
    # Configure base logger
    loglevel = getattr(logging, options.loglevel.upper())
    logging.basicConfig(level=loglevel, format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
    
    # Call command line main function
    main(options)
    sys.exit(0)