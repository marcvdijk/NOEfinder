#!/usr/bin/env python

import os, sys, copy, re, imp
import logging
import ConfigParser

from   math      import sqrt

__docformat__ = 'restructuredtext'
__version__   = '{major:d}.{minor:d}.{micro:d}'.format(major=1, minor=0, micro=0)
__author__    = 'Marc van Dijk'
__status__    = 'release beta1'
__date__      = '10 March 2016'
__license__   = "Apache 2.0"
__rootpath__  = os.path.dirname(__file__)

USAGE = """
  Making better use of unused peaks (workingtitle)

Author:     {author}
Version:    {version}
License:    {license}

Required input files:
- peakfile: H-C-H noesy peakfile in Xeasy format
- Sparky project: a Sparky .proj file
- PDB: PDB structure file

  ./NOEfinder.py -p <peakfile> -y <Sparky project> -s <PDB file>
  
Use the -h/--help option for more information on program options

FilterPeaks is designed to resolve peaks that have not been used
for the structure calculation process of a typical protein.
These peaks may contain additional protein assignements but also
protein-ligand NOE peaks. The latter are often not a target of 
structure calculation programs because they focus on resolving 
the protein structure but they contain valuble information on 
the orientation of the ligand.

FilterPeaks performs this task by:
- Loading unassigned peaks from a 3D NOE experiment (Xeasy format).
- Find carbond-proton assignment pairs for unassigned peaks
  matching carbon (13C) chemical shift +/- tolerance and proton (1H)
  chemical shift +/- tolerance in w2 and w3 dimensions respectivly.
- 

  # The first unassigned peak is the number 2 with
  # 1.626  60.334   3.163 ppm
  # The first ppm refer to unknown, 60.334 and 3.16 refer to a CH of protein.
  # I'm gone on that peak using sparky and check possible assignments.
  # Sparky suggests Lys106 Hd* or Lys135Hd*, protein peaks, or thy1ch3/thy10ch3/thy38ch3/thy9ch3, DNA peaks, for 1.62 ppm.
  # 60.334 ppm/3.163ppm belong to Thr114 Ca-Ha.
  # The distance Thr114Ha-Lys106Hd* and Thr114Ha-Lys135Hd* is too large then it could be a noe between
  # thr114Ha and thy1ch3/thy10ch3/thy38ch3/thy9ch3

""".format(author=__author__, version=__version__, license=__license__)

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
CONFIG_FILE     = 'config_filterpeaks'

def import_settings(pyfile):
  
  """
  Import NOEfinder settings from Python style settings file.
  The importer expects to find a dictionary named 'settings' in
  the Python settings file. Only that dictionary is imported.
  
  :param pyfile: Name or path to NOEfinder settings file.
  :return: settings
  :rtype: dict
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
  :ptype residue: string
  :return: translated residue(s)
  :rtype: list
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
          print key, residue, match.group()
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
  :ptype residue: string
  :return: translated atom(s)
  :rtype: list
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
    Default record fields are defined in the self._base
    dictionary wich is specific to the class that 
    inherits from the _Common base class.
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
  Class that controls a dictionary style 
  storage of PDB style structures. 
  Inherits from the _Common baseclass.
  """
  
  def __init__(self):
    
    self._id    = 1
    self._store = {}
    self._base  = STRUCTURE_DICT  
      
class Assignments(_Common):
  
  """
  Class that controls a dictionary style 
  storage of NMR assignments. 
  Inherits from the _Common baseclass.
  """
  
  def __init__(self):
    
    self._id    = 1
    self._store = {}
    self._base  = ASSIGNMENT_DICT
     
class Peaks(_Common):
  
  """
  Class that controls a dictionary style 
  storage of NMR peaks. 
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
    :ptype peak_id: int
    """
    
    assert peak_id in self._store, logging.error('No peak with ID: {0}'.format(peak_id))    
    return all([self._store[peak_id][ass] != 0 for ass in ('ass1','ass2','ass3')])
          
def calculatedistance(vector1, vector2):
  
  """
  Calculate the distance between two vectors (atoms)
  Returns the distance squared or None if vectors not aligned.
  
  :param vector1: First atom coordinate (x,y,z)
  :ptype vector1: list
  :param vector2: Second atom coordinate (x,y,z)
  :ptype vector2: list
  :return: squared distance between two atoms
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
  :ptype proj_file: string
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
  :ptype peak_file: string
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
  :ptype pdb_file: string
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

def get_ch_couple(peaks, assignments, chem_shift_ctol=0.5, chem_shift_ptol=0.05):
  
  """
  Find carbond-proton assignment pairs for unassigned peaks
  matching carbon chemical shift +/- tolerance and proton 
  chemical shift +/- tolerance in w2 and w3 dimensions respectivly.
  
  :param peaks: peaks
  :ptype peaks: peak object
  :param assignments: assignments
  :ptype assignments: assignemnt object
  :param chem_shift_ctol: carbon chemical shift tolerance range in ppm
  :ptype chem_shift_ctol: float
  :param chem_shift_ctol: Proton chemical shift tolerance range in ppm
  :ptype chem_shift_ctol: float
  :return: pairs of carbon,proton assignment ID's stored as dictionary
           of peak id,assignment pair key/value pairs
  :rtype: dict
  """
  
  logging.info('Looking for Carbon-Proton pairs. C chemical shift tolerance: {0:.3f}, H chemical shift tolerance: {1:.3f}'.format(chem_shift_ctol, chem_shift_ptol))
  
  # For all unassigned peaks, search assignments for residues
  # that match the carbon shift +/- chem_shift_ctol
  ch_couple = {}
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
          
          ch_couple[peak['id']] = (id2, id3)
          logging.debug('Peak {0}: carbon-proton assignment pair found {1}-{2}'.format(peak['id'],ass2['id'],ass3['id']))
    
      if not peak['id'] in ch_couple:
        logging.info("No C-H couple found for peak %i" % peak['id'])
  
  return ch_couple

def find_unused(ch_couple, peaks, assignments, structure, dist_cutoff=6.0, chem_shift_ctol=0.5, chem_shift_ptol=0.05):
  
  """
  For each previously identified Carbon-Proton pair do:
  - Check if the carbon coordinates can be retrieved from the structure.
  - If the carbon coordinates where not found, protons in the w1 dimension cannot be resolved. 
    Store the C-H couple as 'noccoor' and continue.
  - If the coordinates where found, check the w1 dimension for unassigned peaks matching carbon
    chemical shift +/- tolerance.
  - If a peak was found, check if the coordinates can be retrieved from the structure.
  - If coordinates where found, check if the distance is below dist_cutoff and store with the 
    C-H couple in 'unused'.
  - If the coordinates where not found, store with the C-H couple in 'distcutoff'. 
  """
  
  logging.info('Looking for Carbon-Proton-Proton pairs. C chemical shift tolerance: {0:.3f}, H chemical shift tolerance: {1:.3f}, Distance cutoff: {2:.3f}'.format(chem_shift_ctol, chem_shift_ptol, dist_cutoff))
  
  dist_cutoff = dist_cutoff**2
  print("      <Carbon dimension (w2)>         <Proton dimension (w3)           <Proton dimension (w1)>")
  print("Peak  ID resi resn atom shift         ID resi resn atom shift          ID resi resn atom shift           Cat.  Distance (A)")
  
  # Using peaks for which CH couple was identified, scan for proton shifts of protein/ligand/noise (w1)
  for pid in sorted(ch_couple.keys()):
    
    peak     = peaks[pid]
    assign_c = assignments[ch_couple[pid][0]]
    assign_h = assignments[ch_couple[pid][1]]
    
    # Look for carbon atom in structure and store coordinates.
    carbon = None
    for aid,atom in structure.iter_items():
      if atom['resn'] == assign_c['resn'] and atom['atom'] == assign_c['atom']:
        carbon = atom['coor']
        break
    
    # If carbon found, continue
    if carbon:
      
      # Look for unused peaks
      for id1,ass1 in assignments.iter_items(): 
        if (peak['w1'] > (ass1['shift'] - chem_shift_ptol) and
            peak['w1'] < (ass1['shift'] + chem_shift_ptol) and
            ass1['nuc'] == '1H'):
          
          found_distance = []
          for aid1,atom1 in structure.iter_items():
            if atom1['resn'] == ass1['resn'] and atom1['atom'] == ass1['atom']:
              dist = calculatedistance(carbon, atom1['coor'])
              if dist <= dist_cutoff:
                found_distance.append(dist)
                
          # No distance found
          dist = 999
          if found_distance:
            dist = sqrt(min(found_distance))
          
          cat = catagorize_assignment(assign_c, assign_h, ass1, dist)
          print("%-4i C: %4i %3s %-3i %4s %-10.3f C-H: %4i %3s %-3i %4s %-10.3f H: %4i %3s %-4i %5s %3.3f %9s %5.3f" % (pid, 
            assign_c['id'], assign_c['resi'], assign_c['resn'], assign_c['atom'], assign_c['shift'],
            assign_h['id'], assign_h['resi'], assign_h['resn'], assign_h['atom'], assign_h['shift'],
            ass1['id'], ass1['resi'], ass1['resn'], ass1['atom'], ass1['shift'], cat, dist))
               
    else:     
      logging.warning('Carbon atom {0} {1}-{2} of carbon-proton pair {3}-{4} (peak {5}) not found in structure'.format(assign_c['atom'],
        assign_c['resi'], assign_c['resn'], assign_c['id'], assign_h['id'], pid))
      
      print("%-4i C: %4i %3s %-3i %4s %-10.3f C-H: %4i %3s %-3i %4s %-10.3f H:  not found in structure   %9s" % (pid, 
        assign_c['id'], assign_c['resi'], assign_c['resn'], assign_c['atom'], assign_c['shift'],
        assign_h['id'], assign_h['resi'], assign_h['resn'], assign_h['atom'], assign_h['shift'], 'NCXX'))

def catagorize_assignment(c, ch, h, dist):
  
  resi = [c['resi'],ch['resi'],h['resi']]
  resn = [c['resn'],ch['resn'],h['resn']]
  
  cat = 'IR'
  if dist < 1.09:
    cat = 'TB'
  elif len(set(resi)) == 1 and len(set(resn)) == 1:
    cat = 'TR'
  
  if all([r in options.protein_residues for r in resi]):
    cat += 'PP'
  else:
    cat += 'PO'
  
  return cat
    
def main(options):
  
  """
  Main program.
  - Parse peaks, assignments and structure files
  - Look for carbon-proton pairs in the assignments
  - Look for missing proton and check distance in structure
  - Print the results
  """
  
  # Parse Xeasy peakfile into peaks class
  peaks = parse_xeasy_peaks(options.peak_file)
    
  # Parse Sparky project file
  assignments = parse_sparky_proj(options.sparky_file)
   
  # Parse PDB structure file into structure class
  structure = parse_pdb(options.pdb_file)
  
  # Look for carbon-proton pairs
  ch_couple = get_ch_couple(peaks, assignments, 
                            chem_shift_ctol=options.chem_shift_ctol,
                            chem_shift_ptol=options.chem_shift_ptol)
  
  # Look for unused peaks in w1 dimension that could be a missing proton in the ch pair
  find_unused(ch_couple, peaks, assignments, structure, 
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
      description='filterPeaks command line options', usage=USAGE,
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument( "-p", "--peaks", action="store", dest="peak_file", type=str, help="Xeasy peak file", required=True)
  parser.add_argument( "-a", "--assignemnts", action="store", dest="assignment_file", type=str, help="Xeasy assignement file")
  parser.add_argument( "-y", "--sparky", action="store", dest="sparky_file", type=str, help="Sparky .proj project file", required=True)
  parser.add_argument( "-s", "--pdb", action="store", dest="pdb_file", type=str, help="PDB structure file", required=True)
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