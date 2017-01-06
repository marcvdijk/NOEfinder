# NOEfinder

Making better use of unused peaks

The NOEfinder program is designed to extract information on intermolecular
contacts involving proteins or proteins and other biomolecules that is 
potentially availale in a collection of unused peaks obtained from a typical
3D 13C-NOESY-HSQC or 15N-NOESY-HSQC NMR experiment.
These peaks may contain additional protein assignments but also 
protein-ligand NOE peaks. The latter are often not a target of structure 
calculation programs because they focus on resolving the protein structure 
but they contain valuable information on the orientation of the ligand.

The program attempts to assign these unused peaks by mimicing parts of the
automatic assignment routines availabe in popular NMR structure calculation
programs using the following execution workflow:

1.  Parse information for unused peaks from a Xeasy peaks file (unused column).
2.  Parse assignment information from a Sparky project file.
3.  Optionally parse the structure of the system under study from a PDB file.
4.  Search for carbon-proton (CH) pairs among the unused peaks by by matching the
    carbon chemical shift with a tolerance of 0.5 ppm and the proton chemical
    shift with a tolerance of 0.05 ppm in the w2 and w3 dimensions respectively.
    Tolerance values can be adjusted.
5.  For the identified CH couple search for a matching unused proton peak in the
    w1 dimension with a proton chemical shift tolerance of 0.05 ppm.
    The identified CHH pair 
6.  If a structure file was parsed, identify the carbon of the carbon-proton
    couple identified in step 4 and the proton of the carbon-proton-proton pair
    identified in step 5. Calculate the distance between the carbon and the
    latter proton in Angstrom.
7.  If the distance is below a predefined cutoff (6A by default) the assignment
    is labeled as potential intermolecular NOE.
8.  Label the CHH pairs as inter- or intra-residue assignment based on the 
    residue numbers. Label the assignment protein-protein (PP) if there are only
    amino-acids involved or protein-other (PO) if one of the partners is not 
    an amino-acid. The requires a correct atom and residue naming scheme to be
    used in the assignment file. NOEfinder includes tools to assist on atom and
    residue translations.
9.  If a structure file was parsed and a distance could be calculated for a CHH
    pair, the information in step 8 can be further specified as through bond (TB),
    through residue (TR) or inter-residue (IR) assignments.
10. Print peak and assignment information for all identified CHH pairs to the 
    terminal standard output in a table format. This table contains the CHH pair
    label assigned in step 8 and 9 and the calculated distance if obtained (999
    otherwise). The table format can be easily parsed and filtered using standard
    UNIX command line tools (e.a. grep,awk,sed).

Input:

- peakfile: H-C-H noesy peakfile in Xeasy format
- Sparky project: a Sparky .proj file
- PDB (obtional): PDB structure file

Run:
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

NOEfinder performs this task by:
* Loading unassigned peaks from a 3D NOE experiment (Xeasy format).
* Find carbond-proton assignment pairs for unassigned peaks
  matching carbon (13C) chemical shift +/- tolerance and proton (1H)
  chemical shift +/- tolerance in w2 and w3 dimensions respectivly.