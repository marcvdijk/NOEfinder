# NOEfinder
Making better use of unused peaks

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

NOEfinder is designed to resolve peaks that have not been used
for the structure calculation process of a typical protein.
These peaks may contain additional protein assignements but also
protein-ligand NOE peaks. The latter are often not a target of 
structure calculation programs because they focus on resolving 
the protein structure but they contain valuble information on 
the orientation of the ligand.

NOEfinder performs this task by:
* Loading unassigned peaks from a 3D NOE experiment (Xeasy format).
* Find carbond-proton assignment pairs for unassigned peaks
  matching carbon (13C) chemical shift +/- tolerance and proton (1H)
  chemical shift +/- tolerance in w2 and w3 dimensions respectivly.