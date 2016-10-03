# filterPeaks
Making better use of unused peaks

FilterPeaks is designed to resolve peaks that have not been used
for the structure calculation process of a typical protein.
These peaks may contain additional protein assignements but also
protein-ligand NOE peaks. The latter are often not a target of 
structure calculation programs because they focus on resolving 
the protein structure but they contain valuble information on 
the orientation of the ligand.

FilterPeaks performs this task by:
* Loading unassigned peaks from a 3D NOE experiment (Xeasy format).
* Find carbond-proton assignment pairs for unassigned peaks
  matching carbon (13C) chemical shift +/- tolerance and proton (1H)
  chemical shift +/- tolerance in w2 and w3 dimensions respectivly.
