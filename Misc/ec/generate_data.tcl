# Load 15 EC molecules from Packmol
mol new ec15.pdb type pdb

# Load original .mol2 to get correct charges/types
mol addfile ec_custom.mol2 type mol2

# Make sure TopoTools is available
package require topotools

# Assign atom types and charges from mol2
topo guessatom element  ;# or use "name" if you want to preserve O1, O2, etc.

# Write LAMMPS data file
topo writelammpsdata ec15.data full
