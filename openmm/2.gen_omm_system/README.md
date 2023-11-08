This second step of LaDyBUGS setup creates an OpenMM system file

Edit these lines in setupOMMSystem.py:
 1) boxsize
 2) nameid list (copy from [scriptname])

Edit generate_LambdaStates.py to correct the nsubs value. Run this script to generate LambdaStates.txt. This file lists out all of the discrete lambda values that will be used in the LaDyBUGS calculation. This list is non unique and can be customized in many different ways.

To scale specific dihedral angles by lambda, set phiscale=True
If phiscale=True in setupOMMSystem.py, take these additional steps to scale specific dihedral angles by lambda:
  - create the dihedrals_to_scale.names.txt file:
        syntax: Substituent# AtomName1 AtomName2 AtomName3 AtomName4 
        ex) "1 C033 C051 C018 C025"      (for scaling the dihedral on substituent 1)
  - run the get_dih_idx.sh script to print a dihedrals_to_scale.txt file (which is read by setupOMMSystem.py


