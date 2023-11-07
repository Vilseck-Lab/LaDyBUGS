This second step of LaDyBUGS setup creates an OpenMM system file

Edit these lines in setupOMMSystem.py:
 1) boxsize
 2) nameid list (copy from [scriptname])


If phiscale=True in setupOMMSystem.py, take these additional steps to scale specific dihedral angles by lambda:
  - create the dihedrals_to_scale.names.txt file:
        syntax: Substituent# AtomName1 AtomName2 AtomName3 AtomName4 
        ex) "1 C033 C051 C018 C025"      (for scaling the dihedral on substituent 1)
  - run the get_dih_idx.sh script to print a dihedrals_to_scale.txt file (which is read by setupOMMSystem.py


