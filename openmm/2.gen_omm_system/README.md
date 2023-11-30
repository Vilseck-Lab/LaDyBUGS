## Generate an OpenMM system file
This second step of LaDyBUGS setup creates an OpenMM system file (`omm_system.xml`).

The `setupOMMSystem.job` file copies in the necessary psf and pdb files from `../1.gen_charmm_psf` and links to the CHARMM prep directory. Prior to submitting `setupOMMSystem.job`, ensure that it loads your correct python environment and features the correct queuing keywords.

We then need to edit these lines in setupOMMSystem.py:
 1) line 10: boxsize (insert the size of the solvated box)
 2) line 66: nameid list (copy the 0-indexed atom number output from get_alchem_idx.py in `../1.gen_charmm_psf`. This is a list of appended values)

Edit `generate_LambdaStates.py` and give it the correct nsubs value. Run this script to generate `LambdaStates.txt`. This file lists out all of the discrete lambda states that will be sampled in the LaDyBUGS calculation. This list is not unique and can be customized in many different ways.

To scale specific dihedral angles by lambda, set `phiscale=True` in `setupOMMSystem.py` <br>
If `phiscale=True` in `setupOMMSystem.py`, take these additional steps to scale specific dihedral angles by lambda:
  - write out the `dihedrals_to_scale.names.txt` file to list dihedral angles that should be scaled and which substituent they belong to
        syntax: "Substituent# AtomName1 AtomName2 AtomName3 AtomName4"
        ex) "1 C033 C051 C018 C025"      (for scaling a dihedral on substituent 1)
  - run the `get_dih_idx.sh` script to print a `dihedrals_to_scale.txt` file. This file will be read by `setupOMMSystem.py`.

Finally, submit the `setupOMMSystem.job` script to slurm (or your local scheduler) to generate an OpenMM system file (`omm_system.xml`).


