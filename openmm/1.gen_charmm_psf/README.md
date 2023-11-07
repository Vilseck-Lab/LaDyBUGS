This step requires the use of CHARMM. CHARMM can be obtained following the instructions listed here:
https://academiccharmm.org/program

After CHARMM is installed, ensure that the prep directory contains all system specific files
Then submit the genpsf.job script to slurm (or your local scheduler) to complete this step of setup

Edit nsubs to contain a space deliminted list of the number of substituents to be sampled at each site. If perturbations are made at only one site, then nsubs contains a single number. For example, for the MUP system in `example`, nsubs is '6', representing that 6 substituents are sampled at one site off the ligand core.

Additional perturbation variables are generated based on the contents of nsubs
