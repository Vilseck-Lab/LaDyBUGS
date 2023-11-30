# Lambda Dynamics with Bias-Updated Gibbs Sampling (LaDyBUGS)
LaDyBUGS is an efficient alchemical free energy method for computing free energy differences between two or more chemical states. LaDyBUGS uses Gibbs sampling to sample alchemical transformations between many chemical states and uses a dynamic bias to encourage equivalent sampling of all discrete lambda states accompanying these transformations. More information about LaDyBUGS can be found here in the LaDyBUGS article.

In this repository, we provide scripts for performing ligand perturbations with LaDyBUGS using OpenMM and CHARMM force fields. Setup scripts are presented in `1.gen_charmm_psf` and `2.gen_omm_system`. Production scripts can be found in `3.ladybugs`. In `example` directories, we provide example files for performing a LaDyBUGS calculation on a series of MUP1 ligands bound to MUP1. 


# Installation
To clone this repository:
`git clone https://github.com/Vilseck-Lab/LaDyBUGS.git`

You can use anaconda to set up a python environment for LaDyBUGS. The follow utilities are needed to run LaDyBUGS in OpenMM. Please follow each program's installation instructions to separately install all tools within a single python environment.

1) *Install OpenMM* <br>
https://openmm.org/

2) *Install PyTorch and FastMBAR* <br>
https://pytorch.org/ <br>
https://fastmbar.readthedocs.io/en/latest/installation.html

3) *Other Requirements* <br>
If CHARMM force fields will be used, some version of CHARMM is needed to generate the corresponding psf and starting pdb files.  <br>
https://academiccharmm.org/program

# For ligand perturbations, create a hybrid topology
Use the msld-py-prep repository to create a CHARMM hybrid ligand topology. A tutorial for this process can be found in that repository: 
https://github.com/Vilseck-Lab/msld-py-prep

Be sure that the files in `prep/toppar/` match the version of cgenff that you use for ligand CGenFF parameterization.

