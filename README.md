# Lambda Dynamics with Bias-Updated Gibbs Sampling (LaDyBUGS)
LaDyBUGS is an efficient alchemical free energy method for computing free energy differences between two or more chemical states. LaDyBUGS uses Gibbs sampling to sample alchemical transformations between chemical states and uses a dynamic bias to encourage the equivalent sampling of all discrete lambda values accompanying these transformations. More information about LaDyBUGS can be found here: [LINK DOI].

In this repository, we provide up-to-date scripts for performing ligand perturbations with LaDyBUGS using OpenMM or pyCHARMM with the CHARMM forcefield. In the `example` directories, we provide a step-by-step tutorial for performing a LaDyBUGS ligand calculation. 


# Installation
To clone this repository:
`git clone https://github.com/Vilseck-Lab/LaDyBUGS.git`

You can use anaconda to set up a python environment for LaDyBUGS.

1) *Install OpenMM* <br>
https://openmm.org/

2) *Install PyTorch and FastMBAR* <br>
https://pytorch.org/ <br>
https://fastmbar.readthedocs.io/en/latest/installation.html

3) *Other Requirements* <br>
Some version of CHARMM is needed to generate the psf and starting pdb files.  <br>
https://academiccharmm.org/program

# For ligand perturbations, create a hybrid topology
Use the msld-py-prep repository to create a CHARMM hybrid topology. A tutorial for this process can be found in that repository. <br>
`https://github.com/Vilseck-Lab/msld-py-prep`

Be sure that the files in "toppar/" match the version of cgenff that you use for the parameterization process.

