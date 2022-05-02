# Lambda Dynamics with Bias-Updated Gibbs Sampling (LaDyBUGS)
LaDyBUGS is an efficient method of free energy calculation by performing Gibbs sampling over a set of discrete lambda values with a dynamic bias. Herein, we provide up-to-date scripts for performing the process using OpenMM using the CHARMM forcefield.

We provide a step-by-step process for performing a LaDyBUGS calculation. 

# Installation
To clone this repository:
`git clone https://github.com/Vilseck-Lab/LaDyBUGS.git`

You can use anaconda to set up a python environment for LaDyBUGS.

*Install OpenMM*
`conda create --name omm100`
`conda activate omm100`
`conda install -c omnia/label/cuda100 -c conda-forge openmm`

*Install FastMBAR*
In the same environment:
`conda install pytorch==1.2.0 torchvision==0.4.0 cudatoolkit=10.0 -c pytorch`
`pip install FastMBAR`

*Other Requirements*
Some version of CHARMM is needed to generate the psf file. 

# Step 0: Create Hybrid Topology
Use the msld-py-prep repository to create a CHARMM hybrid topology. A tutorial for this process can be found in that repository.
`https://github.com/Vilseck-Lab/msld-py-prep`

Be sure that the files in "toppar/" match the version of cgenff that you use for the parameterization process.

# Step 1: Create Charmm PSF and PDB
Take the build.name directory, and rename it to "prep".
Copy the relevant sections of "prep/name.inp" into "prep/lig_specific.str" (look at the corresponding example file here for the format).
Update the "genpsf.inp" file to load the protein, solvent, and corresponding ions.
Change "box" variable in "genpsf.inp" to match the size of your system.

You can then run "genpsf.inp" through CHARMM to generate the patch files ("patch.pdb", "patch.psf"), an example execution script "genpsf.job" is given.

# Step 2: Create OpenMM patch, generate lambda states
Change the "nsubs" value in "generate_LambdaStates.py" to match the number of substituents you are working with.

Make sure that the "boxsize" variable in "setupOMMSystem.py" matches your system size.

You can execute the setup script with "setupOMMSystem.job". This job execution script will also copy over the requisite files needed to generate the omm patch.

# Step 3: Run LaDyBUGS
Copy "patch.pdb" "patch.psf" "LambdaStates.txt" and "omm_patch.xml" into the "template" directory.
Make sure that in "sampleSystem.py", the "nsubs" variable matches the number of substituents you are working with.

The default sampling parameters sample for 25 ns, with 100 MD steps between Gibbs sampling, 1000 Gibbs sampling steps between FastMBAR updates, and 125 FastMBAR updates overall. This can be changed in "sampleSystem.py".

The default bias equations use a flat bias of 10 kcal/mol for the first set of Gibbs Sampling steps, and an exponential bias afterwards that uses the FastMBAR updates. The mathematical parameters for these functions can also be changed in "sampleSystem.py".

By default, LaDyBUGS takes a snapshot every time an endstate is sampled. These snapshots can be found in the "frames/" directory.

The execution script "sub2q.job" submits an array of replicas that all copy from the "template" directory. 

# Obtaining Final Free Energy Values
When the LaDyBUGS jobs are completed, the configurations from all replicas can be pooled together with a final FastMBAR calculation. The script for running such a calculation is "fastmbar_final.py", and the exection script associated with it is "fastmbar_final.job".

Be sure to adjust the following values in "fastmbar_final.py" to match your run:
- num_lambstates 
- num_endstates 
- temp
- prod_start (should be 0 unless you want to use the first X GS loops for equilibration)

The final energies should be listed in "final_results.txt".
