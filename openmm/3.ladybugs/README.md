After the charmm psf and openmm xml files have been generated, a LaDyBUGS simulation can begin.

The following files should be copied into the template directory before running sub2q.job:
    omm_system.xml, patch.pdb, patch.psf, LambdaStates.txt, phiIDK.txt

Edit sampleSystem.py and insert the correct value of nsubs. Other file names (e.g., patch.pdb, patch.psf, etc.)
can be changed to match the current system's needs.

Edit the following variables to change the sampling behavior of LaDyBUGS:
numStepPerCycle = 100      ## number of Molecular Dynamics steps per Gibbs Sampler step
numCycle = 1000            ## number of Gibbs Sampler steps between bias updates
numBiasLoops = 75          ## number of bias updates with FastMBAR (occurs every numCycle steps)
                           ## total length of sampling = numStepPerCycle * numCycle * numBiasLoops * stepsize (in fs)

For example, changing numBiasLoops to 125 will run a 25 ns simulation

Submit the sub2q.job script to slurm to start a LaDyBUGS simulation
