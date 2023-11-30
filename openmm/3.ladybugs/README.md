## Production LaDyBUGS
After the CHARMM psf, pdb, and OpenMM xml files have been generated, a LaDyBUGS simulation can begin.

The following files should be copied into the template directory before running `sub2q.job`:<br>
- omm_system.xml
- patch.pdb
- patch.psf
- LambdaStates.txt
- phiIDK.txt (if dihedrals will be scaled by lambda)

Edit `sampleSystem.py` and insert the correct value of nsubs. Other file names (e.g., patch.pdb, patch.psf, etc.)
can be changed to match the current system's needs (if they differ from these default names).

Within `sampleSystem.py`, the following variables can be edited to change the sampling behavior of LaDyBUGS: <br>
- numStepPerCycle = 100  (This equals the number of Molecular Dynamics steps performed per Gibbs Sampler step)
- numCycle = 1000  (This equals the number of Gibbs Sampler steps performed between bias updates)
- numBiasLoops = 75  (This equals the number of bias updates performed with FastMBAR (occurs after every numCycle steps))

The total length of LaDyBUGS sampling = numStepPerCycle * numCycle * numBiasLoops * stepsize (in fs). For example, changing numBiasLoops to 125, but leaving numCycle = 1000 and numStepPerCycle = 100, will run a 25 ns simulation.

Finally, submit the `sub2q.job` script to slurm (or your local scheduler) to start a LaDyBUGS simulation. `sub2q.job` currently runs three duplicate calculations via a job array, but this can be customized to match local needs.
