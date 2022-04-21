#! /usr/bin/env python

# Non-Equilibrium Gibbs Sampler Lambda Dynamics + MBAR in OpenMM

from simtk.openmm.app import Simulation, CharmmPsfFile, PDBFile, StateDataReporter
from simtk.openmm import XmlSerializer, MonteCarloBarostat, LangevinIntegrator, Platform
from simtk.unit import kelvin, femtoseconds, picoseconds, bar, kilocalorie_per_mole
from sys import stdout, exit, stderr
from functions import cleanMEM, sampleLambda, update_exp_biases, calcFastMBAR
import numpy as np

## requisite files
xml_file = './omm_patch.xml'
pdb_file = './patch.pdb'
psf_file = './patch.psf'
lstate_file = './LambdaStates.txt'

## MD parameters
kb = 1.987204259 * 0.001 # kcal/(mol*K)
temp = 298.15 # temp in kelvin
pressure = 1.01325 #pressure in bar
stepsize = 2 * femtoseconds

## negs parameters
numStepPerCycle = 100    ## number of Molecular Dynamics steps per numCycle run
numCycle = 1000         ## number of Gibb Sampler Steps (500 = 500ps when numStepPerCycle is 500)
numBiasLoops = 25         ## number of MBAR loops to complete

#exponential bias function variables # exp_bias = exp_mult * exp_base ** (counts[i] - lowest_count)
pre_exp_bias = 100
exp_base = 2.0     ## number to be exponentiated
exp_mult = 1
exp_start = 1      # the first loop (0-ordered) where exp_biasing is used  


## read system
f = open(xml_file,'r')
xml = f.read()
f.close()
system = XmlSerializer.deserialize(xml)

## force group you want energies for
fgE=20

## parameters
T = temp * kelvin
beta = (temp * kb) ** -1.0
#stepsize = 2 * femtoseconds

## MonteCarloBarostat
mcbar = MonteCarloBarostat(pressure*bar,T,25)
system.addForce(mcbar)

## read coordinate from pdb file
psf = CharmmPsfFile(psf_file)
pdb = PDBFile(pdb_file)

## build OMM context object
fricCoef = 10/picoseconds
integrator = LangevinIntegrator(T, fricCoef, stepsize)

platform = Platform.getPlatformByName('CUDA')
#platform = Platform.getPlatformByName('CPU')
simulation=Simulation(pdb.topology,system,integrator,platform)
simulation.context.setPositions(pdb.positions)

## Load in pre-built lambda states
nsubs=[6] # MOVE TO TOP
endstates=[i for i in range(nsubs[0])] # specific to a single site system
Lstates=np.loadtxt(lstate_file)
print("\nThere are",Lstates.shape[0],"lambda states\n")
print("These are the endstates:",endstates,"\n")


## open output file for writing
foutput=open('output','w')
foutput.write("Step energy lambda_index\n")
foutput.flush()


## pick a random lambda to start with
lambda_idx = np.random.randint(Lstates.shape[0])
print("Starting Lambda IDX =",lambda_idx)

# -----------------------------------
##### Minimize and Equilibrate
ibuff=0
for site in range(len(nsubs)):
    for sub in range(nsubs[site]):
        simulation.context.setParameter("lambda"+str(ibuff+1),Lstates[lambda_idx,ibuff]) 
        ibuff+=1


# print beforemin.pdb
state = simulation.context.getState(getPositions = True)
fp=open("beforemin.pdb",'w')
PDBFile.writeFile(psf.topology,state.getPositions(),fp)
fp.close()

state = simulation.context.getState(getEnergy = True)
energy = state.getPotentialEnergy().in_units_of(kilocalorie_per_mole)
print("Initial energy is",energy)
for i in range(5):
    simulation.minimizeEnergy(maxIterations=200)
    state = simulation.context.getState(getEnergy = True)
    energy = state.getPotentialEnergy().in_units_of(kilocalorie_per_mole)
    print("Energy after Min Loop",i,"is",energy)
print("")

# print aftermin.pdb    
state = simulation.context.getState(getPositions = True)
fp=open("aftermin.pdb",'w')
PDBFile.writeFile(psf.topology,state.getPositions(),fp)
fp.close()

## # calc Es at every lambda state
## print("AfterMin Es:",)
## for lam in range(Lstates.shape[0]):
##     ibuff=0
##     for site in range(len(nsubs)):
##         for sub in range(nsubs[site]):
##             simulation.context.setParameter("lambda"+str(ibuff+1),Lstates[lam,ibuff])
##             ibuff+=1
##     state = simulation.context.getState(getEnergy = True, groups = {fgE})
##     energyFG = state.getPotentialEnergy().in_units_of(kilocalorie_per_mole)
##     print(lam, energyFG,)
## print("\n")

foutput=open("output",'w')
foutput.write("Step energy lambda_index\n")

# initial equilibrium
ibuff=0
for site in range(len(nsubs)):
    for sub in range(nsubs[site]):
        simulation.context.setParameter("lambda"+str(ibuff+1),Lstates[lambda_idx,ibuff])
        ibuff+=1
#integrator.step(5000) # 10 ps of sampling
simulation.step(5000) # 10 ps of sampling

state = simulation.context.getState(getEnergy = True, groups = {fgE})
energyFG = state.getPotentialEnergy().in_units_of(kilocalorie_per_mole)
foutput.write("EQ %f %d\n" % (energyFG.value_in_unit(kilocalorie_per_mole),lambda_idx))

# print pdb    
state = simulation.context.getState(getPositions = True)
fp=open("dcd/traj.EQ.pdb",'w')
PDBFile.writeFile(psf.topology,state.getPositions(),fp)
fp.close()

# calc Es at every lambda state
biases=np.zeros(Lstates.shape[0])
energies=np.zeros(Lstates.shape[0])
#print("AfterEQ  Es:",)
for lam in range(Lstates.shape[0]):
    ibuff=0
    for site in range(len(nsubs)):
        for sub in range(nsubs[site]):
            simulation.context.setParameter("lambda"+str(ibuff+1),Lstates[lam,ibuff])
            ibuff+=1
    state = simulation.context.getState(getEnergy = True, groups = {fgE})
    energies[lam] = state.getPotentialEnergy().in_units_of(kilocalorie_per_mole).value_in_unit(kilocalorie_per_mole)

# EQ energies & lambda states are written to disk, but are subsequently overwritten once production begins!
#printState(0,Lstates[lambda_idx],energies) # print info to disk

# sample lambdas at EQ coords
energies=(energies+biases)*beta # to give reduced energies with biases
res=sampleLambda(energies)      # res returns a dictionary of 'prob' (probabilities) and 'idx' (the next lambda index)


# -----------------------------------
##### Start Prod Sampling

# initialize counts
counts=np.zeros(Lstates.shape[0])

# set new lambda & increment bias
#simulation.reporters.append(PDBReporter('output.pdb',numStepPerCycle))
simulation.reporters.append(StateDataReporter(stdout,numStepPerCycle,step=True,
                            potentialEnergy=True,temperature=True,volume=True))

# Wang-Landau Lambda Dynamics
time=0  # step counter
for loop in range(numBiasLoops):
    # open bias, count, lambda, and energy files at the start of the loop
    bfp=open('lambda_energies/bias.'+str(loop)+'.txt','w')
    cfp=open('lambda_energies/count.'+str(loop)+'.txt','w')
    efp=open('lambda_energies/energy.'+str(loop)+'.txt','w')
    lfp=open('lambda_energies/lambda.'+str(loop)+'.txt','w')

    for i in range(numCycle):
        # set up new lambda state
        ibuff=0
        for site in range(len(nsubs)):
            for sub in range(nsubs[site]):
                simulation.context.setParameter("lambda"+str(ibuff+1),Lstates[res['idx'],ibuff])
                ibuff+=1

        # sample for "n" MD steps
        simulation.step(numStepPerCycle)

        # save the energy
        state = simulation.context.getState(getEnergy = True, groups = {fgE})
        energyFG = state.getPotentialEnergy().in_units_of(kilocalorie_per_mole)
        foutput.write("%d %f %d [" % (i,energyFG.value_in_unit(kilocalorie_per_mole),res['idx']))
        ##for lam in range(Lstates.shape[0]):
        ##    foutput.write(" %s" % (str(res['prob'][lam])))
        foutput.write(" ]\n")
        foutput.flush()

        # print coords to pdb for all endstates   
        if res['idx'] in endstates:
            state = simulation.context.getState(getPositions = True)
            #fp=open('dcd/state.'+str(res['idx'])+'.traj.'+str(loop)+'.'+str(i)+'.pdb','w')
            fp=open('dcd/state.'+str(res['idx'])+'.time.'+str(time)+'.pdb','w')
            PDBFile.writeFile(psf.topology,state.getPositions(),fp)
            fp.close()

        # save counts and biases per step
        localcounts=np.zeros(Lstates.shape[0])
        localcounts[res['idx']]+=1
        for cc in range(Lstates.shape[0]):
            cfp.write("%d " % (localcounts[cc])) # stepwise count
        cfp.write("\n")
        for bb in range(Lstates.shape[0]):      
            bfp.write("%f " % (biases[bb]))     # biases
        bfp.write("\n")

        # add bias before sampling lambda
        counts[res['idx']]+=1
        time+=1
        if loop < exp_start:
            biases[res['idx']]+=pre_exp_bias
        else:
            update_exp_biases(biases, bias_fill, counts, exp_base, exp_mult)

        # calc Es at every lambda state
        energies=np.zeros(Lstates.shape[0])
        for lam in range(Lstates.shape[0]):
            ibuff=0
            for site in range(len(nsubs)):
                for sub in range(nsubs[site]):
                    simulation.context.setParameter("lambda"+str(ibuff+1),Lstates[lam,ibuff])
                    ibuff+=1
            state = simulation.context.getState(getEnergy = True, groups = {fgE})
            energies[lam] = state.getPotentialEnergy().in_units_of(kilocalorie_per_mole).value_in_unit(kilocalorie_per_mole)

        # save a copy of the energies if it's the last step
        if i+1 == numCycle:
            final_energies=np.copy(energies)

        # print energies 
        #printState4(loop,Lstates[res['idx']],energies) # print info to disk
        for ee in range(len(energies)):
            efp.write("%f " % (energies[ee]))
        efp.write("\n")
    
        # write out the lambdas
        lambda_state=Lstates[res['idx']]
        for ll in range(len(lambda_state)):
            lfp.write("%.2f " % (lambda_state[ll]))
        lfp.write("\n")

        # flush the output
        cfp.flush()
        bfp.flush()
        efp.flush()
        lfp.flush()

        # sample lambdas at XYZ coords
        energies=(energies+biases)*beta # to give reduced energies with biases
        res=sampleLambda(energies)      # res returns a dictionary of 'prob' (probabilities) and 'idx' (the next lambda index)
     
        print("Finished step",i,"of loop",loop)
    
    # loop finished; save checkpoint
    simulation.saveCheckpoint('checkpoints/save.'+str(loop)+'.chck') 
    # loop finished; close open files
    bfp.close()
    cfp.close()
    efp.close()
    lfp.close()

    # print biases and counts at the end of production
    tfp=open('lambda_energies/total_BCs.'+str(loop)+'.txt','w')
    for BC in range(len(biases)):
        tfp.write("%d %f %d\n" % (BC,biases[BC],counts[BC]))
    tfp.close()

    # calculate FastMBAR free energies, load those in as new biases, and resample lambdas again
    newbiases=calcFastMBAR(Lstates.shape[0],beta,numCycle,(loop+1))
    if loop + 1 < exp_start : 
        biases=np.copy(newbiases)
    else:
        bias_fill = np.copy(newbiases)    # for use in setting floor of exponential biases
        update_exp_biases(biases, bias_fill, counts, exp_base, exp_mult) # this function changes the values in 'biases'
    
    energies=final_energies           # unmodified energies from last step
    energies=(energies+biases)*beta   # get reduced energies
    res=sampleLambda(energies)        # sample lambda with NEW biases

    # no clean up for non-CUDA FastMBAR
    # # clean up memory used by FastMBAR
    cleanMEM()
    

    # end output priting for this loop iteration
    foutput.write("###### BIAS LOOP %d FINISHED ######\n" % (loop))

foutput.close()

exit()

