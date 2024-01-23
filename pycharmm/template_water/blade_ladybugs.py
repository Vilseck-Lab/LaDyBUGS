##
## pyCHARMM script 
##
## ** LaDyBUGS with pyCHARMM+BLADE **
##

import os
import sys
import glob
import numpy as np
import pandas
from functions import cleanMEM, sampleLambda, update_exp_biases, calcFastMBAR

##############################################
# Load pyCHARMM libraries

import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.scalar as scalar
from pycharmm.lib import charmm as libcharmm


##############################################
# Set up global parameters

# variables
sysname  = 'water'
box      = ??.??            # angstrom - reqs user input
pmegrid  = ??.??           # angstrom - reqs user input
temp     = 298.15           # kelvin
builddir = './prep'

# nonbonded conditions
nb_fswitch = False          # normal fswitching functions
nb_pme   = True             # normal PME

# dynamics conditions
blade    = True
no_save  = True             # True = do NOT save dcd or res files in prod sampling

# dynamics variables
cpt_on   = True             # run with CPT for NPT?
timestep = 0.002            # ps
nequil   = 50000            # equil for 100 ps
nsavc    = 10000            # dcd save frequency

# LaDyBUGS sampling parameters  (This example will run for 25 ns)
# total length of sampling = nprod * nGSloop * nBiasloop * stepsize (in ps)
nprod     = 100             # number of Molecular Dynamics steps per Gibbs Sampler step     (~numStepPerCycle)
nGSloop   = 1000             # number of Gibbs Sampler loops between bias updates w/FastMBAR (~numCycle)
nBiasloop = 125             # number of bias updates with FastMBAR                          (~numBiasLoops)

# LaDyBUGS exponential bias function variables;  exp_bias = exp_mult * exp_base ** (counts[i] - lowest_count)
pre_exp_bias = 10.0
exp_base     = 2.0          # number to be exponentiated
exp_mult     = 1.0
exp_start    = 1            # the first nBiasloop update (0-ordered) where exp_biasing is used  
kb = 1.987204259 * 0.001    # kcal/(mol*K)
beta = 1/(temp * kb)


# ld variables
nblocks=np.loadtxt('./nblocks',dtype='int')
nsubs=np.loadtxt('./nsubs',dtype='int',ndmin=1)
nsites=len(nsubs)

# Load in (pre-built) lambda states; assumes endstates == first nblock lambda states
lstate_file = './LambdaStates.txt'
endstates=[i for i in range(nsubs[0])] # specific to a single site system
Lstates=np.loadtxt(lstate_file)
print("\nLDB: There are",Lstates.shape[0],"lambda states\n")
print("LDB: These are the endstates:",endstates,"\n")


## pick a random lambda to start with
#lambda_idx = np.random.randint(Lstates.shape[0])
lambda_idx = 0
print("LDB: Starting Lambda IDX =",lambda_idx)


# ldb LIGAND perturbations
lpert={}                    # dict of nested dicts; index = site # (zero indexed)
lpert[0]={'subs':['1','2','3','4','5','6'],  # site1_sub# (0-index = site 1)
         'segid':'LIG',
         'resid':'1'}


##############################################
# Read in toppar files, coordinate files, etc.

# toppar files
#settings.set_bomb_level(-2)
pycharmm.lingo.charmm_script('bomblev -2')
pycharmm.lingo.charmm_script('stream '+builddir+'/toppar.str')
read.rtf(builddir+'/core.rtf', append=True)
read.prm(builddir+'/full_ligand.prm', append=True, flex=True)

#settings.set_bomb_level(0)
pycharmm.lingo.charmm_script('set builddir = {}'.format(builddir))


# read in ligand core
read.sequence_pdb(builddir+'/core.pdb')
gen.new_segment(lpert[0]['segid'], setup_ic=True)   # (SEGID, FIRST, LAST, [options])
read.pdb(builddir+'/core.pdb',resid=True)


# load in patches for alchem residues (lig)
pycharmm.lingo.charmm_script('ic generate')
for site in range(nsites):
    for sub in range(nsubs[site]): 
        read.rtf(builddir+'/site'+str(site+1)+'_sub'+lpert[site]['subs'][sub]+'_pres.rtf', append=True)
        pycharmm.lingo.charmm_script('''
patch p{}_{} {} {} setup
read coor pdb resid name {}
ic param
ic build
'''.format(str(site+1),
           lpert[site]['subs'][sub],
           lpert[site]['segid'],
           lpert[site]['resid'],
           builddir+'/site'+str(site+1)+'_sub'+lpert[site]['subs'][sub]+'_frag.pdb')
)

# read in LonePair sites (if applicable)
pycharmm.lingo.charmm_script('stream '+builddir+'/lpsites.inp')

# define ligand substituent selections
#pycharmm.lingo.charmm_script('stream '+builddir+'/ligand_selections.str')
for site in range(nsites):
    lpert[site]['select']=[]
    for sub in range(nsubs[site]): 
        # append charmm variable name for substituent selection
        lpert[site]['select'].append('site'+str(site+1)+'_sub'+lpert[site]['subs'][sub])
        # extract alchem patch atoms from patch file
        sub_atoms=''
        rtffile = builddir+'/site'+str(site+1)+'_sub'+lpert[site]['subs'][sub]+'_pres.rtf'
        for line in open(rtffile,'r'):
            if line[0:4] == 'ATOM': sub_atoms=sub_atoms+line.split()[1].upper()+' '
        sub_atoms=sub_atoms[:-1] # remove trailing space
        atoms_in_sub = pycharmm.SelectAtoms().by_res_and_type(lpert[site]['segid'],lpert[site]['resid'],sub_atoms)
        select.store_selection(lpert[site]['select'][sub],atoms_in_sub) # saves the identified atoms in the charmm variable


# delete angles and dihedrals between alchem groups
pycharmm.lingo.charmm_script('auto angle dihe')
#settings.set_bomb_level(-1)
pycharmm.lingo.charmm_script('bomblev -1')
for site in range(nsites):
    for sub1 in range(nsubs[site]):
        for sub2 in range(sub1+1,nsubs[site]):
            pycharmm.lingo.charmm_script('dele connectivity sele {} show end sele {} show end'
            .format(lpert[site]['select'][sub1],lpert[site]['select'][sub2]))


# read in solvent and ions
read.sequence_pdb(builddir+'/solvent.36.583.pdb')
gen.new_segment('WT00', setup_ic=True, angle=False, dihedral=False)
read.pdb(builddir+'/solvent.36.583.pdb',resid=True)

pycharmm.lingo.charmm_script('print coor sele .not. init end')
#settings.set_bomb_level(0)
pycharmm.lingo.charmm_script('bomblev 0')

# write out psf, crd, pdb files
write.psf_card('patch.psf')
write.coor_card('patch.crd')
write.coor_pdb('patch.pdb')


##############################################
# Create water box & periodic images

crystal.define_cubic(box)
crystal.build(14.0)

# center at 0.0 for blade and regular charmm; center at boxhalf for omm!
image.setup_segment(0.0, 0.0, 0.0, 'LIG')  # list out one by one
image.setup_residue(0.0, 0.0, 0.0, 'WT00')

pycharmm.lingo.charmm_script('coor copy comp') # for truncated system


##############################################
# Set up BLOCK module for MSLD

# check that the system's net charge is 0
pycharmm.lingo.charmm_script('set charge = ?cgtot')
netQ = pycharmm.lingo.get_charmm_variable('CHARGE')


# STARTING BLOCK module
# ** the block module must be passed by pycharmm.lingo as one complete unit. 
#    It CANNOT be divided into parts
# ** Therefore, multiple strings are created and passed at once to lingo
blockplusone = nblocks + 1
# initialize block
block_init='''
!! BLOCK setup
BLOCK {}
   clear
END
BLOCK {}
'''.format(blockplusone,blockplusone)
# load blocks
block_call=''
lpert[site]['block']=[]  # could also use a dict
ii=2
for site in range(nsites):
    for sub in range(nsubs[site]):
        block_call+='Call {} sele {} show end\n'.format(ii,lpert[site]['select'][sub])
        lpert[site]['block'].append(ii)
        ii+=1

knoe = 118.4  
block_parm='''
! scat on
! scat k {}
! cats sele atom ?segid ?resid ?atomname .or. [list of atom names to cat]

qldm theta
lang temp {}
soft w14
pmel ex ! no pme for LDB vs OMM comparison

ldin 1 1.0  0.0  5.0  0.0  5.0'''.format(knoe,temp)
# ldin lines
block_ldin=''
sitestr=''
ibuff=0
for site in range(nsites):
    for sub in range(nsubs[site]):
        block_ldin+='ldin {} {:.4f} 0.0 5.0 0 5.0\n'.format(lpert[site]['block'][sub],Lstates[lambda_idx,ibuff])
        sitestr+=str(site+1)+' '
        ibuff+=1
# add in exclusions with adex
block_adex=''
for site in range(nsites):
    for sub1 in range(nsubs[site]):
        for sub2 in range(sub1+1,nsubs[site]):
            block_adex+='adex {} {}\n'.format(lpert[site]['block'][sub1],lpert[site]['block'][sub2])

# msld parameters
block_msld='''
!! remove printing of lambdas in restart file
RLFR OFF 

!!rmla bond thet dihe impr
rmla bond thet impr
msld 0  {} ffix
msma
'''.format(sitestr)

# msld variable biases
block_varb='''
ldbi  0       ! no biasing potential
END
'''

pycharmm.lingo.charmm_script('''
{}
{}
{}
{}
{}
{}
{}'''.format(block_init,block_call,block_parm,block_ldin,block_adex,block_msld,block_varb))

def block_change_Lstate(lambda_state_idx):
    ''' Update the lambda state within BLOCK -
        No other changes made; all exclusions, settings, etc. 
        are maintained from BLOCK initialization
    '''
    # ldin lines
    new_block_ldin='ldin 1 1.0  0.0  5.0  0.0  5.0\n'
    ibuff=0
    for site in range(nsites):
        for sub in range(nsubs[site]):
            new_block_ldin+='ldin {} {:.4f} 0.0 5.0 0 5.0\n'.format(lpert[site]['block'][sub],Lstates[lambda_state_idx,ibuff])
            ibuff+=1
    # pass changes to charmm
    pycharmm.lingo.charmm_script('''
    BLOCK
    {}
    msma
    END'''.format(new_block_ldin))

    return


##############################################
# Set NonBonded settings & SP energy calc

cutnb = 14.0
cutim = cutnb
ctofnb = 12.0
ctonnb = 10.0

## nbond switching
## use a dictionary so that it becomes easy to switch between w/ vs w/o PME
nbonds_dict = {'cutnb':cutnb,'cutim':cutim,
           'ctonnb':ctonnb,'ctofnb':ctofnb,
           'atom':True,'vatom':True,
           'cdie':True,'eps':1.0,
           'inbfrq':-1, 'imgfrq':-1}

if nb_pme:
    nbonds_dict['switch']   = True
    nbonds_dict['vfswitch'] = True
    nbonds_dict['ewald']    = True
    nbonds_dict['pmewald']  = True
    nbonds_dict['kappa']    = 0.32
    nbonds_dict['fftx']     = pmegrid
    nbonds_dict['ffty']     = pmegrid
    nbonds_dict['fftz']     = pmegrid
    nbonds_dict['order']    = 4

elif nb_fswitch:
    nbonds_dict['fswitch']  = True
    nbonds_dict['vfswitch'] = True
    nbonds_dict['ewald']    = False
    nbonds_dict['pmewald']  = False

else: 
    print("NonBonded Parameter Error - both pme and switch options are false")
    pycharmm.lingo.charmm_script('stop')

nbonds=pycharmm.NonBondedScript(**nbonds_dict)
nbonds.run()
energy.show()

##############################################
# Minimize the system

if len(glob.glob('minimized.crd')) == 1:    # read in minimized coordinates
    pycharmm.lingo.charmm_script('read coor card resid name minimized.crd')
else:                                       # minimized & write out psf, crd, pdb files
    minimize.run_sd(nstep=250,nprint=50,step=0.005,tolenr=1e-3,tolgrd=1e-3)
    write.psf_card('minimized.psf')
    write.coor_card('minimized.crd')
    write.coor_pdb('minimized.pdb')

energy.show()

#minimize.run_abnr(nstep=250,nprint=50,tolenr=1e-3,tolgrd=1e-3)
#energy.show()



##############################################
# Set up and run Dynamics

# dynamics conditions
if blade:
    useblade = 'prmc pref 1 iprs 100 prdv 100'  # only for CPT - don't use any of these for not CPT (just blade keyword)
    gscal = 0.1
    ntrfrq=0
    leap = True
    openmm = False
else: 
    print("MSLD can only be run with BLADE - exiting...")
    pycharmm.lingo.charmm_script('stop')

# set shake
shake.on(bonh=True,fast=True,tol=1e-7)
dyn.set_fbetas(np.full((psf.get_natom()),gscal,dtype=float))

# initialize blade
pycharmm.lingo.charmm_script('energy blade')

# set up output directories
if not os.path.isdir('res'): os.system('mkdir res')
if not os.path.isdir('dcd'): os.system('mkdir dcd')
if not os.path.isdir('lambda_energies'): os.system('mkdir lambda_energies')
if not os.path.isdir('mbar_results'): os.system('mkdir mbar_results')

# set up dynamics dictionary of parameters
dynamics_dict = {'cpt':cpt_on,'leap':True,'langevin':True,
    'nsavl':nprod,  # frequency for saving lambda values in lamda-dynamics
    'timestep':timestep,'ntrfrq':ntrfrq,
    'firstt':temp,'finalt':temp,'tstruct':temp,'tbath':temp,
    'iasors': 1,'iasvel':1,'iscvel': 0,'iscale': 0,
    'ihtfrq':0,'ieqfrq':0,'ichecw': 0,
    'inbfrq':0,'imgfrq':0,'ihbfrq':0,'ilbfrq':0,
    'echeck': -1}

#   !! INBFRQ and IMGFRQ should be 0 (ZERO) to avoid calling UPDATE everytime DYNA restarts below !!!!

if blade:
    dynamics_dict['omm']   = False
    dynamics_dict['blade'] = useblade


# MD equilibration
dcd_file = pycharmm.CharmmFile(file_name='dcd/{}.{}.dcd'.format(sysname,0), 
               file_unit=1,formatted=False,read_only=False)
res_file = pycharmm.CharmmFile(file_name='res/{}.{}.res'.format(sysname,0),
               file_unit=2,formatted=True,read_only=False)
lam_file = pycharmm.CharmmFile(file_name='res/{}.{}.lam'.format(sysname,0), 
               file_unit=3,formatted=False,read_only=False)

dynamics_dict['start']  = True
dynamics_dict['nstep']  = nequil
dynamics_dict['iunrea'] = -1
dynamics_dict['iunwri'] = res_file.file_unit
dynamics_dict['iuncrd'] = dcd_file.file_unit
dynamics_dict['iunldm'] = lam_file.file_unit
dynamics_dict['nsavc']  = nsavc   # use nsavc for equilibration; nprod for production (if saving dcd)
dynamics_dict['isvfrq'] = nsavc   # use nsavc for equilibration; nprod for production (if saving res)
dynamics_dict['nprint'] = nsavc   # Frequency to write to output
dynamics_dict['iprfrq'] = nsavc   # Frequency to calculate averages

equil_dyn = pycharmm.DynamicsScript(**dynamics_dict)
equil_dyn.run()

dcd_file.close()
res_file.close()
lam_file.close()

write.coor_pdb('dcd/{}.{}.pdb'.format(sysname,0)) # write out final frame

###############################
### LaDyBUGS sampling starts:
# - calculate the energy of every lambda state after equil
# - pick a new lambda value to start production sampling with

# calc Es at every lambda state
biases=np.zeros(Lstates.shape[0])
energies=np.zeros(Lstates.shape[0])
for lam in range(Lstates.shape[0]):
    block_change_Lstate(lam)
    pycharmm.lingo.charmm_script('energy blade abic')
    energies[lam] = pycharmm.lingo.get_energy_value('ENER')


# sample lambdas at EQ coords
energies=(energies+biases)*beta # to give reduced energies with biases
res=sampleLambda(energies)      # res returns a dictionary of 'prob' (probabilities) and 'idx' (the next lambda index)
print("LDB: (After equilibration) The new lambda state is",res['idx'])


# LaDyBUGS Production
dynamics_dict['start']   = False
dynamics_dict['restart'] = True
dynamics_dict['nstep']   = nprod

dynamics_dict['nsavc']   = nprod  # frequency for saving dcd coordinates
dynamics_dict['isvfrq']  = nprod  # Frequency to save restart file
dynamics_dict['nprint']  = nprod  # Frequency to write to output
dynamics_dict['iprfrq']  = nprod  # Frequency to calculate averages

dynamics_dict['iasors']  = 0      # 0 scale vel; >0 assign vel
dynamics_dict['iasvel']  = 0      # 0 == assign vel using comparison coordinate values
dynamics_dict['iscvel']  = 0      # 0 single scale factor
dynamics_dict['iscale']  = 0      # 0 no scaling; >0 scale by that #

dynamics_dict['abic']    = True   # assume blade is current

prod_dyn = pycharmm.DynamicsScript(**dynamics_dict)   # load new dynamics settings

# io check
if no_save:
    if dynamics_dict['restart']:
        dynamics_dict['restart'] = False
        dynamics_dict['start']   = True


# Lambda Dynamics with Bias Updated Gibbs Sampling
time=0  # step counter
counts=np.zeros(Lstates.shape[0]) # initialize counts
for loop in range(nBiasloop):
    # open bias, count, lambda, and energy files at the start of the loop
    bfp=open('lambda_energies/bias.'+str(loop)+'.txt','w')
    cfp=open('lambda_energies/count.'+str(loop)+'.txt','w')
    efp=open('lambda_energies/energy.'+str(loop)+'.txt','w')
    lfp=open('lambda_energies/lambda.'+str(loop)+'.txt','w')

    for i in range(nGSloop):
        # set up new lambda state
        block_change_Lstate(res['idx'])

        # sample for "nprod" MD steps
        lam_file = pycharmm.CharmmFile(file_name='res/{}.{}.lam'.format(sysname,time+1), 
                   file_unit=3,formatted=False,read_only=False)
        if no_save: # do not save dcd or res files
            dynamics_dict['iunrea'] = -1
            dynamics_dict['iunwri'] = -1
            dynamics_dict['iuncrd'] = -1
            dynamics_dict['iunldm'] = lam_file.file_unit
        else:
            dcd_file = pycharmm.CharmmFile(file_name='dcd/{}.{}.dcd'.format(sysname,time+1), 
                       file_unit=1,formatted=False,read_only=False)
            res_file = pycharmm.CharmmFile(file_name='res/{}.{}.res'.format(sysname,time+1), 
                       file_unit=2,formatted=True,read_only=False)
            prv_rest = pycharmm.CharmmFile(file_name='res/{}.{}.res'.format(sysname,time), 
                       file_unit=4,formatted=True,read_only=False)
            dynamics_dict['iunrea'] = prv_rest.file_unit
            dynamics_dict['iunwri'] = res_file.file_unit
            dynamics_dict['iuncrd'] = dcd_file.file_unit
            dynamics_dict['iunldm'] = lam_file.file_unit

        prod_dyn = pycharmm.DynamicsScript(**dynamics_dict)   # load new dynamics settings
        prod_dyn.run()

        lam_file.close()
        if not no_save:
            dcd_file.close()
            res_file.close()
            prv_rest.close()

        # print coords to pdb for all endstates   
        if res['idx'] in endstates:
            write.coor_pdb('dcd/state.{}.time.{}.pdb'.format(str(res['idx']),str(time))) # write out final frame

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
        old_level = settings.set_verbosity(1) # sets prnlev 1
        energies=np.zeros(Lstates.shape[0])
        for lam in range(Lstates.shape[0]):
            block_change_Lstate(lam)
            pycharmm.lingo.charmm_script('energy blade abic')
            energies[lam] = pycharmm.lingo.get_energy_value('ENER')
        settings.set_verbosity(old_level) # back to prnlev [default]

        # save a copy of the energies if it's the last step
        if i+1 == nGSloop:
            final_energies=np.copy(energies)

        # print energies 
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
     
        print("LDB: Finished step",i,"of loop",loop)
    

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
    newbiases=calcFastMBAR(Lstates.shape[0],beta,nGSloop,(loop+1))
    if loop + 1 < exp_start : 
        biases=np.copy(newbiases)
    else:
        bias_fill = np.copy(newbiases)    # for use in setting floor of exponential biases
        update_exp_biases(biases, bias_fill, counts, exp_base, exp_mult) # this function changes the values in 'biases'
    
    energies=final_energies           # unmodified energies from last step
    energies=(energies+biases)*beta   # get reduced energies
    res=sampleLambda(energies)        # sample lambda with NEW biases

    # clean up memory used by FastMBAR (not needed for non-CUDA FastMBAR)
    cleanMEM()
    
    # end output printing for this loop iteration
    print("LDB: ###### BIAS LOOP %d FINISHED ######\n" % (loop))


# clean up restraints
#cons_harm.turn_off()
#pycharmm.lingo.charmm_script('cons harm clear')


##############################################
# FINISHED

pycharmm.lingo.charmm_script('stop')


