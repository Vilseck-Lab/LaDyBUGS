__author__ = "JV, MTR"
__date__ = "2020/04/09, 2021/06/02"

import numpy as np
import torch
from FastMBAR import FastMBAR
#from guppy import hpy
#from memory_profiler import profile
#import sys

# ----------------------------------------------------------------------------
def cleanMEM():
    """
    Clean up any memory allocated by torch on the GPU
    """
    torch.cuda.empty_cache()
    return
# ----------------------------------------------------------------------------
def sampleLambda(reducedEnergy):
    reducedEnergy = np.array(reducedEnergy)
    minReducedEnergy = np.min(reducedEnergy)
    reducedEnergy = reducedEnergy - minReducedEnergy
    prob = np.exp(-reducedEnergy)
    prob = prob / np.sum(prob)
    sampledLambdaIdx = np.random.choice(len(reducedEnergy),1,p=prob)
    res = {}
    res['idx'] = sampledLambdaIdx[0]
    res['prob'] = prob
    return res
# ----------------------------------------------------------------------------
def update_exp_biases(biases, bias_fill, counts, exp_base, exp_mult):
    lowest_count = np.amin(counts)
    for i in range(len(biases)):
        exp_bias = exp_mult * exp_base ** (counts[i] - lowest_count)
        biases[i] = bias_fill[i] + exp_bias
# ----------------------------------------------------------------------------
def calcFastMBAR(lambstates,beta,numCycle,numBiasLoops):
    print("Running FastMBAR after",numCycle*numBiasLoops,"steps")

    # load in the counts, biases, and energies
    ener=np.zeros((numCycle*numBiasLoops,lambstates), dtype=np.float64)
    N_k=np.zeros(lambstates, dtype=int)
 
    ibuff = 0
    # record the index that was sampled at each GS step
    list_sampled = [] # gives the index of lambda state that was sampled at each GS step
    for loop in range(numBiasLoops):
        # Assemble the counts
        filename='./lambda_energies/count.'+str(loop)+'.txt'
        N_tmp = np.loadtxt(filename,dtype=int)
        N_tmp = N_tmp.sum(0)
        N_k = N_k + N_tmp
        # Assemble the energies
        filename='./lambda_energies/energy.'+str(loop)+'.txt'
        ener[loop*numCycle:(loop+1)*numCycle]=np.loadtxt(filename)

#   change energies into "reduced" energies
#   ener must also be transposed
    #kb = 1.987204259 * 10^-3
    #beta = (kb*temp)^-1
    u_kn = ener.T * beta

    # Unsort sampled vs unsampled states for new FastMBAR analysis 
    config=numCycle*numBiasLoops
    Usampled=np.zeros((lambstates,config))
    Uunsampled=np.zeros((lambstates,config))
    ibuff=0 # rolling index for Usampled
    jbuff=0 # rolling index for Uunsampled
    SamIDX={} # old index -> new index
    USamIDX={} # old index -> new index

    Nsampled=[]
    for n in range(N_k.shape[0]):
        if N_k[n] == 0: # then it is an unsampled state
            USamIDX[n]=jbuff
            Uunsampled[jbuff,:]=u_kn[n,:]
            jbuff+=1
        else: #N_k[n] == 1 # it is a sampled state
            SamIDX[n]=ibuff
            Usampled[ibuff,:]=u_kn[n,:]
            ibuff+=1
            Nsampled.append(N_k[n]) # keeps counts
    Nsampled=np.asarray(Nsampled)

    # trim off empty values
    num_unsampled=len(USamIDX.keys())
    if num_unsampled > 0:
        Usampled=Usampled[0:Usampled.shape[0]-num_unsampled,:]
        Uunsampled=Uunsampled[0:num_unsampled]

    do_cuda=True
    do_cuda_batch=False
    do_bootstrap=False
    # construct a FastMBAR object with the energy matrix and the number of configuration array
    fastmbar = FastMBAR(energy = Usampled, num_conf = Nsampled, cuda=do_cuda, cuda_batch_mode=do_cuda_batch, bootstrap=do_bootstrap) 

    # calculate free energies by solving the MBAR equations
    F = fastmbar.F
    if do_bootstrap:
        Fstdev=fastmbar.F_std

    # Calculate Unsampled/Perturbed states
    if num_unsampled > 0:
        UF,UFstdev = fastmbar.calculate_free_energies_of_perturbed_states(Uunsampled)

    # Now unmix the separated sampled and unsampled states back into their original order
    # if every state is sampled at least once, Fcomplete will the the same as F
    SamIDXkeys=list(SamIDX.keys())
    Fcomplete=np.zeros((lambstates))
    #Fcomplete=np.zeros((lambstates*config))
    for n in range(N_k.shape[0]):
        if n in SamIDXkeys:
            Fcomplete[n]=F[SamIDX[n]]
        else:
            Fcomplete[n]=UF[USamIDX[n]]
    if do_bootstrap:
        #SDcomplete=np.zeros((lambstates*config))
        SDcomplete=np.zeros((lambstates))
        for n in range(N_k.shape[0]):
            if n in SamIDXkeys:
                SDcomplete[n]=Fstdev[SamIDX[n]]
            else:
                SDcomplete[n]=UFstdev[USamIDX[n]]

    # remove "reduction"
    rF = Fcomplete*1/beta
    if do_bootstrap:
        rSD = SDcomplete*1/beta

    # make first value the relative 0.0 (if it isn't already)
    rF=rF-rF[0]

    # remove "reduction"
    rF = Fcomplete*1/beta
    if do_bootstrap:
        rSD = StDev*1/beta
 
    # make first value the relative 0.0 (if it isn't already)
#    if N_k[0] == 0:
    rF=rF-rF[0]
 
    # # debug print
    # print("CHECK rF",rF.shape)
    # print(rF)
 
    # print all FEs to a file and new_biases for additional runs
    fp=open("mbar_results/mbar."+str(loop)+".txt",'w')
    for l in range(lambstates):
        if do_bootstrap:
            fp.write("%d %.6f %.3f\n" % (l+1,rF[l],rSD[l]))
        else:
            fp.write("%d %.6f\n" % (l+1,rF[l]))
    fp.close()
 
    # calculate new biases
    newbiases=rF*-1

    #gp=open("new_biases.txt",'w')
    #for l in range(lambstates):
    #    gp.write("%.6f\n" % (rF[l]*-1))
    #gp.close()

    # finished
    return newbiases

