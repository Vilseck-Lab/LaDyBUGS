import numpy as np
import torch
from FastMBAR import *
import glob,sys

# THIS SCRIPT ASSUMES THAT ALL STATES ARE SAMPLED

# find production run directories - given as cmd ln argu
dir_list = glob.glob(sys.argv[1]+'*')
#print(dir_list)

# import lambda states
num_lambstates=0
for line in open(dir_list[0]+'/LambdaStates.txt'): num_lambstates+=1

# find nsubs
num_endstates = nsubs=int(np.loadtxt(dir_list[0]+'/nsubs',dtype=int))

# output file name
output_file_name = sys.argv[1]+'_final_results.txt'

temp = 298.15 #  in K, make sure this matches your production runs
kb = 0.0019872041        # kcal/(mol*K)


# refers to each Gibbs sampling loop in LaDyBUGS
prod_start = 0
# find total number of mbar loops in LDB
prod_end = len(glob.glob(dir_list[0]+'/mbar_results/mbar.*.txt'))

# verbose output
print("automatically determined variables:")
print("dir_list",dir_list)
print("num_lambstates",num_lambstates)
print("num_endstates",num_endstates)
print("prod_end",prod_end)

# get energies and counts from all runs in directory_list
def gather_MBAR_inputs(directory_list, num_lambstates, prod_start, prod_end):
    global_counts = []
    global_ener = []
    for directory in directory_list:
        count_file_list = glob.glob(directory+'/lambda_energies/energy.*.txt')
        numBiasLoops = len(count_file_list)
        if prod_end == -1:
            prod_end = numBiasLoops
        N_k = np.zeros(num_lambstates, dtype=int)
        for loop in range(prod_start, prod_end):
            filename = directory+'/lambda_energies/count.'+str(loop)+'.txt'
            N_tmp = np.loadtxt(filename, dtype=int)
            if loop == prod_start:
                numCycle = len(N_tmp) # assumes the num of GS steps never changes per bias loop
                ener = np.zeros((numCycle*(int(prod_end)-int(prod_start)), num_lambstates), dtype=np.float64)
                print(ener.shape)
            N_tmp = N_tmp.sum(0)
            N_k = N_k + N_tmp
            ener_file = directory+'/lambda_energies/energy.'+str(loop)+'.txt'
            ener[(loop-prod_start)*numCycle:(loop+1-prod_start)*numCycle] = np.loadtxt(ener_file)
        
        if len(global_counts) == 0:
            global_counts = np.copy(N_k)
        else:
            global_counts = global_counts + N_k

        if len(global_ener) == 0:
            global_ener = np.copy(ener)
        else:
            global_ener = np.concatenate((global_ener, ener), axis = 0)
    return(global_counts, global_ener)

def calcFastMBAR(beta,N_k, ener):
    do_cuda=True
    do_cuda_batch=True 
    block_size=1

    U_k = ener.T * beta
    
    # == run FastMBAR: (i) get exact fe value
    do_bootstrap=False
    fastmbar = FastMBAR(energy = U_k, num_conf = N_k, cuda=do_cuda, cuda_batch_mode=do_cuda_batch, bootstrap=do_bootstrap, bootstrap_block_size = block_size, bootstrap_num_rep = 150) 
    F = fastmbar.F
    
    # == run FastMBAR: (ii) get bootstrapped uncertainties
    do_bootstrap=True
    fastmbar = FastMBAR(energy = U_k, num_conf = N_k, cuda=do_cuda, cuda_batch_mode=do_cuda_batch, bootstrap=do_bootstrap, bootstrap_block_size = block_size, bootstrap_num_rep = 150) 
    if do_bootstrap:
        F_std=fastmbar.F_std
    
    # remove "reduction"
    rF = F / beta
    if do_bootstrap:
        rSD = F_std / beta

    # normalize to 1st lambda state
    rF = rF-rF[0] 

    if do_bootstrap == True:    
        return(rF, rSD)
    else:
        return(rF)

beta = (temp * kb) ** -1.0

# main operations
N_k, ener = gather_MBAR_inputs(dir_list, num_lambstates, prod_start, prod_end)
rF, rSD = calcFastMBAR(beta, N_k, ener)

# print output
fp=open(output_file_name, 'w')
for l in range(num_endstates):
        fp.write("%d %.6f %.3f\n" % (l+1,rF[l],rSD[l]))
        #fp.write("%d %.6f\n" % (l+1,rF[l])) # with no bootstrapping
fp.close()


print('FINISHED')

