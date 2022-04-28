#! /usr/bin/env python

import numpy as np

#
# generate the LambdaStates file for a ligand system 
# with substituent changes at 1 site only
#

nsubs=[6]                   # number of substituent changes
output="LambdaStates.txt"   # name of output file

#num_lambdas = 11
#lambs=np.linspace(1,0,num_lambdas)
lambs=[0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00]

# for the endstate strings:
ends=[]
for site in range(len(nsubs)):
    ends.append([])
    print("Endstates for site",site+1)
    for i in range(nsubs[site]):
        ends[site].append([])
        for j in range(nsubs[site]):
            if i == j:
                print(" ",lambs[-1],end="")
                ends[site][-1].append(lambs[-1])
            else:
                print(" ",lambs[0],end="")
                ends[site][-1].append(lambs[0])
        print("")
print("")

# do the intermediate states
intlam=lambs[1:-1]

nsites=[]
for site in range(len(nsubs)):
    nsites.append([])
    nsites[site].append(ends[site])

for ll in range(len(intlam)):
    bb = -1*(ll+1)

    ints=[]
    for site in range(len(nsubs)):
        ints.append([])
        #print("Int State for site",site+1)
        for i in range(nsubs[site]):
            for j in range(i+1,nsubs[site]):
                ints[site].append([])
                for c in range(nsubs[site]):
                    if c == i:
                        #print(" ",intlam[bb],end="")
                        ints[site][-1].append(intlam[bb])
                    elif c == j:
                        #print(" ",intlam[ll],end="")
                        ints[site][-1].append(intlam[ll])
                    else: # i > j:
                        #print(" ",0.00,end="")
                        ints[site][-1].append(0.00)
                #print("")
        nsites[site].append(ints[site])
    #print("")

#print(len(nsites))
#print([len(x) for x in nsites])
#print([len(y) for x in nsites for y in x])

# Now print it all out together
fp=open('LambdaStates.txt','w')

# get a combination of states: all of nsubs[0] with all of nsubs[1]
site1_states=[y for x in nsites[0] for y in x]
#print(site1_states)
print("A total of",len(site1_states),"discrete lambda states were generated")
for line in site1_states:
    for value in line:
        fp.write(" %.2f" % (value))
    fp.write("\n")

fp.close()



