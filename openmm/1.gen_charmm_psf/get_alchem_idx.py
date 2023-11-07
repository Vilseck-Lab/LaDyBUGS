#! /usr/bin/env python

###
### Read the charmm output file and the psf file to get the
### atom indices for all alchemical groups
###

outfile = 'prep/msld_sele.str'
patfile = 'patch.psf'

## read the output file to pull out the alchemical groups
alchem_grps=[]  # group numbers ranging from 1...n
alchem_info={}  # info about each alchem group
alchem_atoms=[] # a list of ALL alchem atoms
grp=0
fp=open(outfile,'r')
line=fp.readline()
while(line): # following lines specific to msld_sele.str format
    if line[0:11] == 'define site':
        grp+=1
        alchem_grps.append(grp)
        alchem_info[grp]={'segid':'LIG','resid':'1','resn':'LIG','atoms':[]}
        while(line):
            if line[0:15] == '   atom @ligseg':
                tmp=line.split()
                alchem_atoms.append(tmp[3])
                alchem_info[grp]['atoms'].append(tmp[3])
            if line[0:13] == '   none ) end': break
            line=fp.readline()
    line=fp.readline()
fp.close()

# # debug
# print(alchem_grps)
# print(alchem_info)
# print(len(alchem_atoms))
# print(alchem_atoms)
# quit()

## read the patch file and look for atoms
at2idx={} # atom to index dictionary
fp=open(patfile,'r')
for a in range(8): # skip reading the first few lines
    line=fp.readline()
chk=0
while(line):
    tmp=line.split()
    if tmp[4] in alchem_atoms:
        # loop through atom lists for the groups; and double check that segid/resid/resn all match
        for g in alchem_grps:
            pchk=chk
            if (tmp[1] == alchem_info[g]['segid']) and (tmp[2] == alchem_info[g]['resid']) and (tmp[3] == alchem_info[g]['resn']):
                for at in alchem_info[g]['atoms']:
                    if tmp[4] == at:
                        chk+=1
                        at2idx[at]=int(tmp[0])
                        break
            if pchk < chk: # atom successfully identified - move on
                break
    if chk == len(alchem_atoms):
        break  # all alchem atoms found
    else:
        line=fp.readline()

#print(chk,len(alchem_atoms))
#print(at2idx)

## print out the alchemical atoms with their indices
print("Atom Name List")
for g in alchem_grps:
    print('name.append([',end="")
    for at in range(len(alchem_info[g]['atoms'])-1):
        print('\''+alchem_info[g]['atoms'][at]+'\',',end="")
    print('\''+alchem_info[g]['atoms'][-1]+'\'',end="")
    print('])')
print("")

print("Atom Index List (1-index matches psf)")
for g in alchem_grps:
    print('nameid.append([',end="")
    for at in range(len(alchem_info[g]['atoms'])-1):
        #print('\''+str(at2idx[alchem_info[g]['atoms'][at]])+'\',',end="")
        print(str(at2idx[alchem_info[g]['atoms'][at]])+',',end="")
    #print('\''+str(at2idx[alchem_info[g]['atoms'][-1]])+'\'',end="")
    print(str(at2idx[alchem_info[g]['atoms'][-1]]),end="")
    print('])')
print("")

print("Atom Index List (0-index matches omm)")
for g in alchem_grps:
    print('nameid.append([',end="")
    for at in range(len(alchem_info[g]['atoms'])-1):
        #print('\''+str(at2idx[alchem_info[g]['atoms'][at]])+'\',',end="")
        print(str(at2idx[alchem_info[g]['atoms'][at]]-1)+',',end="")
    #print('\''+str(at2idx[alchem_info[g]['atoms'][-1]])+'\'',end="")
    print(str(at2idx[alchem_info[g]['atoms'][-1]]-1),end="")
    print('])')
print("")

## print bonded lists for visualizing OMM-WLLD pdb trajectories in PYMOL
op=open('../fix_pymol_bonds.txt','w')
chainID='B' # hardcoded
for g1 in range(len(alchem_grps)):
    for g2 in range(g1+1,len(alchem_grps)):
        for at1 in alchem_info[alchem_grps[g1]]['atoms']:
            for at2 in alchem_info[alchem_grps[g2]]['atoms']:
                op.write("unbond chain %s and resid %s and name %s, chain %s and resid %s and name %s\n" % (chainID,alchem_info[alchem_grps[g1]]['resid'],at1,chainID,alchem_info[alchem_grps[g2]]['resid'],at2))
op.close()

## print "hide state" scripts for visualizing things in PYMOL
for g in range(len(alchem_grps)):
    op=open('../hide_state.'+str(g)+'.txt','w')
    for at in alchem_info[alchem_grps[g]]['atoms']:
        op.write("hide everything, chain %s and resid %s and name %s\n" % (chainID,alchem_info[alchem_grps[g]]['resid'],at))
    op.close()

# finished

