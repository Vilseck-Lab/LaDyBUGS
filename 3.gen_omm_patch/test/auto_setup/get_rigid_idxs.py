#! /usr/bin/env python2

# Try to translate the rigidbb file into atom indices

# open files
rigid="rigidbb.txt"
pdb="patch.pdb"

# read rigid; assume no blank lines in rigid and that everything is in order
ratoms={}
keylist=[]
fp=open(rigid,'r')
for line in fp:
    tmp=line.split()
    key=tmp[0]+"-"+tmp[1]+'-'+tmp[2]
    keylist.append(key)
    ratoms[key]=tmp[3:]
    
fp.close()

# read pdb
i=-1 # index counter
idxlist=[]
fp=open(pdb,'r')
gp=open("rigid_idx.txt",'w')
for line in fp:
    tmp=line.split()
    if tmp[0] == 'ATOM':
        i+=1
        key=tmp[10]+'-'+tmp[4]+'-'+tmp[3]
        if key in keylist:
            if tmp[2] in ratoms[key]:
                idxlist.append(i)
                print i,tmp[10],tmp[4],tmp[3],tmp[2]
                gp.write("%d %s %s %s %s\n" % (i,tmp[10],tmp[4],tmp[3],tmp[2]))
fp.close()
gp.close()

