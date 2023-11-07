import glob
import math
from sys import argv


#syntax -> python find_close_waters.py <ligand> <solvent_box> <output_solvent_box> 

lig_lines = [i.strip('\n').split() for i in open(argv[1], 'r').readlines()]
lig_coords = [i[5:8] for i in lig_lines if len(i) == 11]
solv_lines = [i.strip('\n').split() for i in open(argv[2], 'r').readlines()]
o = open(argv[3], 'w')




def rms(i, j):
    diff = float(i) - float(j)
    return diff ** 2

def length(x, y):
    list(x)
    list(y)
    return math.sqrt(rms(x[1], y[1]) + rms(x[2], y[2]) + rms(x[0], y[0]))


def make_pdb_line(atom):
    prt='ATOM'
    prt=prt+"{0: >7}".format(atom[1])
    if len(atom[2]) == 4:
        prt=prt+"{0: >5}".format(atom[2])
    else:
        prt=prt+"  "+"{0: <3}".format(atom[2])
    prt=prt+" "+"{0: <4}".format(atom[3])
    if int(atom[4]) <= 9999:
        prt=prt+"{0: >5}".format(atom[4])
        prt=prt+"{0: >12}".format('{:.3f}'.format(float(atom[5])))
    elif int(atom[4]) <= 99999 and int(atom[4]) > 9999:
        prt=prt+" "+"{0: >5}".format(atom[4])
        prt=prt+"{0: >11}".format('{:.3f}'.format(float(atom[5])))
    else:
        print('ERROR, RESIDUE OVERFLOW')
        exit(1)
    prt=prt+"{0: >8}".format('{:.3f}'.format(float(atom[6])))
    prt=prt+"{0: >8}".format('{:.3f}'.format(float(atom[7])))
    prt=prt+"  1.00"
    prt=prt+"  0.00"
    prt=prt+"{0: >10}".format(atom[10])
    prt=prt+'\n'
    return(prt)


water_res = None
mol_buffer = []
for line in solv_lines:
    if len(line) < 11:
        continue
    if line[4] != water_res: #write the previous water once we've gotten to the next one
        if len(mol_buffer) == 3: # only do so if all 3 atoms passed the distance test
            for mol_line in mol_buffer:
                o.write(make_pdb_line(mol_line))
        water_res = line[4]
        mol_buffer = []
    too_close_to_water = 0 # stopping 1 atom from being added stops the entire water from being added
    for coord in lig_coords:
        if length(coord, line[5:8]) < 3:
            too_close_to_water = 1
            break
    if too_close_to_water == 0:
        mol_buffer.append(line)

#print last water if possible
if len(mol_buffer) == 3:
    for mol_line in mol_buffer:
        o.write(make_pdb_line(mol_line))

o.write('END\n')
"""
docked_atoms = []

for pose in pose_list:
    lines = [i.strip('\n').split() for i in open(pose, 'r').readlines()]
    for line in lines:
        if line[0] == 'ATOM':
            docked_atoms.append(line[5:8])

prot_lines = [i.strip('\n').split() for i in open(protein, 'r').readlines()]

close_residues = []

for line in prot_lines:
    if line[0] == 'ATOM':
        res = line[4]
        coords = line[5:8]
        for atom in docked_atoms:
            if length(coords,atom) < 5:
                if res not in close_residues:
                    close_residues.append(res)

print(close_residues)

o = open(out_file, 'w')

o.write('set endvar = '+str(len(close_residues))+'\n')
o.write('\n')

for i,res in enumerate(close_residues):
    o.write('set var'+str(i+1)+' = '+str(res)+'\n')
"""
