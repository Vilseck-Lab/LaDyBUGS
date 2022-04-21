import os
from fe_setup import functions_firstStep as ffs  
from fe_setup import overlay_mcs 
#import msld_chk
#from charmm_core_align import charmm_align
from fe_setup.charmm_align import charmm_align

from fe_setup.msld_tools import msld_chk
from fe_setup.msld_tools import hacked_msld_mcs
from fe_setup.msld_tools import msld_crn
from fe_setup.msld_tools import msld_prm
from fe_setup.msld_tools import msld_lig_str
import glob

from fe_setup.find_close_waters import find_close_waters
#from setupSystem_NBFIXs import setup_NBFIXs
import re

root_name = 'lig_'
CHARMM_ALIGN = False
genpsf_file = 'genpsf.inp'

# for mol2 creation
#FML, how many mol2 generators do I need?
SPORES_exec = ffs.set_exec('/nfs/home/mrobo/SPORES_64bit')
SPORES_exec = SPORES_exec + ' --mode settypes '
unicon_exec = ffs.set_exec('/nfs/home/mrobo/unicon_1.4.0/unicon')
antechamber_exec = ffs.set_exec('/nfs/home/mrobo/anaconda3/envs/amber/bin/antechamber')
babel_exec = ffs.set_exec('/nfs/home/mrobo/openbabel-master/install/bin/obabel')

#for atom typing
cgenff_exec = ffs.set_exec('/nfs/home/mrobo/silcsbio.2020.1/cgenff/cgenff')

#use only if using cgenff w/ MATCH charges
MATCH_exec = ffs.set_exec('/nfs/home/mrobo/MATCH_RELEASE/MATCH/scripts/MATCH.pl')
MATCH_exec = MATCH_exec + ' -forcefield top_all36_cgenff_new '


#read in the pdb of the reference ligand
#ref_ligand = 'ex_00.pdb'
ref_ligand = 'start_struct.pdb'
with open(ref_ligand, 'r') as r_file:
    ref_lines = [i.strip('\n').split() for i in r_file.readlines()]

inp_pdbs = [ref_lines]

#step1: generate pdbs from smiles
new_smiles_file = './new_mols.txt'
if os.path.exists(new_smiles_file):
    new_mols_file = new_smiles_file
    with open(new_mols_file, 'r') as nm_file:
        new_mols = [i.strip('\n').split()[0] for i in nm_file.readlines()]
        for mol in new_mols:
            pdb_lines = ffs.gen_pdb_block(mol)
            pdb_lines = pdb_lines.split('\n')
            pdb_lines = [i.split() for i in pdb_lines]
            inp_pdbs.append(pdb_lines)

#read in additional pdbs
#pdb_list = ['ex_01.pdb', 'ex_02.pdb', 'ex_03.pdb', 'ex_04.pdb', 'ex_05.pdb']
pdb_list = ['p-tail_Cl.pdb',  'p-tail_F.pdb',  'p-tail_N.pdb',  'p-tail_OC.pdb',  'p-tail_O.pdb']
if len(pdb_list) > 0:
    for pdb in pdb_list:
        with open(pdb, 'r') as pdb_file:
            pdb_lines = [i.strip('\n').split() for i in pdb_file.readlines()]
            inp_pdbs.append(pdb_lines)

ele_counts = {}

molecules = [] 
"""
a list of dictionaries providing the identifiers for each molecule
# dictionary = {
'res_name' = ...
'seg_name' = ...
'base_name' = ...
'formatted_pdb_lines' = ...
}
"""


#  test_0 will be the reference ligand
# step 2: uniquely label atoms, then create pdb, mol2, prm, rtf, and str files
for i,pdb_lines in enumerate(inp_pdbs):
    new_dict = {}
    pad_num = str(i).zfill(2)
    res_name = 'L'+pad_num
    seg_name = 'LG'+pad_num
    base_name = root_name+pad_num
    
    # generate pdb, ger updated ele dict
    new_pdb_lines, ele_counts = ffs.renumber_pdb(pdb_lines, base_name, 
                    res_name, seg_name, ele_counts) 
    ffs.write_pdb(new_pdb_lines, base_name)
    
    #append new dictionary
    fields = ['res_name', 'seg_name', 'base_name', 'formatted_pdb_lines']
    values = [res_name, seg_name, base_name, new_pdb_lines] 
    new_dict = dict(zip(fields, values))
    molecules.append(new_dict)

#step 3: create mol2 files, from there create prm, rtf, and str files
for mol_dict in molecules:
    base_name = mol_dict['base_name']
    res_name = mol_dict['res_name']
    #MOL2 generation SUCKS!
    try:
        ffs.gen_antechamber_mol2(base_name, antechamber_exec)
        ffs.cgenff_MATCH_type(base_name, res_name, cgenff_exec, MATCH_exec)
    except:
        try:
            ffs.gen_SPORES_mol2(base_name, SPORES_exec)
            ffs.cgenff_MATCH_type(base_name, res_name, cgenff_exec, MATCH_exec)
        except:
            try:
                ffs.gen_unicon_mol2(base_name, unicon_exec)
                ffs.cgenff_MATCH_type(base_name, res_name, cgenff_exec, MATCH_exec)
            except:
                try:
                    ffs.gen_babel_mol2(base_name, babel_exec)
                    ffs.cgenff_MATCH_type(base_name, res_name, cgenff_exec, MATCH_exec)
                except:
                    raise SystemError('Atom typing failed!. Check how each method failed individually!')
    
#step 4: align the generated ligands (lig_01, lig_02 ...) using CHARMM
if CHARMM_ALIGN == True:
    for mol_dict in molecules[1:]:
        ref_name = molecules[0]['base_name']
        ref_lig_lines = molecules[0]['formatted_pdb_lines']
        gen_name = mol_dict['base_name']
        gen_lig_lines = mol_dict['formatted_pdb_lines']
        #step 3: Find the core atoms, use them to align the ligands in CHARMM
        mol_list = [ref_name, gen_name] # ref ligand is root_name+'00'
        core_atom_correspondence = overlay_mcs.MsldMCS(mol_list, groupRingFrags=True, verbose=False)
        print('providing core!')
        print(core_atom_correspondence)

        iterator = 0
        aligned_list = [] # we can get more than one MCS patt match in one (or both) ligands
        for ref_core in core_atom_correspondence[0]:     # to distinguish between all combinations
            for gen_core in core_atom_correspondence[1]: # we align them in CHARMM and get the lower energy pair
                attempt = charmm_align(ref_core, gen_core, ref_lig_lines, gen_lig_lines, ref_name, gen_name, str(iterator))
                iterator = iterator + 1
                aligned_list.append(attempt) # attempt = [energy, iterator, name_of_aligned_pdb]

        aligned_list.sort(key=lambda x:x[0]) # get alignment that had lowest energy in CHARMM
        best_mcs_pdb = aligned_list[0][2]    # get pdb of best alignment
        
        best_mcs_it = aligned_list[0][1]     # get MCS combination of best alignment
        best_ref_core = best_mcs_it // len(core_atom_correspondence[1]) # for a system with 2 ref cores and 3 gen cores
        best_ref_core = core_atom_correspondence[0][best_ref_core]      # interators 0 1 2 will correspond to a 
        best_gen_core = best_mcs_it % len(core_atom_correspondence[1])  # single ref core, and iterators 2 3   
        best_gen_core = core_atom_correspondence[1][best_gen_core]      # will correspond to a single gen core 
        
        os.rename(base_name+'.pdb', 'orig_'+base_name+'.pdb')
        os.rename(best_mcs_pdb, base_name+'.pdb')   
        
        #gen_spores_mol2(base_name, SPORES_exec)
        #gen_unicon_mol2(base_name, unicon_exec)
        ffs.gen_antechamber_mol2(base_name, antechamber_exec)

#step 5: use the libraries from py_prep

sysname = "ship1"                              # name of future output files
mol_list = [i['base_name'] for i in molecules] # list of mol2 file names
mcsout = 'MCS_for_MSLD.txt'               # MCS output filename
outdir = 'prep/'                 # MSLD output directory
otherdir = 'other/'
cgenff=True                               # Are CGenFF/ParamChem parameters being used?

#use these lists of lists only if you need to explicitly include atoms in Frag or Core
#inFrag=[[],['N001', 'C008'],['N004','C00M'],['N007', 'C019'],['N00A', 'C01W'], ['N00E', 'C02I'],['N00J','C034']]
#inFrag = [[], ['C008', 'O000'], ['C00M', 'O002']]
#inFrag = [[], ['O000'], ['O002']]
inFrag = [[]]
inCore=[[]]

## (5.1) Check molfile and toppar files before getting started
msld_chk.MsldCHK(mol_list)
reflig = hacked_msld_mcs.MsldMCS(mol_list,mcsout,inFrag,inCore,cutoff=1.0,debug=False)

## (5.2) Perform Charge-Renormalization
## To manually change what atoms are in the core/fragments, manually change mcsout
## or use "inFrag" and "inCore" nested lists above. "Anchor atoms" (and connected Hs)
## are automatically included in each fragment unless specifically stated to be "inCore"
msld_crn.MsldCRN(mcsout,outdir,inFrag,inCore,ChkQChange=False,verbose=False,debug=False)

## (5.3) Write Ligand Parameters & the lig-specific.str file
msld_prm.MsldPRM(outdir,cgenff,verbose=False,debug=False)
msld_lig_str.write_str(sysname,outdir,cgenff)
#exit()
#step 6: get all coordinates of all atoms in the ligand and then
#delete any waters that collide with any ligand atoms
all_lig_coords = []
all_ligand_files = glob.glob(outdir+'/site*_sub*.pdb')
all_ligand_files.append(outdir+'core.pdb')
for pdb_file in all_ligand_files:
    lig_lines = [i.strip('\n').split() for i in open(pdb_file, 'r').readlines()]
#    lig_coords = [i[5:8] for i in lig_lines if len(i) == 10]
    lig_coords = [[i[2],float(i[5]),float(i[6]),float(i[7])] for i in lig_lines if len(i) == 11] #need atom name for consistency
    for coord in lig_coords:
        all_lig_coords.append(coord)
find_close_waters(all_lig_coords, 'solv_10.0.pdb', outdir+'fitted-solv_10.0.pdb')

#step 7: generate the patch.{pdb,psf} files using CHARMM
#make sure that trunc.inp has all of the information (solvent, protein, ions) that you need
#exit()
#step 8: generate the omm_patch.xml file 
# grab box size from genpsf_file
with open(genpsf_file, 'r') as gpsf_file:
    for line in gpsf_file:
        regex_val = re.search('set box =', line)
        if regex_val != None:
            boxline = line
            break
boxsize = float(boxline.split()[3])

#run setupNBFIXEs
#This function will read the files inside of outdir, and create the omm_patch.xml file
#setup_NBFIXs(boxsize, outdir)
