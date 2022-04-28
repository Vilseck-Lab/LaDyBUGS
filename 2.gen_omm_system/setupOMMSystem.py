#! /usr/bin/env python

from simtk.openmm.app import PDBFile, CharmmPsfFile, CharmmParameterSet, CutoffPeriodic, HBonds
from simtk.openmm import NonbondedForce, XmlSerializer, CustomNonbondedForce
from simtk.unit import angstrom, nanometer, kilocalorie_per_mole, kilojoule_per_mole, elementary_charge, sqrt
from sys import stdout, exit, stderr
from glob import glob
#import numpy as np

def setup_NBFIXs(boxsize, prep_root, trunc=False):
    
    #### (1) Load Charmm PSF/PDB structure files & set up the system ####
    pdb = PDBFile('./patch.pdb')
    psf = CharmmPsfFile('./patch.psf')
    box=boxsize*angstrom
    psf.setBox(box,box,box,90,90,90)
    
    site2_check = len(glob(prep_root+'/site2*'))
    if site2_check != 0:
        raise ValueError('support for multi-site has not been developed yet!')

    # get list of all parameter files
    # hard coded charmm topology/parameter files
    param_list = [prep_root+'/toppar/top_all36_prot.rtf',prep_root+'/toppar/par_all36m_prot.prm',
                    prep_root+'/toppar/top_all36_na.rtf',prep_root+'/toppar/par_all36_na.prm',
                    prep_root+'/toppar/top_all36_carb.rtf',prep_root+'/toppar/par_all36_carb.prm',
                    prep_root+'/toppar/top_all36_cgenff.rtf',prep_root+'/toppar/par_all36_cgenff.prm',
                    prep_root+'/toppar/top_water.rtf',prep_root+'/toppar/par_water.prm']
   
    # use these lines for a generic set of parameter files
    #param_list = glob(prep_root+'/toppar/*.rtf')
    #param_list = param_list + glob(prep_root+'/toppar/*.prm')

    # ligand specific parameter files, doesn't support site2 quite yet
    param_list = param_list + [prep_root+'/core.rtf', prep_root+'/full_ligand.prm']
    param_list = param_list + glob(prep_root+'/site1_sub*_pres.rtf')

    params = CharmmParameterSet(*param_list)

    # create the system
    system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic,
                              solventDielectric=1.0, removeCMMotion=False,
                              nonbondedCutoff=1.2*nanometer, constraints=HBonds)

    # periodic boundary condition
    v1 = (box,0,0)
    v2 = (0,box,0)
    v3 = (0,0,box)
    system.setDefaultPeriodicBoxVectors(v1,v2,v3)

    #find site1_sub*.pdb files:
    sub_files = glob(prep_root+'/site1_sub*_frag.pdb')
    sub_files.sort() # need to ensure consistent order for fragments
    if len(sub_files) >= 10:
        raise ValueError('support for 10 fragments has not been developed yet!\n Maybe implement zero-padding in msld_tools')
    

    ## Set up the selections 
    nsites=1                                          # number of sites
    resid1=1                                          # resid number
    name=[[]]                                         # name[0] = environment (emtpy)
    sub2site = [0]
    # get atom names form sub_files, collect them into a list, append list to name
    for sub in sub_files:
        sub2site.append(1)
        atom_names = []
        with open(sub, 'r') as s:
            for line in s:
                if len(line.split()) > 9:
                    atom = line.split()[2] # the 3rd column of the pdb is the atom name
                    atom_names.append(atom.upper())
        name.append(atom_names)
    
    #alternative method of making sub2site (needs numpy)
    #sub2site = [0] 
    #sub2site = sub2site + list(np.ones(len(sub_files),np.int8))
    
    print('####### SITE NAMES #######')
    print(name)
    print('sub2site')
    print(sub2site)
    
    sele=[[] for i in range(0,len(name))]             #sele[0] = env indices
    sele_all=[]
    sele_env=[]
    sele_lam=[]

    for at in pdb.topology.atoms():
      found=0
      for i in range(1,len(name)):
        if at.name in name[i] and int(at.residue.id)==resid1:
          # print ("{} in {}".format(at.name,name[i]))
          sele[i].append(at.index)
          sele_all.append(at.index)
          sele_lam.append(at.index)
          found=1
        else:
          # print ("{} not in {}".format(at.name,name[i]))
          pass
      if found==0:
        sele[0].append(at.index)
        sele_all.append(at.index)
        sele_env.append(at.index)


    ## get list of atom types in nbfix_types
    nblist=[]
    for key in params.nbfix_types.keys():
        for i in range(len(key)):
            if key[i] not in nblist:
                nblist.append(key[i])

    ## get nonbonded parameters from NonbondedForce
    for force in system.getForces():
        if type(force) == NonbondedForce:
            nbf = force
            break
    epsilons = []
    sigmas = []
    charges = []
    sites = []
    subs = []
    for i in range(system.getNumParticles()):
        [q, sig, eps] = nbf.getParticleParameters(i)
        charges.append(q)
        sigmas.append(sig)
        epsilons.append(eps)
        foundcount=0
        for j in range(0,len(sele)):
          if i in sele[j]:
            subs.append(j)
            sites.append(sub2site[j])
            foundcount+=1
        if foundcount != 1:
          msg="found atom "+str(i)+" in "+str(foundcount)+" selections"
          print(msg)
          quit()

    # remove the NBFIX CustomNonbondedForce psf.createSystem makes
    for i in range(system.getNumForces()):
      if type(system.getForce(i)) == CustomNonbondedForce:
        system.removeForce(i)
        break

    ## make a Discrete2DFunction for vdW interactions
    ## (i) make a list of unique atom types, and assign each a unique "type"
    type_list = []
    type_num = {}
    type_params = {}
    i=0
    for at in psf.atom_list:
        if not (at.attype in type_list):
            type_list.append(at.attype)
            type_params[at.attype] = {'Rmin':at.type.rmin,'Eps':at.type.epsilon, \
                                      'Rmin_14':at.type.rmin_14,'Eps_14':at.type.epsilon_14}
            type_num[at.attype] = i
            i+=1

    ## (ii) generate 2DFunction w/o NBFIXs
    len_conv = angstrom.conversion_factor_to(nanometer)
    ene_conv = kilocalorie_per_mole.conversion_factor_to(kilojoule_per_mole)
    lj12a = [0 for i in range(len(type_list)*len(type_list))]
    lj06b = lj12a[:]
    for i in range(len(type_list)):
        for j in range(len(type_list)):
            Eps_ij = sqrt(type_params[type_list[i]]['Eps']*type_params[type_list[j]]['Eps'])*ene_conv
            Rmin_ij = (type_params[type_list[i]]['Rmin']+type_params[type_list[j]]['Rmin'])*len_conv
            lj12a[i+len(type_list)*j] = Eps_ij * Rmin_ij**12
            #lj12a[i+len(type_list)*j] = sqrt(Eps_ij) * Rmin_ij**6
            lj06b[i+len(type_list)*j] = 2 * Eps_ij * Rmin_ij**6

    ## (iii) back correct for NBFIX changes
    for at1 in range(len(psf.atom_list)):
        for at2 in range(len(psf.atom_list)):
            # does at1 have an nbfix?
            if len(psf.atom_list[at1].type.nbfix) > 0:
                # if yes, is at2 a part of that nbfix?
                if psf.atom_list[at2].attype in psf.atom_list[at1].type.nbfix.keys():
                    # if yes, update lj12a and lj06b lists
                    Rmin_nb=psf.atom_list[at1].type.nbfix[psf.atom_list[at2].attype][0]
                    Eps_nb=psf.atom_list[at1].type.nbfix[psf.atom_list[at2].attype][1]
                    type1 = type_num[psf.atom_list[at1].attype]
                    type2 = type_num[psf.atom_list[at2].attype]
                    lj12a[type1+len(type_list)*type2] = (Eps_nb*ene_conv) * (Rmin_nb*len_conv)**12
                    #lj12a[type1+len(type_list)*type2] = sqrt(Eps_nb*ene_conv) * (Rmin_nb*len_conv)**6
                    lj06b[type1+len(type_list)*type2] = 2 * (Eps_nb*ene_conv) * (Rmin_nb*len_conv)**6

    ## (iv) Repeat for 1,4 Interactions
    lj12a_14 = [0 for i in range(len(type_list)*len(type_list))]
    lj06b_14 = lj12a_14[:]
    for i in range(len(type_list)):
        for j in range(len(type_list)):
            Eps_ij = sqrt(type_params[type_list[i]]['Eps_14']*type_params[type_list[j]]['Eps_14'])*ene_conv
            Rmin_ij = (type_params[type_list[i]]['Rmin_14']+type_params[type_list[j]]['Rmin_14'])*len_conv
            lj12a_14[i+len(type_list)*j] = Eps_ij * Rmin_ij**12
            #lj12a_14[i+len(type_list)*j] = sqrt(Eps_ij) * Rmin_ij**6
            lj06b_14[i+len(type_list)*j] = 2 * Eps_ij * Rmin_ij**6

    ## back correct for NBFIX changes
    for at1 in range(len(psf.atom_list)):
        for at2 in range(len(psf.atom_list)):
            # does at1 have an nbfix?
            if len(psf.atom_list[at1].type.nbfix) > 0:
                # if yes, is at2 a part of that nbfix?
                if psf.atom_list[at2].attype in psf.atom_list[at1].type.nbfix.keys():
                    # if yes, update lj12a_14 and lj06b_14 lists
                    Rmin_nb=psf.atom_list[at1].type.nbfix[psf.atom_list[at2].attype][2]
                    Eps_nb=psf.atom_list[at1].type.nbfix[psf.atom_list[at2].attype][3]
                    type1 = type_num[psf.atom_list[at1].attype]
                    type2 = type_num[psf.atom_list[at2].attype]
                    lj12a_14[type1+len(type_list)*type2] = (Eps_nb*ene_conv) * (Rmin_nb*len_conv)**12
                    #lj12a_14[type1+len(type_list)*type2] = sqrt(Eps_nb*ene_conv) * (Rmin_nb*len_conv)**6
                    lj06b_14[type1+len(type_list)*type2] = 2 * (Eps_nb*ene_conv) * (Rmin_nb*len_conv)**6


    ## lj12a and lj06b should be ready for use as Discrete2DFunctions in the below CustomNonbondedForces


    ## Set up the nonbonded forces
    ## force group
    force_group_0 = 10  # environment
    force_group_1 = 20  # alchemical
    force_group_2 = 30  # restraints

    bondlist = []
    for bond in psf.bond_list:
        bondlist.append([bond.atom1.idx, bond.atom2.idx])

    # force equations, matches the CHARMM/OpenMM Interface (J. Comp. Chem. 2009, 30, 1545-1615)
    # environment
    formulaenv='(step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6)+step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3)-step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3)+step(Ron-r)*ch*(r1+eadd)+(step(r-Ron)*step(Roff-r)-step(r-Ron)*step(Ron-r))*ch*( r1* (acoef - s2*(bcoef + s2*(cover3 + dover5*s2))) +const));ch=138.935456*charge1*charge2;const  = bcoef*Roff-acoef/Roff+cover3*off3+dover5*off5;dover5 = dcoef/5.0;dcoef  = 2.0*denom;ccoef  = 3.0*cover3;cover3 = -(c2onnb+c2ofnb)*denom;acoef  = off4*(c2ofnb-3.0*c2onnb)*denom;bcoef  = 6.0*onoff2*denom;eadd   = (onoff2*(Roff-Ron)-(off5-on3*c2onnb)/5.0)*8.0*denom;denom  = 1.0/(c2ofnb-c2onnb)^3;off4   = c2ofnb*c2ofnb;off5   = off3*c2ofnb;onoff2 = c2onnb*c2ofnb;cr6  = ccnbb*ofdif3*rjunk3;cr12 = ccnba*ofdif6*rjunk6;rjunk3 =r3-recof3;rjunk6 = tr6-recof6;r3 = r1*tr2;r1 = sqrt(tr2);tr6 = tr2 * tr2 * tr2;tr2 = 1.0/s2;s2 = r*r;ccnbb = lj06b(type1,type2);ccnba = lj12a(type1,type2);onoff3 = recof3/on3;onoff6 = recof6/on6;ofdif3 = off3/(off3 - on3);ofdif6 = off6/(off6 - on6);recof3 = 1.0/off3;on6 = on3*on3;on3 = c2onnb*Ron;recof6 = 1.0/off6;off6 = off3*off3;off3 = c2ofnb*Roff;c2ofnb = Roff*Roff;c2onnb = Ron*Ron;'

    # soft core for alchemical atoms, matches function used for MSLD (J. Phys. Chem. B 2017, 121, 15, 3626–3635)
    softcore = "rp=step(r-rcsoft)*r+(1-step(r-rcsoft))*rcsoft*(-0.5*(r/(rcsoft+delta(1-scale)))^4+(r/(rcsoft+delta(1-scale)))^3+0.5);rcsoft = 2.0*sqrt(0.04)*(1.0-scale)"
    formulasc='(step(Ron-rp)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6)+step(rp-Ron)*step(Roff-rp)*(cr12*rjunk6 - cr6*rjunk3)-step(rp-Ron)*step(Ron-rp)*(cr12*rjunk6 - cr6*rjunk3)+step(Ron-rp)*ch*(r1+eadd)+(step(rp-Ron)*step(Roff-rp)-step(rp-Ron)*step(Ron-rp))*ch*( r1* (acoef - s2*(bcoef + s2*(cover3 + dover5*s2))) + const));ch=138.935456*charge1*charge2;const  = bcoef*Roff-acoef/Roff+cover3*off3+dover5*off5;dover5 = dcoef/5.0;dcoef  = 2.0*denom;ccoef  = 3.0*cover3;cover3 = -(c2onnb+c2ofnb)*denom;acoef  = off4*(c2ofnb-3.0*c2onnb)*denom;bcoef  = 6.0*onoff2*denom;eadd   = (onoff2*(Roff-Ron)-(off5-on3*c2onnb)/5.0)*8.0*denom;denom  = 1.0/(c2ofnb-c2onnb)^3;off4   = c2ofnb*c2ofnb;off5   = off3*c2ofnb;onoff2 = c2onnb*c2ofnb;cr6  = ccnbb*ofdif3*rjunk3;cr12 = ccnba*ofdif6*rjunk6;rjunk3 = r3-recof3;rjunk6 = tr6-recof6;r3 = r1*tr2;r1 = sqrt(tr2);tr6 = tr2 * tr2 * tr2;tr2 = 1.0/s2;s2 = rp*rp;ccnbb = lj06b(type1,type2);ccnba = lj12a(type1,type2);onoff3 = recof3/on3;onoff6 = recof6/on6;ofdif3 = off3/(off3 - on3);ofdif6 = off6/(off6 - on6);recof3 = 1.0/off3;on6 = on3*on3;on3 = c2onnb*Ron;recof6 = 1.0/off6;off6 = off3*off3;off3 = c2ofnb*Roff;c2ofnb = Roff*Roff;c2onnb = Ron*Ron;'

    # 1,4 Interactions (for both environment and alchemial atoms)
    formula14='(step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6)+step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3)-step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3)+step(Ron-r)*ch*(r1+eadd)+(step(r-Ron)*step(Roff-r)-step(r-Ron)*step(Ron-r))*ch*( r1* (acoef - s2*(bcoef + s2*(cover3 + dover5*s2))) + const));ch=138.935456*chargeprod;const  = bcoef*Roff-acoef/Roff+cover3*off3+dover5*off5;dover5 = dcoef/5.0;dcoef  = 2.0*denom;ccoef  = 3.0*cover3;cover3 = -(c2onnb+c2ofnb)*denom;acoef  = off4*(c2ofnb-3.0*c2onnb)*denom;bcoef  = 6.0*onoff2*denom;eadd   = (onoff2*(Roff-Ron)-(off5-on3*c2onnb)/5.0)*8.0*denom;denom  = 1.0/(c2ofnb-c2onnb)^3;off4   = c2ofnb*c2ofnb;off5   = off3*c2ofnb;onoff2 = c2onnb*c2ofnb;cr6  = ccnbb*ofdif3*rjunk3;cr12 = ccnba*ofdif6*rjunk6;rjunk3 = r3-recof3;rjunk6 = tr6-recof6;r3 = r1*tr2;r1 = sqrt(tr2);tr6 = tr2 * tr2 * tr2;tr2 = 1.0/s2;s2 = r*r;ccnbb = lj06b_14;ccnba = lj12a_14;onoff3 = recof3/on3;onoff6 = recof6/on6;ofdif3 = off3/(off3 - on3);ofdif6 = off6/(off6 - on6);recof3 = 1.0/off3;on6 = on3*on3;on3 = c2onnb*Ron;recof6 = 1.0/off6;off6 = off3*off3;off3 =c2ofnb*Roff;c2ofnb = Roff*Roff;c2onnb = Ron*Ron;'

    # define lambda scaling
    scaling1 = "scale1=step(0.5-abs(subs1-0))"
    scaling2 = "scale2=step(0.5-abs(subs2-0))"
    for i in range(1,len(sele)):
      scaling1=scaling1+"+step(0.5-abs(subs1-"+str(i)+"))*lambda"+str(i)
      scaling2=scaling2+"+step(0.5-abs(subs2-"+str(i)+"))*lambda"+str(i)
    scaling=scaling1+";"+scaling2
    scaling="samesub=step(0.5-abs(subs1-subs2));"+scaling
    scaling="samesite=step(0.5-abs(sites1-sites2));"+scaling
    scaling="diffsite=step(abs(sites1-sites2)-0.5);"+scaling
    scaling="scale=samesite*samesub*scale1+diffsite*scale1*scale2;"+scaling

    #### Setup Custom Nonbonded Forces ####
    # Env-Env
    nonbond=[]
    nonbond.append(CustomNonbondedForce(formulaenv))
    nonbond[0].addTabulatedFunction('lj12a',Discrete2DFunction(len(type_list),len(type_list),lj12a))
    nonbond[0].addTabulatedFunction('lj06b',Discrete2DFunction(len(type_list),len(type_list),lj06b))
    nonbond[0].addPerParticleParameter('charge')
    nonbond[0].addPerParticleParameter('type')
    nonbond[0].addGlobalParameter("Roff", 1.2)
    nonbond[0].addGlobalParameter("Ron", 1.0)
    for i in range(len(psf.atom_list)):
        nonbond[0].addParticle([charges[i],type_num[psf.atom_list[i].attype]])
    nonbond[0].addInteractionGroup(sele_env,sele_env)
    nonbond[0].setUseSwitchingFunction(False)
    nonbond[0].setCutoffDistance(1.2*nanometer)
    nonbond[0].createExclusionsFromBonds(bondlist, 3)
    excl={}
    for i in range(0,nonbond[0].getNumExclusions()):    
      tmp = nonbond[0].getExclusionParticles(i)         
      for j in range(0,2):
        if tmp[j] in excl:
          excl[tmp[j]].append(tmp[1-j])
        else:
          excl[tmp[j]]=[tmp[1-j]]
    for i in range(1,len(sele)):
      for j in range(i+1,len(sele)):
        if i!=j and sub2site[i]==sub2site[j]:
          for ii in sele[i]:
            for jj in sele[j]:
              if not (ii in excl[jj]):
                nonbond[0].addExclusion(ii,jj)
    #           else:
    #             print("Redundant exclusion: "+str(ii)+" "+str(jj))
    nonbond[0].setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    nonbond[0].setForceGroup(force_group_0)
    system.addForce(nonbond[0])

    # Alchemical Atoms/Groups
    nonbond.append(CustomNonbondedForce("scale * "+formulasc+";"+softcore+";"+scaling))
    nonbond[1].addTabulatedFunction('lj12a',Discrete2DFunction(len(type_list),len(type_list),lj12a))
    nonbond[1].addTabulatedFunction('lj06b',Discrete2DFunction(len(type_list),len(type_list),lj06b))
    nonbond[1].addPerParticleParameter('charge')
    nonbond[1].addPerParticleParameter('type')
    nonbond[1].addPerParticleParameter('subs')
    nonbond[1].addPerParticleParameter('sites')
    nonbond[1].addGlobalParameter("Roff", 1.2) 
    nonbond[1].addGlobalParameter("Ron", 1.0)
    for i in range(1,len(sele)):
      nonbond[1].addGlobalParameter("lambda"+str(i), 0.000)
    for i in range(len(psf.atom_list)):
        nonbond[1].addParticle([charges[i],type_num[psf.atom_list[i].attype],subs[i],sites[i]])
    nonbond[1].addInteractionGroup(sele_lam,sele_all)
    nonbond[1].setUseSwitchingFunction(False)
    nonbond[1].setCutoffDistance(1.2*nanometer)
    nonbond[1].createExclusionsFromBonds(bondlist, 3)
    excl={}
    for i in range(0,nonbond[1].getNumExclusions()):
      tmp = nonbond[1].getExclusionParticles(i)
      for j in range(0,2):
        if tmp[j] in excl:
          excl[tmp[j]].append(tmp[1-j])
        else:
          excl[tmp[j]]=[tmp[1-j]]
    for i in range(1,len(sele)):
      for j in range(i+1,len(sele)):
        if i!=j and sub2site[i]==sub2site[j]:
          for ii in sele[i]:
            for jj in sele[j]:
              if not (ii in excl[jj]):
                nonbond[1].addExclusion(ii,jj)
    #           else:
    #             print("Redundant exclusion: "+str(ii)+" "+str(jj))
    nonbond[1].setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    nonbond[1].setForceGroup(force_group_1)
    system.addForce(nonbond[1])

    # Env-Env 1,4-Interactions
    nonbond14=[]
    nonbond14.append(CustomBondForce(formula14))
    nonbond14[0].addPerBondParameter("chargeprod")
    nonbond14[0].addPerBondParameter("lj12a_14")
    nonbond14[0].addPerBondParameter("lj06b_14")
    nonbond14[0].addGlobalParameter("Roff", 1.2)
    nonbond14[0].addGlobalParameter("Ron", 1.0)
    for i in range(nbf.getNumExceptions()):
        tmp = nbf.getExceptionParameters(i)
        if tmp[0] in sele_env and tmp[1] in sele_env and \
           (tmp[2].value_in_unit(elementary_charge**2) != 0 or tmp[4].value_in_unit(kilojoule/mole) != 0):
            # replace epsilons and sigmas with types
            type1=type_num[psf.atom_list[tmp[0]].attype]
            type2=type_num[psf.atom_list[tmp[1]].attype]
            nonbond14[0].addBond(tmp[0], tmp[1], [tmp[2].value_in_unit(elementary_charge**2), \
                                 lj12a_14[type1+len(type_list)*type2], lj06b_14[type1+len(type_list)*type2]])
    nonbond14[0].setForceGroup(force_group_0)
    #nonbond14[0].setForceGroup(25)
    system.addForce(nonbond14[0])

    # Alchemial 1,4-Interactions
    nonbond14.append(CustomBondForce("scale * "+formula14+";"+scaling))
    nonbond14[1].addPerBondParameter("chargeprod")
    nonbond14[1].addPerBondParameter("lj12a_14")
    nonbond14[1].addPerBondParameter("lj06b_14")
    nonbond14[1].addPerBondParameter("subs1")
    nonbond14[1].addPerBondParameter("subs2")
    nonbond14[1].addPerBondParameter("sites1")
    nonbond14[1].addPerBondParameter("sites2")
    nonbond14[1].addGlobalParameter("Roff", 1.2)
    nonbond14[1].addGlobalParameter("Ron", 1.0)
    for i in range(1,len(sele)):
      nonbond14[1].addGlobalParameter("lambda"+str(i), 0.000)
    for i in range(nbf.getNumExceptions()):
        tmp = nbf.getExceptionParameters(i)
        if (tmp[0] in sele_lam or tmp[1] in sele_lam) and \
           (tmp[2].value_in_unit(elementary_charge**2) != 0 or tmp[4].value_in_unit(kilojoule/mole) != 0):
            type1=type_num[psf.atom_list[tmp[0]].attype]
            type2=type_num[psf.atom_list[tmp[1]].attype]
            nonbond14[1].addBond(tmp[0], tmp[1], [tmp[2].value_in_unit(elementary_charge**2), \
                                 lj12a_14[type1+len(type_list)*type2], lj06b_14[type1+len(type_list)*type2], \
                                 subs[tmp[0]],subs[tmp[1]],sites[tmp[0]],sites[tmp[1]]])
    nonbond14[1].setForceGroup(force_group_1)
    #nonbond14[1].setForceGroup(27)
    system.addForce(nonbond14[1])


    # remove original OMM Nonbonded Force
    for i in range(system.getNumForces()):
      if type(system.getForce(i)) == NonbondedForce:
        system.removeForce(i)
        break

    if trunc==True:
        ## Add harmonic restraints in for truncated peptide caps
        # read in rigid_atom_indices
        rigid_idx=[]
        rigid_dict={}
        fp=open('rigid_idx.txt','r')   # user generated list of the format "ID# SEGID RESID RESN ATOMNAME" (ex. "0 PA01 16 ASN N")
        for line in fp:
            tmp=line.split()
            rigid_idx.append(int(tmp[0]))
            rigid_dict[int(tmp[0])]={'SEGID':tmp[1],'RESID':int(tmp[2]),'RESN':tmp[3],'ATNAME':tmp[4]}
        fp.close()
        # make the harmonic restraint
        cons_harm=CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
        cons_harm.addPerParticleParameter('k')
        cons_harm.addPerParticleParameter('x0')
        cons_harm.addPerParticleParameter('y0')
        cons_harm.addPerParticleParameter('z0')
        for i in rigid_idx:
            # check for correct idx match
            if ((psf.atom_list[i].name == rigid_dict[i]['ATNAME']) and \
                (psf.atom_list[i].system == rigid_dict[i]['SEGID']) and \
                (psf.atom_list[i].residue.resname == rigid_dict[i]['RESN']) and \
                (psf.atom_list[i].residue.idx == rigid_dict[i]['RESID'])):
                # add the particle
                k=10*kilocalorie/(mole*angstrom*angstrom)  # Similar to what is used in CHARMM
                k=k*psf.atom_list[i].mass._value           # mass weighted (without the mass unit)
                k.value_in_unit(kilojoule/(mole*nanometer*nanometer))  # OMM units
                x0=pdb.positions[i][0]
                y0=pdb.positions[i][1]
                z0=pdb.positions[i][2]
                cons_harm.addParticle(i,[k,x0,y0,z0])
            else:
                # throw an error
                print("Error in matching props for",i)
                if psf.atom_list[i].name != rigid_dict[i]['ATNAME']:
                    print(psf.atom_list[i].name, ' DOES NOT MATCH ', rigid_dict[i]['ATNAME'])
                if psf.atom_list[i].system != rigid_dict[i]['SEGID']:
                    print(psf.atom_list[i].system, ' DOES NOT MATCH ', rigid_dict[i]['SEGID'])
                if psf.atom_list[i].residue.resname != rigid_dict[i]['RESN']:
                    print(psf.atom_list[i].residue.resname, ' DOES NOT MATCH ', rigid_dict[i]['RESN'])
                if psf.atom_list[i].residue.idx != rigid_dict[i]['RESID']:
                    print(psf.atom_list[i].residue.idx, ' DOES NOT MATCH ', rigid_dict[i]['RESID'])

        cons_harm.setForceGroup(force_group_2)
        system.addForce(cons_harm)


    ## save the system 
    xml = XmlSerializer.serialize(system)
    f = open("./omm_system.xml", 'w')
    f.write(xml)
    f.close()

    print("system set up")

   
if __name__ == "__main__":
    setup_NBFIXs(83, 'prep/', trunc=True)
