__author__ = "JV"
__date__ = "2023"

import numpy as np
from simtk.openmm import PeriodicTorsionForce

def get_TorsionForce(system):
    """ pull out PeriodicTorsionForce for easier reference """
    for force in system.getForces(): 
        if type(force) == PeriodicTorsionForce:
            phif = force # identify the torsion force object
            break
    return phif

def read_phiIDK(nsubs):
    """ read in phiIDK.txt (assumes nsubs is already defined) """
    phiIDK=[]
    for s in range(nsubs[0]): # initialize variables ** ASSUMES 1 SITE
        phiIDK.append({})
    for line in open('phiIDK.txt','r'):
        tmp=line.split()
        phiIDK[int(tmp[0])][int(tmp[1])]=float(tmp[2])
    return phiIDK

def scale_phis(nsubs,Lrow,Lstates,phif,phiIDK):
    """ scale phi dihedral angles by lambda(i)
        !!! Assumes we only have 1 SITE (not multisite) !!!
        in most cases Lrow corresponds to res['idx']
    """
    jbuff=0
    #for site in range(len(nsubs)): # to make this MS would require site & sub definitions in phiIDK.txt
    #    for sub in range(nsubs[site]): # and site+sub nested lists/dicts in the phiIDK variable
    for sub in range(nsubs[0]):
        for k in phiIDK[sub].keys(): # the keys are the phi indices
            scaledk = phiIDK[sub][k]*Lstates[Lrow,jbuff]
            bphi = phif.getTorsionParameters(k)
            phif.setTorsionParameters(k,bphi[0],bphi[1],bphi[2],bphi[3],bphi[4],bphi[5],scaledk)
        jbuff+=1
    return


