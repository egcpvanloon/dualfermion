################################################################################
#
# dualfermion: second-order dual fermion perturbation theory
#
# Copyright (C) 2020, E.G.C.P. van Loon
# 
# Based on TRIQS: a Toolbox for Research in Interacting Quantum Systems
# 
# dualfermion is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# dualfermion is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# dualfermion. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

"""
Test in the atomic limit
A system consisting of two dimers (4 sites) is solved in two ways
1) ED using Pomerol on the 4-site cluster
2) ED using Pomerol on the 2-site dimer + dual fermion perturbation theory
In the limit of small interdimer hopping t1 (compared to U), 
method (2) should agree well with (1), since DF is exact to second-order in t1/U
"""


import dualfermion

import pytriqs.gf as gf

import pytriqs.operators as ops
import pytriqs.operators.util as op
import pytriqs.utility.mpi as mpi
import pytriqs.archive as ar
from pytriqs.archive import HDFArchive
from itertools import product

from pomerol2triqs import PomerolED

import numpy as np

dptkeys = ['ksi_delta','t1','verbosity','calculate_sigma']

parms = {
        # Solver parameters
        'n_iw' : 20,
        # Physical parameters
        'U'   : 10.,
        'J'   : 0.0,
        't1'  : 0.01,
        't0'  : 4.,
        'beta': 20.,
        'calculate_sigma' : True,
        'measure_G2_iw_ph': True,
        "measure_G2_n_bosonic": 9,
        "measure_G2_n_fermionic": 9,
        "verbosity":2,
        }

parms["N_x"]=2
parms["N_y"]=1
parms["N_z"]=1
parms["ksi_delta"]=1.0

# Chemical potential depends on the definition of H(k) that is used
parms['chemical_potential'] = -3. 

n_orbs = 2                             # Number of orbitals
off_diag=True
spin_names = ['up','dn']                   # Outer (non-hybridizing) blocks
orb_names = ['%s'%i for i in range(n_orbs)]  # Orbital indices
gf_struct = op.set_operator_structure(spin_names,orb_names,off_diag=off_diag) 
#Ntot = sum(ops.n(sigma,i) for sigma,i in product(spin_names,range(n_orbs)))

# Quadratic part of the atomic Hamiltonian
t_intra_atomic = parms["t0"] * np.array([[0,1],[1,0]])

# Construct the DF2 program
X = dualfermion.Dpt(beta=parms['beta'],gf_struct=gf_struct,n_iw=parms['n_iw'],n_iw2=parms["measure_G2_n_fermionic"],n_iW =parms["measure_G2_n_bosonic"],N_x=parms['N_x'],N_y=parms['N_y'],N_z=parms['N_z'])

def Hk_f(k):
    return -2 *parms['t1'] * (np.cos(k[0]))*np.eye(2)

for spin,_ in X.Hk:
    for k in X.Hk.mesh:
        X.Hk[spin][k] = Hk_f(k.value)

g2_blocks = set([("up", "up"), ("up", "dn"), ("dn", "up"),("dn","dn")])
index_converter = {(sn, o) : ("loc", int(o), "down" if sn == "dn" else "up")
                   for sn, o in product(spin_names, orb_names)}

#print index_converter,type(index_converter)
#raise Exception("")

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)
N = sum(ops.n(sn, o) for sn, o in product(spin_names, orb_names))
H = (
            parms["U"] * (ops.n('up','0') * ops.n('dn','0') )
            +parms['t0']*(
                  ops.c_dag('up','0')*ops.c('up','1') + ops.c_dag('up','1')*ops.c('up','0')
                 +ops.c_dag('dn','0')*ops.c('dn','1') + ops.c_dag('dn','1')*ops.c('dn','0')
                 )
            - parms['chemical_potential']*N
   )

#####
#
# Reference: 4 site cluster, calculate only G, not G2
#
#####
def calc_reference():

    ref_orbs = ['%s'%i for i in range(n_orbs*parms['N_x'])]
    ref_gf_struct = op.set_operator_structure(spin_names,ref_orbs,off_diag=off_diag) 
    ref_index_converter = {(sn, o) : ("loc", int(o), "down" if sn == "dn" else "up")
                   for sn, o in product(spin_names, ref_orbs)}
    #print ref_index_converter,ref_orbs    
    ref_ed = PomerolED(ref_index_converter, verbose = True)
    ref_N =sum(ops.n(sn, o) for sn, o in product(spin_names, ref_orbs))
    ref_H = (
            parms["U"] * (ops.n('up','0') * ops.n('dn','0')+ops.n('up','2') * ops.n('dn','2') )
            +parms['t0']*(
                  ops.c_dag('up','0')*ops.c('up','1') + ops.c_dag('up','1')*ops.c('up','0')
                 +ops.c_dag('dn','0')*ops.c('dn','1') + ops.c_dag('dn','1')*ops.c('dn','0')
                 +ops.c_dag('up','2')*ops.c('up','3') + ops.c_dag('up','3')*ops.c('up','2')
                 +ops.c_dag('dn','2')*ops.c('dn','3') + ops.c_dag('dn','3')*ops.c('dn','2')   
                 )
            -2.*parms['t1']*(
                  ops.c_dag('up','0')*ops.c('up','2') + ops.c_dag('up','2')*ops.c('up','0')
                 +ops.c_dag('dn','0')*ops.c('dn','2') + ops.c_dag('dn','2')*ops.c('dn','0')
                 +ops.c_dag('up','1')*ops.c('up','3') + ops.c_dag('up','3')*ops.c('up','1')
                 +ops.c_dag('dn','1')*ops.c('dn','3') + ops.c_dag('dn','3')*ops.c('dn','1') 
                 )
            - parms['chemical_potential']*ref_N
            )
    # Run the solver
    ref_ed.diagonalize(ref_H)
    # Compute G(i\omega)
    ref_G_iw = ref_ed.G_iw(ref_gf_struct, parms['beta'], parms['n_iw'])
    return ref_G_iw
            
ref_G_iw = calc_reference()
ref = ref_G_iw['up']


# Initial hybridization
X.Delta << 0.

# Obtain bath sites from Delta and create H_ED
H_ED = H
        
# Run the solver
ed.diagonalize(H_ED)
# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, parms['beta'], parms['n_iw'])

if parms["measure_G2_iw_ph"]:
    common_g2_params = {'gf_struct' : gf_struct,
                'beta' : parms['beta'],
                'blocks' : g2_blocks,
                'n_iw' : parms['measure_G2_n_bosonic']}
    G2_iw = ed.G2_iw_inu_inup(channel = "PH",
                                block_order = "AABB",
                                n_inu = parms['measure_G2_n_fermionic'],
                                **common_g2_params)
    
    X.G2_iw << G2_iw

    
# Some intermediate saves
gimp = G_iw['up']    

if mpi.is_master_node():

    # Run the dual perturbation theory
    X.gimp << G_iw # Load G from impurity solver
    dpt_parms = {key:parms[key] for key in parms if key in dptkeys}
    X.run(**dpt_parms)


    G_k = HDFArchive("G_k.h5",'r')['G_k']['up']


    inv = np.linalg.inv
    t1 = parms['t1']
    t0 = parms['t0']
    # Calculate lattice self-energy
    for ki,k in enumerate(G_k.mesh[1]):
        selfenergylist = []
        ref_list = []
        simp_list = []

        kx = k.value[0]
        tkx= -2*np.cos(kx)*t1
        for wn in G_k.mesh[0]:
            wval = wn.value.imag
            G = G_k[(wn,k)]
            
            G0k  = inv(np.matrix([[wn.value+parms['chemical_potential'] -tkx,-t0],[-t0,wn.value+parms['chemical_potential']-tkx]]))
            invG0k  = np.matrix([[wn.value+parms['chemical_potential'] -tkx,-t0],[-t0,wn.value+parms['chemical_potential']-tkx]])
            
            S = invG0k- inv(G) 
            selfenergylist.append( [wval,]+list(S.A.flatten()) )
            
            # Reference data, convert to 2x2 matrix with correct k
            sign = 1-2*ki
            ref_G = ref[(wn)][0:2,0:2] + sign * ref[(wn)][0:2,2:4]
            ref_S = invG0k - inv(ref_G)
            ref_list.append( [wval,]+list(ref_S.A.flatten()) )
            
            # calculate impurity self-energy as well
            # This was done using t1=0, so tkx=0
            g = gimp[(wn)]
            invg0= np.matrix([[wn.value+parms['chemical_potential'],-t0],[-t0,wn.value+parms['chemical_potential']]])
            Simp = invg0 - inv(g)
            simp_list.append( [wval,]+list(Simp.A.flatten()) )
            
            
        s = np.array(selfenergylist)
        r = np.array(ref_list)
        s2= np.array(simp_list)
        
        for column in range(1,5):
            myDat=(s[:,column]-s2[:,column]).real # Dual Fermion data
            myRef=(r[:,column]-s2[:,column]).real # Reference data
        
            assert np.max(np.abs(myRef-myDat)) < 1e-3
