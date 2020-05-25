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

ref_spin = 'up' # used for the comparison

# with bath_prefactor = 1, the bath is switched on
# with bath_prefactor = 0, the bath is switched off
bath_prefactor =1

def run_test(t1,filename):
    dptkeys = ['verbosity','calculate_sigma','calculate_sigma1','calculate_sigma2']

    parms = {
            # Solver parameters
            'n_iw' : 100,
            # Physical parameters
            'U'   : 0.5,
            't1'  : t1,
            'beta': 10,
            # DMFT loop control parameters
            'calculate_sigma' : True,
            'calculate_sigma1' : True,
            'calculate_sigma2' : True,
            'measure_G2_iw_ph': True,
            "measure_G2_n_bosonic": 10,
            "measure_G2_n_fermionic": 10,
            "verbosity":4,
            }

    parms["N_x"]=2
    parms["N_y"]=1
    parms["N_z"]=1
    parms["ksi_delta"]=1.0

    # Chemical potential depends on the definition of H(k) that is used
    parms['chemical_potential_bare'] = 0.
    parms['chemical_potential'] = parms['U']/2.+parms['chemical_potential_bare']

    n_orbs = 1                             # Number of orbitals
    off_diag=True
    spin_names = ['up','dn']                   # Outer (non-hybridizing) blocks
    orb_names = ['%s'%i for i in range(n_orbs)]  # Orbital indices
    gf_struct = op.set_operator_structure(spin_names,orb_names,off_diag=off_diag) 

    #####
    #
    # Reference: 4 site cluster, calculate only G, not G2
    #
    #####
    def calc_reference():

        ref_orbs = ['%s'%i for i in range(n_orbs*parms['N_x']*parms['N_y'])]
        ref_gf_struct = op.set_operator_structure(spin_names,ref_orbs,off_diag=off_diag) 
        ref_index_converter = {(sn, o) : ("loc", int(o), "down" if sn == "dn" else "up")
                    for sn, o in product(spin_names, ref_orbs)}
        print ref_index_converter,ref_orbs    
        ref_ed = PomerolED(ref_index_converter, verbose = True)
        ref_N =sum(ops.n(sn, o) for sn, o in product(spin_names, ref_orbs))
        #  2 3
        #  0 1
        ref_H = (
                parms["U"] * (
                    ops.n('up','0') * ops.n('dn','0')
                    +ops.n('up','1') * ops.n('dn','1')
                    )
                -2.*parms['t1']*(
                    ops.c_dag('up','0')*ops.c('up','1') + ops.c_dag('up','1')*ops.c('up','0')
                    +ops.c_dag('dn','0')*ops.c('dn','1') + ops.c_dag('dn','1')*ops.c('dn','0')
                    )
                - parms['chemical_potential']*ref_N
                )
        # Run the solver
        ref_ed.diagonalize(ref_H)
        # Compute G(i\omega)
        ref_G_iw = ref_ed.G_iw(ref_gf_struct, parms['beta'], parms['n_iw'])
        return ref_G_iw
                
    ref_G_iw = calc_reference()
    ref=ref_G_iw[ref_spin]

    # Construct the DF2 program
    X = dualfermion.Dpt(beta=parms['beta'],gf_struct=gf_struct,n_iw=parms['n_iw'],n_iw2=parms["measure_G2_n_fermionic"],n_iW =parms["measure_G2_n_bosonic"],N_x=parms['N_x'],N_y=parms['N_y'],N_z=parms['N_z'])

    def Hk_f(k):
        return -2 *parms['t1'] * (np.cos(k[0]))*np.eye(1)

    for spin,_ in X.Hk:
        for k in X.Hk.mesh:
            X.Hk[spin][k] = Hk_f(k.value)

    g2_blocks = set([("up", "up"), ("up", "dn"), ("dn", "up"),("dn","dn")])
    index_converter = {(sn, o) : ("loc", int(o), "down" if sn == "dn" else "up")
                    for sn, o in product(spin_names, orb_names)}
    
    
    # 1 Bath degree of freedom
    # Level of the bath sites
    epsilon = [-parms['chemical_potential_bare'],]
    index_converter.update({("B%i_%s" % (k, sn), 0) : ("bath" + str(k), 0, "down" if sn == "dn" else "up")
                            for k, sn in product(range(len(epsilon)), spin_names)})

    # Make PomerolED solver object
    ed = PomerolED(index_converter, verbose = True)
    N = sum(ops.n(sn, o) for sn, o in product(spin_names, orb_names))
    H_loc = (
                parms["U"] * (ops.n('up','0') * ops.n('dn','0') )
                - parms['chemical_potential']*N
    )
    
    
    # Bath Hamiltonian: levels
    H_bath = sum(eps*ops.n("B%i_%s" % (k, sn), 0)
                for sn, (k, eps) in product(spin_names, enumerate(epsilon)))

    # Hybridization Hamiltonian 
    # Bath-impurity hybridization
    V = [-2 *bath_prefactor* t1,]
    H_hyb = ops.Operator()
    for k, v in enumerate(V):
        H_hyb += sum(        v   * ops.c_dag("B%i_%s" % (k, sn), 0) * ops.c(sn, '0') +
                    np.conj(v)  * ops.c_dag(sn, '0') * ops.c("B%i_%s" % (k, sn), 0)
                    for sn in spin_names)

    for name,g0 in X.Delta:
        X.Delta[name] << gf.inverse( gf.iOmega_n+parms['chemical_potential_bare'])*bath_prefactor**2*4*t1**2 

                
    # Obtain bath sites from Delta and create H_ED
    H_ED = H_loc+H_bath +H_hyb
            
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
    gimp = G_iw[ref_spin]


    # Run the dual perturbation theory
    X.gimp << G_iw # Load G from impurity solver
    dpt_parms = {key:parms[key] for key in parms if key in dptkeys}
    X.run(**dpt_parms)
        
    if mpi.is_master_node():        
        arch = ar.HDFArchive(filename,'w')
        arch["parms"] = parms
        arch["gimp"] = gimp
        arch["ref"] = ref


def plot(filename):
    
    if mpi.is_master_node():        
        arch = ar.HDFArchive(filename,'r')
        parms = arch["parms"]
        sigma_k = HDFArchive("sigma_k.h5",'r')['sigma_k'][ref_spin]
        G_k = HDFArchive("G_k.h5",'r')['G_k'][ref_spin]


        inv = np.linalg.inv
        t1 = parms['t1']
        gimp=arch["gimp"]
        ref=arch["ref"]

        # Calculate lattice self-energy
        for ki,k in enumerate(sigma_k.mesh[1]):
            selfenergylist = []
            ref_list = []
            simp_list = []

            kx = k.value[0]
            tk= -2*t1*(np.cos(kx))
            def factor(ii):
                return np.exp(1j*kx*ii)


            for wn in sigma_k.mesh[0]:
                wval = wn.value.imag
                if wval<0: continue
                S = sigma_k[(wn,k)]+parms['chemical_potential']
                
                invG0k  = wn.value+parms['chemical_potential'] -tk
                #S = invG0k-1./G_k[(wn,k)]
                
                selfenergylist.append( [wval,]+list(S.flatten()) )
                
                # Reference data, convert from 2x2 matrix to complex number with correct k
                ref_G = np.sum(factor(ii)* ref[(wn)][0,ii] for ii in range(2))
                ref_S = invG0k - 1./ref_G
                ref_list.append( [wval,]+list(ref_S.flatten()) )
                #ref_list.append( [wval,]+list(ref_G.flatten()) )
                
                # calculate impurity self-energy as well, from runfile['G_iw']
                g = gimp[(wn)]
                
                # This was done using TODO check factors Delta
                invg0= wn.value+parms['chemical_potential']-bath_prefactor**2*4*t1**2/(wn.value+parms['chemical_potential_bare'])
                Simp = invg0 - 1./g
                simp_list.append( [wval,]+list(Simp.flatten()) )
                
            s = np.array(selfenergylist)
            r = np.array(ref_list)
            s2= np.array(simp_list)
            for column in range(1,2):
                assert np.max(np.abs(s[:5,column]-r[:5,column])) < 1e-4

if __name__ == "__main__":
    for t1 in [0.02,]:
        filename = 'arch.h5'
        run_test(t1,filename)
        plot(filename)

