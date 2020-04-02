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

from pytriqs.archive import HDFArchive
from pytriqs.gf import *

import triqs_df2.plothelper as plothelper

import numpy as np

import itertools

from pytriqs.lattice import BravaisLattice, BrillouinZone
from pytriqs.gf import Gf, MeshProduct, MeshBrillouinZone, MeshImFreq

n_k = 32
n_w = 20
t=1
beta = 10.

BL = BravaisLattice([(1, 0, 0), (0, 1, 0)]) # Two unit vectors in R3
BZ = BrillouinZone(BL)

kmesh = MeshBrillouinZone(BZ, n_k=n_k)
wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_w)

g0 = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=[])  # g0(k,omega), scalar valued

def eps(k):
    return -2 * t* (np.cos(k.value[0]) + np.cos(k.value[1]))

# NB : loop is a bit slow in python ...
for k in g0.mesh[1]:
    for w in g0.mesh[0]:
        g0[w,k] = 1/(w - eps(k))

#name = "gd_k"
#G_k = HDFArchive(name+".h5",'r')[name]


spin = 'up'

def test_colormap():
    w0 =[wn for wn in g0.mesh[0] if wn.value.imag>0][0]
    
    c = lambda k : g0[(w0,k)].imag
    x = lambda k : k.value[0]    
    y = lambda k : k.value[1]
    
    plothelper.mesh_plot(g0.mesh[1],x,y,c)

def test_epsilonk():

    def c(k):
        return eps(k)
    def x(k):
        return k.value[0]    
    def y(k):
        return k.value[1]
    
    
    plothelper.mesh_plot(g0.mesh[1],x,y,c)

def test_wedge():

    w0 =[wn for wn in g0.mesh[0] if wn.value.imag>0][0]
    def c(k):
        return g0[(w0,k)].imag
    def x(k):
        return k.value[0]
    
    def y(k):
        return k.value[1]
    
    def check(k):
        return k.value[0] <= np.pi and k.value[1] <=k.value[0]
    
    klist = [k for k in g0.mesh[1] if check(k) ]
    
    plothelper.mesh_plot(klist,x,y,c)

def test_hs():
        
    def on_line(K,A,B,err=10**-6):
        """ Is K on the line between A and B? """
        if np.linalg.norm(K-A)<err: return True
        return abs(np.dot(K-A,B-A)/(np.linalg.norm(K-A)*np.linalg.norm(B-A)) - 1) < err and np.linalg.norm(K-A) <= np.linalg.norm(B-A)
        
    
    def make_section(A,B):
        for k in g0.mesh[1]:
            K = np.array(k.value[0:2])
            if on_line(K,A,B):  yield k
            
    def ord_section(a,b):
        A = np.array(a[0:2])
        B = np.array(b[0:2])
        def key_helper(X):
            return np.linalg.norm(np.array(X.value[0:2])-A)
        return sorted(make_section(A,B),key=key_helper)
        
    def make_path(hs_points):
        sections = [ord_section(hs_points[n],hs_points[(n+1)%len(hs_points)]) for n in range(len(hs_points))]
            
        for k in itertools.chain(*sections) : yield k

    hs_points = [(0,0),(np.pi,0),(np.pi,np.pi)]
    klist = [(n,k) for n,k in enumerate(make_path(hs_points)) ]
    
    def y(k):
        return eps(k[1])
    def x(k):
        return k[0]
    
    pass
    
    plothelper.mesh_plot(klist,x,y)
    
        
if __name__ == '__main__':
    test_colormap()    
    test_epsilonk()    
    test_wedge()    
    
    test_hs()
