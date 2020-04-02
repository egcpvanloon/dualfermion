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

import matplotlib.pyplot as plt


def mesh_plot(klist,x,y,c=None,fig=None,**kwargs):
    """
        Makes a plot
        The general structure is that klist contains mesh points k,
        x,y and c should be functions that map k to a real number

        
        klist can be a list or generator
    """
    
    if fig is None:
        fig=plt.figure()
        

    x_list=[ x(k) for k in klist]
    y_list=[ y(k) for k in klist]

    if c is not None:
        c_list=[ c(k) for k in klist]
    else:
        c_list = None
        
        
    plt.scatter(x_list,y_list,c=c_list,**kwargs)        


# Examples:
# epsilon(k) colormap in 2d
# epsilon(k) 2d cross-section in 3d
# epsilon(k) on Gamma - X - M - Gamma path
#
# G(w0,k) colormap on the irreducible wedge
