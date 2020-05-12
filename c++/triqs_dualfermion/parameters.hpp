/*******************************************************************************
 *
 * dualfermion: second-order dual fermion perturbation theory
 * 
 * Copyright (C) 2020, E.G.C.P. van Loon
 * 
 * Based on TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * dualfermion is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * dualfermion is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * dualfermion. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "./types.hpp"
#include <boost/mpi.hpp>
#include <triqs/operators/many_body_operator.hpp>

namespace triqs_dualfermion {

  using namespace triqs::operators;
  using indices_map_t = std::map<triqs::operators::indices_t, triqs::operators::indices_t>;

  // All the arguments of the dpt_core constructor
  struct constr_parameters_t {
    
    /// Inverse temperature
    double beta;

    ///block structure of the gf
    gf_struct_t gf_struct;

    /// Number of Matsubara frequencies for gf<imfreq, matrix_valued>
    int n_iw = 128;
    int n_iW = 31;
    int n_iw2 = 32;
    
    int N_x = 1;
    int N_y = 1;
    int N_z = 1;
    
    /// Write constr_parameters_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_parameters_t const &sp);

    /// Read constr_parameters_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_parameters_t &sp);
    
    /// Order of block indices in the definition of G^2.
    block_order G2_block_order = block_order::AABB;    
    
  };
  
  // All the arguments of the run function
  struct run_parameters_t {

    /// Verbosity level
    /// default: 3 on MPI rank 0, 0 otherwise.
    int verbosity = ((boost::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes
    
    
    gf_struct_t sigmad_subset;
    
    //double t1 = 1.;
    
    //double ksi_delta = 1.;
    
    bool delta_initial = false;
    
    bool calculate_sigma = false;

    bool calculate_sigma1 = true;
    bool calculate_sigma2 = true;

    
    run_parameters_t() {}

    /// Write run_parameters_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, run_parameters_t const &sp);

    /// Read run_parameters_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, run_parameters_t &sp);

  };

  /// A struct combining both constr_params_t and run_params_t
  struct params_t : constr_parameters_t, run_parameters_t {
    params_t(constr_parameters_t constr_parameters_, run_parameters_t run_parameters_)
      : constr_parameters_t(constr_parameters_), run_parameters_t(run_parameters_) {}
  };
}
