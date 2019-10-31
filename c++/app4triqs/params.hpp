/*******************************************************************************
 *
 * app4triqs: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
 *
 * app4triqs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * app4triqs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * app4triqs. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "./types.hpp"

namespace app4triqs {

  /// The parameters for the solver construction
  struct constr_params_t {

    /// Number of tau points
    int n_tau = 5001;

    /// Number of Matsubara frequencies
    int n_iw = 500;

    /// Inverse temperature
    double beta;

    /// Block structure of the gf
    gf_struct_t gf_struct;

    /// Number of block indeces for the Green function
    int n_blocks() const { return gf_struct.size(); }

    /// Names of block indeces for the Green function
    auto block_names() const {
      std::vector<std::string> v;
      for (auto const &bl : gf_struct) v.push_back(bl.first);
      return v;
    }

    /// Write constr_params_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp);

    /// Read constr_params_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t &cp);
  };

  /// The parameters for the solve function
  struct solve_params_t {

    // ----------- System Specific -----------

    /// Interaction Hamiltonian
    many_body_operator h_int;

    // ----------- Solver Specific -----------

    /// Maximum running time in seconds (-1 : no limit)
    int max_time = -1;

    /// Verbosity
    int verbosity = mpi::communicator().rank() == 0 ? 3 : 0;

    /// Perform post processing
    bool post_process = true;

    /// Write constr_params_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp);

    /// Read constr_params_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp);
  };

  /// A struct combining both constr_params_t and solve_params_t
  struct params_t : constr_params_t, solve_params_t {
    params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
       : constr_params_t(constr_params_), solve_params_t(solve_params_) {}
  };

} // namespace app4triqs
