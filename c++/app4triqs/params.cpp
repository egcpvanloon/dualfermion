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
#include "./params.hpp"

namespace app4triqs {

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "n_tau", cp.n_tau);
    h5_write(grp, "n_iw", cp.n_iw);
    h5_write(grp, "beta", cp.beta);
    h5_write(grp, "gf_struct", cp.gf_struct);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t &cp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "n_tau", cp.n_tau);
    h5_read(grp, "n_iw", cp.n_iw);
    h5_read(grp, "beta", cp.beta);
    h5_read(grp, "gf_struct", cp.gf_struct);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "h_int", sp.h_int);
    h5_write(grp, "max_time", sp.max_time);
    h5_write(grp, "verbosity", sp.verbosity);
    h5_write(grp, "post_process", sp.post_process);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    auto grp = h5group.open_group(subgroup_name);
    // Take care! Do not read verbosity as they should be different based on mpi rank
    h5_read(grp, "h_int", sp.h_int);
    h5_read(grp, "max_time", sp.max_time);
    h5_read(grp, "post_process", sp.post_process);
  }

} // namespace app4triqs
