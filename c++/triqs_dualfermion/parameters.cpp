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

#include "./parameters.hpp"

#include <triqs/utility/itertools.hpp>
using triqs::utility::enumerate;

namespace triqs_dualfermion {

  // -- pair<string, string>
  
  inline void h5_write(triqs::h5::group h5group, std::string name, std::pair<std::string, std::string> const &pair) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.create_group(name);
    h5_write(grp, "0", std::string(pair.first));
    h5_write(grp, "1", std::string(pair.second));
  }

  inline void h5_read(triqs::h5::group h5group, std::string name, std::pair<std::string, std::string> &pair) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.open_group(name);
    assert(grp.get_all_subgroup_names().size() == 2);
    h5_read(grp, "0", pair.first);
    h5_read(grp, "1", pair.second);
  }

  // -- set<pair<string, string>>
  
  inline void h5_write(triqs::h5::group h5group, std::string name, std::set<std::pair<std::string, std::string>> const &pair_set) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.create_group(name);
    for( auto [idx, pair] : enumerate(pair_set) ) {
      h5_write(grp, std::to_string(idx), pair);
    }
  }

  inline void h5_read(triqs::h5::group h5group, std::string name, std::set<std::pair<std::string, std::string>> &pair_set) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.open_group(name);
    for( auto sgrp_name : grp.get_all_subgroup_names() ) {
      std::pair<std::string, std::string> pair;
      h5_read(grp, sgrp_name, pair);
      pair_set.insert(pair);
    }
  }
  
  void h5_write(triqs::h5::group h5group, std::string name, constr_parameters_t const &cp) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.create_group(name);
    h5_write(grp, "beta", cp.beta);
    h5_write(grp, "gf_struct", cp.gf_struct);
    h5_write(grp, "n_iw", cp.n_iw);
    h5_write(grp, "n_iW", cp.n_iW);
    h5_write(grp, "n_iw2", cp.n_iw2);
    h5_write(grp, "N_x", cp.N_x);
    h5_write(grp, "N_y", cp.N_y);
    h5_write(grp, "N_z", cp.N_z);    
  }

  void h5_read(triqs::h5::group h5group, std::string name, constr_parameters_t &cp) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.open_group(name);
    h5_read(grp, "beta", cp.beta);
    h5_read(grp, "gf_struct", cp.gf_struct);
    h5_read(grp, "n_iw", cp.n_iw);
    h5_read(grp, "n_iW", cp.n_iW);
    h5_read(grp, "n_iw2", cp.n_iw2);
    h5_read(grp, "N_x", cp.N_x);
    h5_read(grp, "N_y", cp.N_y);
    h5_read(grp, "N_z", cp.N_z);
  }

  void h5_write(triqs::h5::group h5group, std::string name, run_parameters_t const &sp) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.create_group(name);

    h5_write(grp, "verbosity", sp.verbosity);
        
    h5_write(grp, "sigmad_subset", sp.sigmad_subset);
    h5_write(grp, "calculate_sigma", sp.calculate_sigma);
    h5_write(grp, "calculate_sigma1", sp.calculate_sigma1);
    h5_write(grp, "calculate_sigma2", sp.calculate_sigma2);
    
  }

  void h5_read(triqs::h5::group h5group, std::string name, run_parameters_t &sp) {
    triqs::h5::group grp = name.empty() ? h5group : h5group.open_group(name);
    
    h5_read(grp, "verbosity", sp.verbosity);    
    
    h5_read(grp, "sigmad_subset", sp.sigmad_subset);
    h5_read(grp, "calculate_sigma", sp.calculate_sigma);
    h5_read(grp, "calculate_sigma1", sp.calculate_sigma1);
    h5_read(grp, "calculate_sigma2", sp.calculate_sigma2);
  }
  
} // namespace triqs_dualfermion
