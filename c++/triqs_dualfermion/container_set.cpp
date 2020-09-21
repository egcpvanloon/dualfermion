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

#include "./container_set.hpp"

namespace triqs_dualfermion {

  /// Function that writes all containers to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, container_set_t const &c) {
    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
  }

  /// Function that reads all containers to hdf5 file
  void h5_read(h5::group h5group, std::string subgroup_name, container_set_t &c) {
    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
  }

} // namespace triqs_dualfermion
