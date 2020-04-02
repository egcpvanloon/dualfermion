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

#include <optional>

#include "types.hpp"

namespace triqs_dualfermion {

  /// Containers for measurements
  struct container_set_t {

    // -- Single particle Green's functions

    // -- Two-particle Green's functions

    /// Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` (three Fermionic frequencies)
    std::optional<G2_iw_t> G2_iw;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set_t const &c);

    /// Function that reads all containers to hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set_t &c);

  }; // struct container_set_t

} // namespace triqs_dualfermion
