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

  /// The collection of all output containers in solver_core
  struct container_set {

    /// Greens function in imaginary time
    g_tau_t G_tau;

    /// Greens function in Matsubara frequencies
    g_iw_t G_iw;

    /// Self-energy in Matsubara frequencies
    g_iw_t Sigma_iw;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set const &c);

    /// Function that reads all containers from hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set &c);
  };

} // namespace app4triqs
