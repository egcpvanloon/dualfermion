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

#define STR(x) #x
#define STRINGIZE(x) STR(x)

#include <triqs/gfs.hpp>
#include <triqs/utility/time_pt.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp> // gf_struct_t

#include <triqs/lattice/brillouin_zone.hpp>
#include <triqs/lattice/bravais_lattice.hpp>
#include <triqs/lattice/gf_mesh_brillouin_zone.hpp>


#include <variant>

namespace triqs_dualfermion {

  using namespace triqs::gfs;
  using namespace triqs::utility;

  using triqs::hilbert_space::gf_struct_t;
  using triqs::utility::time_pt;

  using gf_struct_t  = triqs::hilbert_space::gf_struct_t;

  // One-particle Green's function types
  using G_iw_t           = block_gf<imfreq, matrix_valued>;
  using iw_k_mesh_t      = cartesian_product<imfreq,brillouin_zone>;
  using G_iw_k_t         = block_gf<iw_k_mesh_t, matrix_valued>;
  using iw_r_mesh_t      = cartesian_product<imfreq,cyclic_lattice>;
  using G_iw_r_t         = block_gf<iw_r_mesh_t, matrix_valued>;
  using G_k_t            = block_gf<brillouin_zone, matrix_valued>;

  // Two-particle Green's function types
  using imfreq_cube_mesh_t = cartesian_product<imfreq, imfreq, imfreq>;
  using G2_iw_t            = block2_gf<imfreq_cube_mesh_t, tensor_valued<4>>;

  /// Order of block indices for Block2Gf objects
  enum class block_order { AABB, ABBA };
  
  inline void h5_write(h5::group h5group, std::string name, block_order const &bo) {
    h5_write(h5group, name, static_cast<int>(bo));
  }

  inline void h5_read(h5::group h5group, std::string name, block_order &bo) {
    int idx;
    h5_read(h5group, name, idx);
    bo = static_cast<block_order>(idx);
  }

} // namespace triqs_dualfermion

namespace triqs {
  namespace gfs {

    /// Function template for block2_gf initialization
    template <typename Var_t>
    block2_gf<Var_t, tensor_valued<4>> make_block2_gf(gf_mesh<Var_t> const &m, triqs::hilbert_space::gf_struct_t const &gf_struct,
                                                      triqs_dualfermion::block_order order = triqs_dualfermion::block_order::AABB) {

      std::vector<std::vector<gf<Var_t, tensor_valued<4>>>> gf_vecvec;
      std::vector<std::string> block_names;

      for (auto const &bl1 : gf_struct) {
        auto &bname  = bl1.first;
        int bl1_size = bl1.second.size();
        block_names.push_back(bname);
        std::vector<std::string> indices1;
        for (auto const &var : bl1.second) visit([&indices1](auto &&arg) { indices1.push_back(std::to_string(arg)); }, var);

        std::vector<gf<Var_t, tensor_valued<4>>> gf_vec;
        for (auto const &bl2 : gf_struct) {
          int bl2_size = bl2.second.size();
          std::vector<std::string> indices2;
          for (auto const &var : bl2.second) visit([&indices2](auto &&arg) { indices2.push_back(std::to_string(arg)); }, var);
          auto I = std::vector<std::vector<std::string>>{indices1, indices1, indices2, indices2};
          switch (order) {
            case triqs_dualfermion::block_order::AABB: gf_vec.emplace_back(m, make_shape(bl1_size, bl1_size, bl2_size, bl2_size), I); break;
            case triqs_dualfermion::block_order::ABBA: gf_vec.emplace_back(m, make_shape(bl1_size, bl2_size, bl2_size, bl1_size), I); break;
          }
        }
        gf_vecvec.emplace_back(std::move(gf_vec));
      }
      return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
    }

  } // namespace gfs
} // namespace triqs
