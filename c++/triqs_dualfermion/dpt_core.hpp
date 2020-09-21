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
#include <triqs/utility/callbacks.hpp>

#include "types.hpp"
#include "container_set.hpp"
#include "parameters.hpp"

#include <triqs/lattice/brillouin_zone.hpp>
#include <triqs/lattice/bravais_lattice.hpp>
#include <triqs/lattice/gf_mesh_brillouin_zone.hpp>


namespace triqs_dualfermion {

  /// Core class of the df2 program
  class dpt_core : public container_set_t {

    double beta;           // inverse temperature
    gf_struct_t gf_struct; // Block structure of the Green function
    int n_iw;
    int n_iw2;
    int n_iW;
    
    int _status;                     // Status of the run upon exit: 0 for clean termination, > 0 otherwise.

    // Return reference to container_set
    //container_set_t &result_set() { return static_cast<container_set_t &>(*this); }
    //container_set_t const &result_set() const { return static_cast<container_set_t const &>(*this); }

    // Single-particle Green's function containers
    G_iw_t _gimp;      // Matsubara Green's function of the impurity model
    G_iw_t _Delta;     // Matsubara hybridization function of the impurity model
    
    G_k_t  _Hk ;       // Stores the dispersion
    
    // Two-particle Green's function container
    G2_iw_t _G2_iw;    // Two-particle Green's function of the impurity model (1 bosonic, 2 fermionic Matsubara frequencies)
                       // N.B.: Loaded directly from TRIQS/CTHYB, not yet connected and amputated
    G2_iw_t vertex;    // Connected, amputated vertex. Data is stored in _gimp and _G2_iw

    // Struct containing the parameters relevant for the dpt construction
    constr_parameters_t constr_parameters;

    // Struct containing the parameters of the last call to the run method
    run_parameters_t run_parameters;

    /// Order of block indices in the definition of G^2.
    block_order G2_block_order = block_order::AABB;    
    
    public:

    /**
     * Construct the dualfermion program
     *
     * @param p Set of parameters specific to the dualfermion program
     */
    CPP2PY_ARG_AS_DICT
    dpt_core(constr_parameters_t const &p);

    // Delete assignement operator because of const members
    dpt_core(dpt_core const &p) = default;
    dpt_core(dpt_core &&p)      = default;
    dpt_core &operator=(dpt_core const &p) = delete;
    dpt_core &operator=(dpt_core &&p) = default;

    /**
     * Solve method that performs dual perturbation theory
     *
     * @param p Set of parameters for the dual perturbation theory
     */
    CPP2PY_ARG_AS_DICT
    void run(run_parameters_t const &p);

    /// Set of parameters used in the construction of the ``dpt_core`` class.
    constr_parameters_t last_constr_parameters() const { return constr_parameters; }

    /// Set of parameters used in the last call to ``run()``.
    run_parameters_t last_run_parameters() const { return run_parameters; }

    /// :math:`\Delta(i\omega)` in imaginary frequencies.
    block_gf_view<imfreq> Delta() { return _Delta; }

    /// :math:`g(i\omega)` in imaginary frequencies.
    block_gf_view<imfreq> gimp() { return _gimp; }
    
    block2_gf_view<imfreq_cube_mesh_t, tensor_valued<4> > G2_iw() { return _G2_iw; }
    
    /// :math:`H(k)` on the Brillouin Zone
    block_gf_view<brillouin_zone> Hk() { return _Hk; }

    /// Status of the ``run()`` on exit.
    int status() const { return _status; }

    static std::string hdf5_format() { return "dualfermion_DptCore"; }

    // Function that writes the dpt_core to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, dpt_core const &s) {
      h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
      h5_write_attribute(grp, "Format", dpt_core::hdf5_format());
      //h5_write(grp, "", s.result_set());
      h5_write(grp, "constr_parameters", s.constr_parameters);
      h5_write(grp, "run_parameters", s.run_parameters);
      h5_write(grp, "gimp_iw", s._gimp);
      h5_write(grp, "G2_iw", s._G2_iw);
      h5_write(grp, "Delta_iw", s._Delta);
      h5_write(grp, "Hk", s._Hk);
    }

    // Function that read all containers to hdf5 file
    CPP2PY_IGNORE
    static dpt_core h5_read_construct(h5::group h5group, std::string subgroup_name) {
      h5::group grp   = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
      auto constr_parameters = h5::h5_read<constr_parameters_t>(grp, "constr_parameters");
      auto s                 = dpt_core{constr_parameters};
      h5_read(grp, "run_parameters", s.run_parameters);
      h5_read(grp, "gimp_iw", s._gimp);
      h5_read(grp, "G2_iw", s._G2_iw);
      h5_read(grp, "Delta_iw", s._Delta);
      h5_read(grp, "Hk", s._Hk);
      return s;
    }
  };
} // namespace triqs_dualfermion
