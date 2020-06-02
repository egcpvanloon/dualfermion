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
#include "./dpt_core.hpp"
#include "./util.hpp"

#include <triqs/utility/callbacks.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/gfs.hpp>
#include <fstream>
#include <triqs/utility/variant.hpp>

#include <triqs/lattice/brillouin_zone.hpp>
#include <triqs/lattice/bravais_lattice.hpp>
#include <triqs/lattice/gf_mesh_brillouin_zone.hpp>
using triqs::lattice::brillouin_zone;
using namespace triqs::gfs;

using triqs::clef::placeholder;

#include <triqs/arrays.hpp>
using triqs::arrays::array;

mpi::communicator world;

namespace triqs_dualfermion {

  dpt_core::dpt_core(constr_parameters_t const &p)
     : constr_parameters(p), beta(p.beta), gf_struct(p.gf_struct),n_iw(p.n_iw),n_iW(p.n_iW),n_iw2(p.n_iw2),N_x(p.N_x),N_y(p.N_y),N_z(p.N_z) {
         
    // Allocate single particle greens functions
    _gimp   = block_gf<imfreq>({beta, Fermion, n_iw}, gf_struct);
    _Delta  = block_gf<imfreq>({beta, Fermion, n_iw}, gf_struct);

    matrix<int> periodization_matrix(3,3) ;
    periodization_matrix() = 0;
    periodization_matrix(0,0) = N_x;
    periodization_matrix(1,1) = N_y;
    periodization_matrix(2,2) = N_z;
    
    const bravais_lattice lat ;
    const brillouin_zone bz(lat) ;
    gf_mesh<brillouin_zone> kmesh( bz, periodization_matrix);

    _Hk     = block_gf<brillouin_zone>({kmesh},gf_struct);
    
    
     gf_mesh<imfreq> mesh_f{beta, Fermion, n_iw2};
     gf_mesh<imfreq> mesh_b{beta, Boson, n_iW};
     gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_bff{mesh_b, mesh_f, mesh_f};
    _G2_iw  = make_block2_gf(mesh_bff, gf_struct, p.G2_block_order);
    vertex  = make_block2_gf(mesh_bff, gf_struct, p.G2_block_order);
  }

  /// -------------------------------------------------------------------------------------------

  void dpt_core::run(run_parameters_t const &run_parameters_) {

    //run_parameters = run_parameters_;
    run_parameters_t params(run_parameters_);
    gf_struct_t sigmad_subset(params.sigmad_subset);

    /* Set-up boolean variables that determine which output is written
     * TODO: make these individual flags that are accessible in python
     * still check for calculate_sigma and world.rank though. The latter avoid simultaneous writing to files
     */                   
    bool output_stdout = params.verbosity >= 2 && (world.rank() == 0);
    
    bool output_sigmak = params.verbosity >= 1 && params.calculate_sigma && (world.rank() == 0);
    bool output_gk = params.verbosity >= 1  && params.calculate_sigma && (world.rank() == 0);
    bool output_gdk = params.verbosity >= 3 && params.calculate_sigma && (world.rank() == 0);
    bool output_gdr = params.verbosity >= 3 && params.calculate_sigma && (world.rank() == 0);
    bool output_sigmadk = params.verbosity >= 3 && params.calculate_sigma && (world.rank() == 0);
    bool output_sigmadr = params.verbosity >= 3 && params.calculate_sigma && (world.rank() == 0);
    bool output_vertex = params.verbosity >= 5 && params.calculate_sigma && (world.rank() == 0);
    
    
    if (output_stdout)
      std::cout << "\n"
                   " ------------------------------\n"
                   " - second-order dual fermion  -\n"
                   " ------------------------------\n\n";
    
                   
    auto kmesh = _Hk[0].mesh(); // Note that Hk also has a spin index, Hk[0] takes Hk for a single spin. 
    
    // Placeholders
    placeholder<0> wn_; // Frequency -- Fermionic
    //placeholder<1> wn1_; // Frequency -- Fermionic 
    //placeholder<2> wn2_; // Frequency -- Fermionic
    //placeholder<3> Wn_; // Frequency -- Bosonic
    
    
    // Create initial guess for Delta
    // Assumes that 1/iw has been loaded into gimp
    // And iw into Delta
    if(params.delta_initial){
        auto gloc  = G_iw_t{{beta,Fermion,n_iw},gf_struct} ;  
    
        // First set to zero
        gloc() = 0. ;
    
        for (auto const &b : range(gloc.size())) {
            for (auto const &iw : gloc[b].mesh()){
                for (auto const &k : kmesh ){
                    gloc[b][iw]  += _gimp[b](iw)*1./(1. -_Hk[b](k) *_gimp[b](iw) );
                }
            }
        }
        gloc() /= double(kmesh.size());
        _Delta(wn_) << _Delta(wn_)-1/gloc(wn_)  ;
        
        if (output_stdout) std::cout<< "Initial Delta has been assigned" << std::endl;
        return;
    }
    


    // Create the container for the bare dual Green's function
    auto gd0 = G_iw_k_t{{ {beta,Fermion,n_iw} ,kmesh},gf_struct} ;  
    
    if (output_stdout) std::cout << "Initializing gd0" <<std::endl;    
    // Load the values of the bare dual Green's function                   
    for (auto const &b : range(gd0.size())) {
      for (const auto &[iw,k] : gd0[b].mesh()){          
          gd0[b][iw,k] = 1./ ( _gimp[b][iw] + _gimp[b][iw]*(_Delta[b][iw] - _Hk[b](k)   )*_gimp[b][iw]) - 1./_gimp[b][iw] ;          
      }
    }
    if (output_gdk) h5_write(h5::file("gd_k.h5",'w'),"gd_k",gd0);
    
    
    /* Initialize empty Sigma_dual in momentum space */
    auto sigmad = G_iw_k_t{{ {beta,Fermion,n_iw} ,kmesh},gf_struct} ;

    
    if(params.calculate_sigma){
    //int s1 =0;    
    //int s2 =0;

    if (output_stdout) std::cout << "Calculating Sigma dual" << std::endl;
    
    /* FT gd0 to real space  */
    //gf_mesh<cyclic_lattice> rmesh( lat, periodization_matrix);
    auto rmesh = make_adjoint_mesh(kmesh) ;
    auto R0 = rmesh[{0,0,0}];
    
    //auto gd_real = G_iw_r_t{{ {beta,Fermion,n_iw} ,rmesh},gf_struct} ;
    //gd_real() = fourier<1>(gd0);
    
    auto gd_real = make_gf_from_fourier<1>(gd0);    

    gf_mesh<imfreq> mesh_f{beta, Fermion, n_iw};    
    gf_mesh<cartesian_product<imfreq, cyclic_lattice>> mesh_fr{mesh_f, rmesh};    
    sigmad_real = block_gf<cartesian_product<imfreq,cyclic_lattice>>({mesh_fr},gf_struct);
    //auto sigmad_real = G_iw_r_t{{ {beta,Fermion,n_iw} ,rmesh},gf_struct} ;
    sigmad_real() = 0;
            
    if (output_stdout) std::cout << "Calculating Sigma dual: assign vertex" << std::endl;

    // Calculate the self-energy
    // Note, need to determine the vertex from G2 and gimp    
    
    
    // TODO:
    // Loops over products could be quite slow (at present), see:
    // https://github.com/TRIQS/triqs/issues/658
    // Code can be made more performant by replacing the loops by individual loops, ugly
    // Hopefully, in the future this goes away
    // In that case, consider rewriting in CLEF notation becomes sufficiently performant
    
    // Calculate vertex first, since it will be used several times. 
    // Be careful about the vertex band labels here so that they can be used safely later on.    
    // Slightly tricky: different Matsubara grids for _gimp and _G2_iw, have to use _gimp[0](n1) instead of _gimp[0][n1]
    
    //TODO: Timing/memory: is it necessary to precalculate the vertex? 
    //      Would it be better to expose this as a separate method in python?
    vertex() = 0.;
    for (auto const &s2 : range(_gimp.size())) {
      for (auto const &s1 : range(_gimp.size())) {
        for (const auto &[iw, n1, n2] : _G2_iw(s1,s2).mesh()){
          for (const auto [Ai, Aj, Ak, Al] : _G2_iw(s1,s2).target_indices()){

            vertex(s1,s2)[iw, n1, n2](Ai, Aj, Ak, Al) 
            = 
            _G2_iw(s1,s2)[iw, n1, n2](Ai, Aj, Ak, Al)
            +beta*kronecker(s1,s2)*kronecker(n1,n2)*_gimp[s1](n1)(Al,Ai)*_gimp[s2](n1+iw)(Aj,Ak) 
            - beta* _gimp[s1](n1)(Aj,Ai)*_gimp[s2](n2)(Al,Ak)* kronecker(iw) ;
                
            
          }
        }
      }
    }   
    if (output_vertex) h5_write(h5::file("vertex.h5",'w'),"vertex",vertex);
    
    if(params.calculate_sigma1){
    if (output_stdout) std::cout << "Calculating Sigma dual: Calculate first order diagram" << std::endl;
    // First-order diagram    
    // Outer loop over s2, the spin of Sigma_dual. This takes values in sigmad_subset =< gf_struct
    // Inner loop over s1, the spin of the internal Green's function
    // R0 is (0,0,0) in real space
    // i,j,k,l are the band indices
    // iw,n1,n2 are the frequency indices
    int s2=0;
    for (auto const &bl : sigmad_subset) {
      int s1=0;  
      for (auto const &bl2 : gf_struct) {
        for (const auto &[iw, n1, n2] : mpi::chunk(vertex(s1,s2).mesh())){
          if(not kronecker(iw)) continue; // Only w=0 contributes to the second-order diagram 
          for (const auto [i, j, k, l] : vertex(s1,s2).target_indices()){
            sigmad_real[s2][n2,R0](l,k) += -vertex(s1,s2)[iw, n1, n2](i, j, k, l) * gd_real[s1](n1,R0)(i,j)/beta ; 
            
          }
        }
        s1++;
      }
      s2++;
    }
    }//sigma1
    
    if(params.calculate_sigma2){
    if (output_stdout) std::cout << "Calculating Sigma dual: Second-order diagram" << std::endl;
    // Second-order diagram
    // Outer loop over s2, the spin of Sigma_dual. This takes values in sigmad_subset =< gf_struct
    // Inner loop over s1, the spin of two of the internal Green's functions
    // R is the spatial argument of Sigma in real space, two Green's functions have R and 1 has -R
    // i,j,k,l are the band indices
    // iw,n1,n2 are the frequency indices
    // There are two spin configurations possible, these are handled separately here. TODO: rewrite order of loops to make this more efficient
    // TODO: Check explicitly that all frequencies are in range in a symmetric way
    auto orbital_indices = vertex(0,0).target_indices();
    for (const auto &[iw, n1, n2] : mpi::chunk(vertex(0,0).mesh(),world)){
      //TESTING: if(not kronecker(iw)) continue; Only zero frequency  
      for (const auto R : rmesh){  
        for (const auto [Ai, Aj, Ak, Al] : orbital_indices){
          for (const auto [Bi, Bj, Bk, Bl] : orbital_indices){                                
            int s2=0;
            //TODO: Pull s2 to the outer level?
            //      Introducing an option to let the loop go over a subset 
            //      of indices, ''upfolding'' to all indices can then be done later
            //      N.b.: Current code assumes that the subset is continuous and starts at 0, since it uses the integer s1,s2 to count blocks. Instead, rewrite to use bl
            for (auto const &bl : sigmad_subset) {
              int s1=0;
              for (auto const &bl2 : gf_struct) {
                // First spin configuration s1 != s2 
                sigmad_real[s2][n2,R](Al,Bk) += 
                -0.5
                *(1-kronecker(s1,s2))
                * vertex(s1,s2)(n1-n2,n2+iw, n2)(Ak, Aj, Ai, Al) 
                * vertex(s1,s2)(n2-n1, n1+iw, n1)(Bi, Bl, Bk, Bj) 
                * gd_real[s1](n2+iw,R)(Ak,Bl)
                * gd_real[s1](n1+iw,-R)(Bi,Aj)
                * gd_real[s2](n1,R)(Ai,Bj)
                /(beta*beta) ;                
                // Second spin configuration: all s1  
                sigmad_real[s2][n2,R](Al,Bk) += 
                -0.5
                * vertex(s1,s2)(iw, n1, n2)(Ai, Aj, Ak, Al) 
                * vertex(s1,s2)(-iw, n1+iw, n2+iw)(Bi, Bj, Bk, Bl) 
                * gd_real[s1](n1,R)(Ai,Bj)
                * gd_real[s1](n1+iw,-R)(Bi,Aj)
                * gd_real[s2](n2+iw,R)(Ak,Bl)
                /(beta*beta) ;
                s1++;                                
              }
              s2++;
            }
          }
        }
      }
    }
    }//sigma2    
    sigmad_real =  mpi::reduce(sigmad_real, world);
    
    //TODO:
    // Perform the ''upfolding'' here? Or in python
    
    if (output_stdout) std::cout << "Calculating Sigma dual: FT" << std::endl;
    sigmad = make_gf_from_fourier<1>(sigmad_real);
    if (output_gdr) h5_write(h5::file("gd_real.h5",'w'),"gd_real",gd_real);
    if (output_sigmadr) h5_write(h5::file("sigmad_real.h5",'w'),"sigmad_real",sigmad_real);
    if (output_sigmadk) h5_write(h5::file("sigmad_k.h5",'w'),"sigmad_k",sigmad);
    }
    else{
      sigmad() = 0;            
    }
    
    /* Calculate G lattice and Gdlattice */
    
    if (output_stdout) std::cout << "Calculating Glat" << std::endl;
    // Create the container for the lattice Green's function
    auto glat = G_iw_k_t{{ {beta,Fermion,n_iw} ,kmesh},gf_struct} ;  
    for (auto const &b : range(glat.size())) {
        for (const auto &[iw,k] : glat[b].mesh()){
            glat[b][iw,k] = 1./( 1./(_gimp[b][iw]+sigmad[b][iw,k]) + _Delta[b][iw] - _Hk[b](k)  );
            gd0[b][iw,k]  = gd0[b][iw,k]/( 1.-sigmad[b][iw,k]*gd0[b][iw,k]);
        }
    }
    // Write the lattice Green's function
    if (output_gk) h5_write(h5::file("G_k.h5",'w'),"G_k",glat);
    
    // Also calculate sigma_lat here from Dyson's equation, since that is what we are frequently interested in
    // TODO: Some of the tests could be written purely in terms of sigma?
    if (output_sigmak){
        auto sigma_k = G_iw_k_t{{ {beta,Fermion,n_iw} ,kmesh},gf_struct} ;  
        for (auto const &b : range(sigma_k.size())) {
            for (const auto &[iw,k] : sigma_k[b].mesh()){
                sigma_k[b][iw,k] = iw - _Hk[b](k)  -1./glat[b][iw,k];
            }
        }
        h5_write(h5::file("sigma_k.h5",'w'),"sigma_k",sigma_k);
    }// output_sigmak

    
    /* Calculate local part of Gd and Gloc */
    auto gloc  = G_iw_t{{beta,Fermion,n_iw},gf_struct} ;  
    auto gdloc = G_iw_t{{beta,Fermion,n_iw},gf_struct} ;
    
    // First set to zero
    gloc() = 0. ;
    gdloc()= 0. ;
    
    // Get local part
    for (auto const &b : range(glat.size())) {
        for (const auto &[iw,k] : glat[b].mesh()){            
            gloc[b][iw]  += glat[b][iw,k];
            gdloc[b][iw] += gd0[b][iw,k];
        }
    }
    gloc() /= double(kmesh.size());
    gdloc() /= double(kmesh.size());
    
    
    /* Update formula        
       Delta_new = Delta_old + ksi * gimp^-1 Gd_loc Gloc^-1
       Now in terms of Gdd = 1/gimp * Gd * 1/gimp
     */
    if (output_stdout) std::cout << "Evaluate new Delta" << std::endl;
    _Delta(wn_) << _Delta(wn_) + gdloc(wn_) * _gimp(wn_)/gloc(wn_);
    
    
    if (output_stdout) std::cout << "dualfermion finished\n" ;
    
    return ;
    
  }
} // namespace triqs_dualfermion
