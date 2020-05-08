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


namespace triqs_dualfermion {

  dpt_core::dpt_core(constr_parameters_t const &p)
     : constr_parameters(p), beta(p.beta), gf_struct(p.gf_struct), n_iw(p.n_iw),n_iW(p.n_iW),n_iw2(p.n_iw2),N_x(p.N_x),N_y(p.N_y),N_z(p.N_z) {
         
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

    if (params.verbosity >= 2)
      std::cout << "\n"
                   " ------------------------------\n"
                   " - second-order dual fermion  -\n"
                   " ------------------------------\n\n";
    
    auto kmesh = _Hk[0].mesh(); // Note that Hk also has a spin index, Hk[0] takes Hk for a single spin. 
    
    // Placeholders
    placeholder<0> wn_; // Frequency -- Fermionic
    placeholder<1> wn1_; // Frequency -- Fermionic 
    placeholder<2> wn2_; // Frequency -- Fermionic
    placeholder<3> Wn_; // Frequency -- Bosonic
    
    
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
        
        if (params.verbosity >= 2) std::cout<< "Initial Delta has been assigned" << std::endl;
        return;
    }
    


    // Create the container for the bare dual Green's function
    auto gd0 = G_iw_k_t{{ {beta,Fermion,n_iw} ,kmesh},gf_struct} ;  
    
    if (params.verbosity >= 2) std::cout << "Initializing gd0" <<std::endl;    
    // Load the values of the bare dual Green's function                   
    for (auto const &b : range(gd0.size())) {
      for (const auto &[iw,k] : gd0[b].mesh()){          
          gd0[b][iw,k] = 1./ ( _gimp[b][iw] + _gimp[b][iw]*(_Delta[b][iw] - _Hk[b](k)   )*_gimp[b][iw]) - 1./_gimp[b][iw] ;          
      }
    }
    if (params.verbosity >= 1) h5_write(h5::file("gd_k.h5",'w'),"gd_k",gd0);
    
    
    /* Initialize empty Sigma_dual in momentum space */
    auto sigmad = G_iw_k_t{{ {beta,Fermion,n_iw} ,kmesh},gf_struct} ;

    
    if(params.calculate_sigma){
    int s1 =0;    
    int s2 =0;

    if (params.verbosity >= 2) std::cout << "Calculating Sigma dual" << std::endl;
    
    /* FT gd0 to real space  */
    //gf_mesh<cyclic_lattice> rmesh( lat, periodization_matrix);
    auto rmesh = make_adjoint_mesh(kmesh) ;
    auto R0 = rmesh[{0,0,0}];
    
    //auto gd_real = G_iw_r_t{{ {beta,Fermion,n_iw} ,rmesh},gf_struct} ;
    //gd_real() = fourier<1>(gd0);
    
    auto gd_real = make_gf_from_fourier<1>(gd0);    
    
    auto sigmad_real = G_iw_r_t{{ {beta,Fermion,n_iw} ,rmesh},gf_struct} ;
    sigmad_real() = 0;
            
    if (params.verbosity >= 2) std::cout << "Calculating Sigma dual: assign vertex" << std::endl;

    // Calculate the self-energy
    // Note, need to determine the vertex from G2 and gimp    
    
    
    // TODO:
    // Loops over products are quite slow (at present), see:
    // https://github.com/TRIQS/triqs/issues/658
    // Code can be made more performant by replacing the loops by individual loops, ugly
    // Hopefully, in the future this goes away
    // Possibly, CLEF notation becomes sufficiently performant
    
    // Calculate vertex first, since it will be used several times. 
    // Be careful about the vertex band labels here so that they can be used safely later on.    
    // Slightly tricky: different Matsubara grids for _gimp and _G2_iw, have to use _gimp[0](n1) instead of _gimp[0][n1]
    
    //TODO: Timing/memory: is it necessary to precalculate the vertex?
    vertex() = 0.;
    for (auto const &s2 : range(_gimp.size())) {
      for (auto const &s1 : range(_gimp.size())) {
        for (const auto &[iw, n1, n2] : _G2_iw(s1,s2).mesh()){
          for (const auto [Ai, Aj, Ak, Al] : _G2_iw(s1,s2).target_indices()){

            //TODO: ensure that all matrix labels are correct                
            vertex(s1,s2)[iw, n1, n2](Ai, Aj, Ak, Al) 
            = 
            _G2_iw(s1,s2)[iw, n1, n2](Ai, Aj, Ak, Al)
            +beta*kronecker(s1,s2)*kronecker(n1,n2)*_gimp[s1](n1)(Al,Ai)*_gimp[s2](n1+iw)(Aj,Ak) 
            - beta* _gimp[s1](n1)(Aj,Ai)*_gimp[s2](n2)(Al,Ak)* kronecker(iw) ;
                
            
          }
        }
      }
    }   
    if (params.verbosity >= 5) h5_write(h5::file("vertex.h5",'w'),"vertex",vertex);
    
    if(params.calculate_sigma1){
    if (params.verbosity >= 2) std::cout << "Calculating Sigma dual: Calculate first order diagram" << std::endl;
    // First-order diagram    
    // Outer loop over s2, the spin of Sigma_dual
    // Inner loop over s1, the spin of the internal Green's function
    // R0 is (0,0,0) in real space
    // i,j,k,l are the band indices
    // iw,n1,n2 are the frequency indices
    s2=0;
    for (auto const &bl : gf_struct) {
      s1=0;  
      for (auto const &bl2 : gf_struct) {
        for (const auto &[iw, n1, n2] : vertex(s1,s2).mesh()){
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
    if (params.verbosity >= 2) std::cout << "Calculating Sigma dual: Second-order diagram" << std::endl;
    // Second-order diagram
    // Outer loop over s2, the spin of Sigma_dual
    // Inner loop over s1, the spin of two of the internal Green's functions
    // R is the spatial argument of Sigma in real space, two Green's functions have R and 1 has -R
    // i,j,k,l are the band indices
    // iw,n1,n2 are the frequency indices
    // There are two spin configurations possible, these are handled separately here. TODO: rewrite order of loops to make this more efficient
    auto orbital_indices = vertex(0,0).target_indices();
    for (const auto &[iw, n1, n2] : vertex(0,0).mesh()){
      for (const auto R : rmesh){  
        for (const auto [Ai, Aj, Ak, Al] : orbital_indices){
          for (const auto [Bi, Bj, Bk, Bl] : orbital_indices){                                
            s2=0;
            for (auto const &bl : gf_struct) {
              // First spin configuration s1=-s2 
              s1= sigma_bar(s2);  
              sigmad_real[s2][n2,R](Al,Bk) += 
              -0.5
              * vertex(s1,s2)(n1-n2,n2+iw, n2)(Ai, Aj, Ak, Al) 
              * vertex(s1,s2)(n2-n1, n1+iw, n1)(Bi, Bj, Bk, Bl) 
              * gd_real[s1](n2+iw,R)(Ai,Bj)
              * gd_real[s1](n1+iw,-R)(Bi,Aj)
              * gd_real[s2](n1,R)(Ak,Bl)
              /(beta*beta) ;                
              // Second spin configuration: sum over s1  
              s1=0;  
              for (auto const &bl2 : gf_struct) {
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
    
    if (params.verbosity >= 2) std::cout << "Calculating Sigma dual: FT" << std::endl;
    sigmad = make_gf_from_fourier<1>(sigmad_real);
    if (params.verbosity >= 1) h5_write(h5::file("gd_real.h5",'w'),"gd_real",gd_real);
    if (params.verbosity >= 1) h5_write(h5::file("sigmad_real.h5",'w'),"sigmad_real",sigmad_real);
    if (params.verbosity >= 1) h5_write(h5::file("sigmad_k.h5",'w'),"sigmad_k",sigmad);
    }
    else{
      sigmad() = 0;            
    }
    
    /* Calculate G lattice and Gdlattice */
    
    if (params.verbosity >= 2) std::cout << "Calculating Glat" << std::endl;
    // Create the container for the lattice Green's function
    auto glat = G_iw_k_t{{ {beta,Fermion,n_iw} ,kmesh},gf_struct} ;  
    for (auto const &b : range(glat.size())) {
        for (const auto &[iw,k] : glat[b].mesh()){
            glat[b][iw,k] = 1./( 1./(_gimp[b][iw]+sigmad[b][iw,k]) + _Delta[b][iw] - _Hk[b](k)  );
            gd0[b][iw,k]  = gd0[b][iw,k]/( 1.-sigmad[b][iw,k]*gd0[b][iw,k]);
        }
    }
    // Write the lattice Green's function
    if (params.verbosity >= 1 && params.calculate_sigma) h5_write(h5::file("G_k.h5",'w'),"G_k",glat);


    
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
    if (params.verbosity >= 2) std::cout << "Evaluate new Delta" << std::endl;
    _Delta(wn_) << _Delta(wn_) + params.ksi_delta * gdloc(wn_) * _gimp(wn_)/gloc(wn_);
    
    
    if (params.verbosity >= 2) std::cout << "dualfermion finished\n" ;
    
    return ;
    
  }
} // namespace triqs_dualfermion
