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

#include <fstream>

inline bool exists(std::string filename) {
  std::ifstream f(filename.c_str());
  bool e = f.good();
  if (e) f.close();
  return e;
}

template <class T>
int kronecker (T a, T b) {
    if(std::complex<double>(a)==std::complex<double>(b)) return 1;
    return 0;
}

template <class T>
int kronecker (T a) {
    if(std::complex<double>(a)==0.) return 1;
    return 0;
}

// Does not work with more than two blocks
// 1 --> 0
// 0 --> 1
// int sigma_bar(int sigma) {
//    return 1 - sigma ;
//}
