/*
 * QIClib (Quantum information and computation library)
 *
 * Copyright (c) 2015 - 2017  Titas Chanda (titas.chanda@gmail.com)
 *
 * This file is part of QIClib.
 *
 * QIClib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QIClib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QIClib.  If not, see <http://www.gnu.org/licenses/>.
 */

namespace qic {

//******************************************************************************

template <typename T1,
          typename TR = typename std::enable_if<
            is_floating_point_var<trait::pT<T1> >::value,
            typename arma::Col<trait::pT<T1> >::template fixed<3> >::type>
TR deficit_reg(const T1& rho1, arma::uword nodal, arma::uvec dim) {
  const auto& rho = _internal::as_Mat(rho1);
  arma::uword party_no = dim.n_elem;
  arma::uword dim1 = arma::prod(dim);

#ifndef QICLIB_NO_DEBUG
  if (rho.n_elem == 0)
    throw Exception("qic::deficit_reg", Exception::type::ZERO_SIZE);

  if (rho.n_rows != rho.n_cols)
    throw Exception("qic::deficit_reg", Exception::type::MATRIX_NOT_SQUARE);

  if (any(dim == 0))
    throw Exception("qic::deficit_reg", Exception::type::INVALID_DIMS);

  if (dim1 != rho.n_rows)
    throw Exception("qic::deficit_reg", Exception::type::DIMS_MISMATCH_MATRIX);

  if (nodal <= 0 || nodal > party_no)
    throw Exception("qic::deficit_reg", "Invalid measured party index");

  if (dim(nodal - 1) != 2)
    throw Exception("qic::deficit_reg", "Measured party is not qubit");
#endif

  auto S_A_B = entropy(rho);

  dim1 /= 2;
  arma::uword dim2(1);
  for (arma::uword i = 0; i < nodal - 1; ++i) dim2 *= dim.at(i);
  arma::uword dim3(1);
  for (arma::uword i = nodal; i < party_no; ++i) dim3 *= dim.at(i);

  arma::Mat<trait::pT<T1> > eye2 =
    arma::eye<arma::Mat<trait::pT<T1> > >(dim1, dim1);
  arma::Mat<trait::pT<T1> > eye3 =
    arma::eye<arma::Mat<trait::pT<T1> > >(dim2, dim2);
  arma::Mat<trait::pT<T1> > eye4 =
    arma::eye<arma::Mat<trait::pT<T1> > >(dim3, dim3);

  typename arma::Col<trait::pT<T1> >::template fixed<3> ret;

  for (arma::uword i = 0; i < 3; ++i) {
    arma::Mat<std::complex<trait::pT<T1> > > proj1 =
      SPM<trait::pT<T1> >::get_instance().proj2.at(0, i + 1);

    arma::Mat<std::complex<trait::pT<T1> > > proj2 =
      SPM<trait::pT<T1> >::get_instance().proj2.at(1, i + 1);

    if (nodal == 1) {
      proj1 = kron(proj1, eye2);
      proj2 = kron(proj2, eye2);
    } else if (party_no == nodal) {
      proj1 = kron(eye2, proj1);
      proj2 = kron(eye2, proj2);
    } else {
      proj1 = kron(kron(eye3, proj1), eye4);
      proj2 = kron(kron(eye3, proj2), eye4);
    }

    arma::Mat<std::complex<trait::pT<T1> > > rho_1 = (proj1 * rho * proj1);
    arma::Mat<std::complex<trait::pT<T1> > > rho_2 = (proj2 * rho * proj2);

    rho_1 += rho_2;
    trait::pT<T1> S_max = entropy(rho_1);
    ret.at(i) = -S_A_B + S_max;
  }
  return ret;
}

//******************************************************************************

}  // namespace qic
