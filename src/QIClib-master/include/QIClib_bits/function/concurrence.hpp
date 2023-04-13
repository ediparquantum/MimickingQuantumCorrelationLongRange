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
            is_floating_point_var<trait::pT<T1> >::value, trait::pT<T1> >::type>

inline TR concurrence(const T1& rho1) {
  const auto& rho = _internal::as_Mat(rho1);

#ifndef QICLIB_NO_DEBUG
  if (rho.n_elem == 0)
    throw Exception("qic::concurrence", Exception::type::ZERO_SIZE);

  if (rho.n_rows != rho.n_cols)
    throw Exception("qic::concurrence", Exception::type::MATRIX_NOT_SQUARE);

  if (rho.n_rows != 4)
    throw Exception("qic::concurrence", Exception::type::NOT_QUBIT_SUBSYS);
#endif

  auto& S2 = SPM<trait::pT<T1> >::get_instance().S.at(2);

  typename arma::Mat<std::complex<trait::pT<T1> > >::template fixed<4, 4> pbar =
    rho * arma::kron(S2, S2) * arma::conj(rho) * arma::kron(S2, S2);

  typename arma::Col<trait::pT<T1> >::template fixed<4> eig =
    arma::sort(arma::real(arma::eig_gen(pbar)));

  for (auto& i : eig) {
    if (i < _precision::eps<trait::pT<T1> >::value)
      i = 0.0;
  }

  return std::max(static_cast<trait::pT<T1> >(0.0),
                  std::sqrt(eig.at(3)) - std::sqrt(eig.at(2)) -
                    std::sqrt(eig.at(1)) - std::sqrt(eig.at(0)));
}

//******************************************************************************

}  // namespace qic
