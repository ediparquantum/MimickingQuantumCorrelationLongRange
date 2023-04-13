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

template <typename T1, typename T2,
          typename TR = typename std::enable_if<
            is_arma_type_var<T1, T2>::value && is_same_pT_var<T1, T2>::value,
            trait::pT<T1> >::type>

inline TR HS_dist(const T1& rho11, const T2& rho12) {
  const auto& rho1 = _internal::as_Mat(rho11);
  const auto& rho2 = _internal::as_Mat(rho12);

#ifndef QICLIB_NO_DEBUG
  if (rho1.n_elem == 0 || rho2.n_elem == 0)
    throw Exception("qic::HS_dist", Exception::type::ZERO_SIZE);

  if (rho1.n_rows != rho1.n_cols || rho2.n_rows != rho2.n_cols)
    throw Exception("qic::HS_dist", Exception::type::MATRIX_NOT_SQUARE);

  if (rho1.n_rows != rho2.n_rows)
    throw Exception("qic::HS_dist", Exception::type::MATRIX_SIZE_MISMATCH);
#endif

  auto rho3 = (rho1 - rho2).eval();
  return std::real(arma::trace(rho3.t() * rho3));
}

//******************************************************************************

template <typename T1, typename T2,
          typename TR = typename std::enable_if<
            is_floating_point_var<trait::pT<T1>, trait::pT<T2> >::value &&
              is_same_pT_var<T1, T2>::value,
            trait::pT<T1> >::type>

inline TR tr_dist(const T1& rho11, const T2& rho12) {
  const auto& rho1 = _internal::as_Mat(rho11);
  const auto& rho2 = _internal::as_Mat(rho12);

#ifndef QICLIB_NO_DEBUG
  if (rho1.n_elem == 0 || rho2.n_elem == 0)
    throw Exception("qic::tr_dist", Exception::type::ZERO_SIZE);

  if (rho1.n_rows != rho1.n_cols || rho2.n_rows != rho2.n_cols)
    throw Exception("qic::tr_dist", Exception::type::MATRIX_NOT_SQUARE);

  if (rho1.n_rows != rho2.n_rows)
    throw Exception("qic::tr_dist", Exception::type::MATRIX_SIZE_MISMATCH);
#endif

  auto rho3 = (rho1 - rho2).eval();
  auto eig1 = arma::eig_sym(rho3);
  return arma::sum(arma::abs(eig1)) * 0.5;
}

//******************************************************************************

template <typename T1, typename T2,
          typename TR = typename std::enable_if<
            is_floating_point_var<trait::pT<T1>, trait::pT<T2> >::value &&
              is_same_pT_var<T1, T2>::value,
            trait::pT<T1> >::type>

inline TR fidelity(const T1& rho11, const T2& rho12) {
  const auto& rho1 = _internal::as_Mat(rho11);
  const auto& rho2 = _internal::as_Mat(rho12);

#ifndef QICLIB_NO_DEBUG
  if (rho1.n_elem == 0 || rho2.n_elem == 0)
    throw Exception("qic::fidelity", Exception::type::ZERO_SIZE);

  if (rho1.n_rows != rho1.n_cols || rho2.n_rows != rho2.n_cols)
    throw Exception("qic::fidelity", Exception::type::MATRIX_NOT_SQUARE);

  if (rho1.n_rows != rho2.n_rows)
    throw Exception("qic::fidelity", Exception::type::MATRIX_SIZE_MISMATCH);
#endif

  auto rho3 = sqrtm_sym((sqrtm_sym(rho1) * rho2 * sqrtm_sym(rho1)).eval());
  return std::pow(std::real(arma::trace(rho3)), 2);
}

//******************************************************************************

template <typename T1, typename T2,
          typename TR = typename std::enable_if<
            is_floating_point_var<trait::pT<T1>, trait::pT<T2> >::value &&
              is_same_pT_var<T1, T2>::value,
            trait::pT<T1> >::type>

inline TR Bures_dist(const T1& rho11, const T2& rho12) {
  const auto& rho1 = _internal::as_Mat(rho11);
  const auto& rho2 = _internal::as_Mat(rho12);

#ifndef QICLIB_NO_DEBUG
  if (rho1.n_elem == 0 || rho2.n_elem == 0)
    throw Exception("qic::bures_dist", Exception::type::ZERO_SIZE);

  if (rho1.n_rows != rho1.n_cols || rho2.n_rows != rho2.n_cols)
    throw Exception("qic::bures_dist", Exception::type::MATRIX_NOT_SQUARE);

  if (rho1.n_rows != rho2.n_rows)
    throw Exception("qic::bures_dist", Exception::type::MATRIX_SIZE_MISMATCH);
#endif
  auto rho3 = sqrtm_sym((sqrtm_sym(rho1) * rho2 * sqrtm_sym(rho1)).eval());
  auto fid = std::real(arma::trace(rho3));
  return std::real(std::sqrt(
    static_cast<std::complex<decltype(fid)> >(2.0 - 2.0 * std::sqrt(fid))));
}

//******************************************************************************

}  // namespace qic
