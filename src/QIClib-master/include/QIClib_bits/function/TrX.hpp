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

#ifdef QICLIB_USE_SERIAL_TRX
// USE SERIAL ALGORITHM

//******************************************************************************

template <typename T1,
          typename TR = typename std::enable_if<
            is_arma_type_var<T1>::value, arma::Mat<trait::eT<T1> > >::type>

inline TR TrX(const T1& rho1, arma::uvec subsys, arma::uvec dim,
              bool is_Hermitian = false) {
  const auto& rho = _internal::as_Mat(rho1);
  bool checkV = (rho.n_cols != 1);

#ifndef QICLIB_NO_DEBUG
  if (rho.n_elem == 0)
    throw Exception("qic::TrX", Exception::type::ZERO_SIZE);

  if (checkV)
    if (rho.n_rows != rho.n_cols)
      throw Exception("qic::TrX",
                      Exception::type::MATRIX_NOT_SQUARE_OR_CVECTOR);

  if (dim.n_elem == 0 || arma::any(dim == 0))
    throw Exception("qic::TrX", Exception::type::INVALID_DIMS);

  if (arma::prod(dim) != rho.n_rows)
    throw Exception("qic::TrX", Exception::type::DIMS_MISMATCH_MATRIX);

  if (dim.n_elem < subsys.n_elem || arma::any(subsys == 0) ||
      arma::any(subsys > dim.n_elem) ||
      subsys.n_elem != arma::unique(subsys).eval().n_elem)
    throw Exception("qic::TrX", Exception::type::INVALID_SUBSYS);
#endif

  if (subsys.n_elem == dim.n_elem) {
    if (checkV)
      return {arma::trace(rho)};
    else
      return {rho.t() * rho};
  }

  if (subsys.n_elem == 0) {
    if (checkV)
      return rho;
    else
      return rho * rho.t();
  }

  _internal::dim_collapse_sys(dim, subsys);
  const arma::uword n = dim.n_elem;
  const arma::uword m = subsys.n_elem;

  arma::uvec keep(n - m);
  arma::uword keep_count(0);
  for (arma::uword run = 0; run < n; ++run) {
    if (!arma::any(subsys == run + 1)) {
      keep.at(keep_count) = run + 1;
      ++keep_count;
    }
  }

  arma::uword dimtrace = arma::prod(dim(subsys - 1));
  arma::uword dimkeep = rho.n_rows / dimtrace;

  arma::uword product[_internal::MAXQDIT];
  product[n - 1] = 1;
  for (arma::uword i = 1; i < n; ++i)
    product[n - 1 - i] = product[n - i] * dim.at(n - i);

  arma::uword productr[_internal::MAXQDIT];
  productr[n - m - 1] = 1;
  for (arma::uword i = 1; i < n - m; ++i)
    productr[n - m - 1 - i] =
      productr[n - m - i] * dim.at(keep.at(n - m - i) - 1);

  arma::Mat<trait::eT<T1> > tr_rho(dimkeep, dimkeep, arma::fill::zeros);

  const arma::uword loop_no = 2 * n;
  constexpr auto loop_no_buffer = 2 * _internal::MAXQDIT + 1;
  arma::uword loop_counter[loop_no_buffer] = {0};
  arma::uword MAX[loop_no_buffer];

  for (arma::uword i = 0; i < n; ++i) {
    MAX[i] = dim.at(i);
    if (arma::any(subsys == (i + 1)))
      MAX[i + n] = 1;
    else
      MAX[i + n] = dim.at(i);
  }
  MAX[loop_no] = 2;

  arma::uword p1 = 0;

  while (loop_counter[loop_no] == 0) {
    arma::uword I(0), J(0), K(0), L(0), n_to_k(0);

    for (arma::uword i = 0; i < n; ++i) {
      if (arma::any(subsys == i + 1)) {
        auto tmp = product[i] * loop_counter[i];
        I += tmp;
        J += tmp;

      } else {
        I += product[i] * loop_counter[i];
        J += product[i] * loop_counter[i + n];
      }

      if (arma::any(keep == i + 1)) {
        K += productr[n_to_k] * loop_counter[i];
        L += productr[n_to_k] * loop_counter[i + n];
        ++n_to_k;
      }
    }

    if ((is_Hermitian && L >= K) || (!is_Hermitian))
      tr_rho.at(K, L) +=
        checkV ? rho.at(I, J) : rho.at(I) * _internal::conj2(rho.at(J));

    ++loop_counter[0];
    while (loop_counter[p1] == MAX[p1]) {
      loop_counter[p1] = 0;
      loop_counter[++p1]++;
      if (loop_counter[p1] != MAX[p1])
        p1 = 0;
    }
  }

  if (is_Hermitian) {
    for (arma::uword L = 0; L < dimkeep; ++L)
      for (arma::uword K = L + 1; K < dimkeep; ++K)
        tr_rho.at(K, L) = _internal::conj2(tr_rho.at(L, K));
  }

  return tr_rho;
}

//******************************************************************************

#else
// USE PARALLEL ALGORITHM

//******************************************************************************

template <typename T1,
          typename TR = typename std::enable_if<
            is_arma_type_var<T1>::value, arma::Mat<trait::eT<T1> > >::type>

inline TR TrX(const T1& rho1, arma::uvec subsys, arma::uvec dim,
              bool is_Hermitian = false) {
  const auto& rho = _internal::as_Mat(rho1);
  bool checkV = (rho.n_cols != 1);

#ifndef QICLIB_NO_DEBUG
  if (rho.n_elem == 0)
    throw Exception("qic::TrX", Exception::type::ZERO_SIZE);

  if (checkV)
    if (rho.n_rows != rho.n_cols)
      throw Exception("qic::TrX",
                      Exception::type::MATRIX_NOT_SQUARE_OR_CVECTOR);

  if (dim.n_elem == 0 || arma::any(dim == 0))
    throw Exception("qic::TrX", Exception::type::INVALID_DIMS);

  if (arma::prod(dim) != rho.n_rows)
    throw Exception("qic::TrX", Exception::type::DIMS_MISMATCH_MATRIX);

  if (dim.n_elem < subsys.n_elem || arma::any(subsys == 0) ||
      arma::any(subsys > dim.n_elem) ||
      subsys.n_elem != arma::unique(subsys).eval().n_elem)
    throw Exception("qic::TrX", Exception::type::INVALID_SUBSYS);
#endif

  if (subsys.n_elem == dim.n_elem) {
    if (checkV)
      return {arma::trace(rho)};
    else
      return {rho.t() * rho};
  }

  if (subsys.n_elem == 0) {
    if (checkV)
      return rho;
    else
      return rho * rho.t();
  }

  _internal::dim_collapse_sys(dim, subsys);
  const arma::uword n = dim.n_elem;
  const arma::uword m = subsys.n_elem;

  arma::uword keep[_internal::MAXQDIT];
  arma::uword keep_count(0);
  for (arma::uword run = 0; run < n; ++run) {
    if (!arma::any(subsys == run + 1)) {
      keep[keep_count] = run + 1;
      ++keep_count;
    }
  }

  arma::uword dimtrace = arma::prod(dim(subsys - 1));
  arma::uword dimkeep = rho.n_rows / dimtrace;

  arma::uword product[_internal::MAXQDIT];
  product[n - 1] = 1;
  for (arma::uword i = 1; i < n; ++i)
    product[n - 1 - i] = product[n - i] * dim.at(n - i);

  arma::Mat<trait::eT<T1> > tr_rho(dimkeep, dimkeep);

  auto worker = [n, m, checkV, &dim, &keep, &subsys, &product,
                 &rho](arma::uword K, arma::uword L) noexcept -> trait::eT<T1> {

    arma::uword Kindex[_internal::MAXQDIT];
    arma::uword Lindex[_internal::MAXQDIT];

    for (arma::uword i = 1; i < n - m; ++i) {
      Kindex[n - m - i] = K % dim.at(keep[n - m - i] - 1);
      Lindex[n - m - i] = L % dim.at(keep[n - m - i] - 1);
      K /= dim.at(keep[n - m - i] - 1);
      L /= dim.at(keep[n - m - i] - 1);
    }
    Kindex[0] = K;
    Lindex[0] = L;

    trait::eT<T1> ret = static_cast<trait::eT<T1> >(0);

    const arma::uword loop_no = m;
    constexpr auto loop_no_buffer = _internal::MAXQDIT + 1;
    arma::uword loop_counter[loop_no_buffer] = {0};
    arma::uword MAX[loop_no_buffer];

    for (arma::uword i = 0; i < m; ++i) {
      MAX[i] = dim.at(subsys.at(i) - 1);
    }
    MAX[loop_no] = 2;
    arma::uword p1 = 0;

    while (loop_counter[loop_no] == 0) {
      arma::uword I(0), J(0), countK(0), countI(0);

      for (arma::uword i = 0; i < n; ++i) {
        if (arma::any(subsys == i + 1)) {
          auto tmp = product[i] * loop_counter[countI];
          I += tmp;
          J += tmp;
          ++countI;
        } else {
          I += product[i] * Kindex[countK];
          J += product[i] * Lindex[countK];
          ++countK;
        }
      }

      ret += checkV ? rho.at(I, J) : rho.at(I) * _internal::conj2(rho.at(J));

      ++loop_counter[0];
      while (loop_counter[p1] == MAX[p1]) {
        loop_counter[p1] = 0;
        loop_counter[++p1]++;
        if (loop_counter[p1] != MAX[p1])
          p1 = 0;
      }
    }

    return ret;
  };

  if (is_Hermitian) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (arma::uword LL = 0; LL < dimkeep; ++LL) {
      for (arma::uword KK = 0; KK <= LL; ++KK)
        tr_rho.at(KK, LL) = worker(KK, LL);
    }

    for (arma::uword LL = 0; LL < dimkeep; ++LL) {
      for (arma::uword KK = LL + 1; KK < dimkeep; ++KK)
        tr_rho.at(KK, LL) = _internal::conj2(tr_rho.at(LL, KK));
    }

  } else {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (arma::uword LL = 0; LL < dimkeep; ++LL) {
      for (arma::uword KK = 0; KK < dimkeep; ++KK)
        tr_rho.at(KK, LL) = worker(KK, LL);
    }
  }
  return tr_rho;
}

//******************************************************************************

#endif

//******************************************************************************

template <typename T1,
          typename TR = typename std::enable_if<
            is_arma_type_var<T1>::value, arma::Mat<trait::eT<T1> > >::type>

inline TR TrX(const T1& rho1, arma::uvec subsys, arma::uword dim = 2,
              bool is_Hermitian = false) {
  const auto& rho = _internal::as_Mat(rho1);

#ifndef QICLIB_NO_DEBUG
  bool checkV = (rho.n_cols != 1);

  if (rho.n_elem == 0)
    throw Exception("qic::TrX", Exception::type::ZERO_SIZE);

  if (checkV)
    if (rho.n_rows != rho.n_cols)
      throw Exception("qic::TrX",
                      Exception::type::MATRIX_NOT_SQUARE_OR_CVECTOR);

  if (dim == 0)
    throw Exception("qic::TrX", Exception::type::INVALID_DIMS);
#endif

  arma::uword n = static_cast<arma::uword>(
    QICLIB_ROUND_OFF(std::log(rho.n_rows) / std::log(dim)));

  arma::uvec dim2(n);
  dim2.fill(dim);
  return TrX(rho, std::move(subsys), std::move(dim2), is_Hermitian);
}

//******************************************************************************

namespace experimental {

//******************************************************************************

template <typename T1,
          typename TR = typename std::enable_if<
            is_arma_type_var<T1>::value, arma::Mat<trait::eT<T1> > >::type>

inline TR TrX(const T1& rho1, const arma::uvec& Sbasis, arma::uvec subsys,
              arma::uvec dim) {
  const auto& rho = _internal::as_Mat(rho1);
  bool checkV = (rho.n_cols != 1);

#ifndef QICLIB_NO_DEBUG
  if (rho.n_elem == 0)
    throw Exception("qic::TrX", Exception::type::ZERO_SIZE);

  if (checkV)
    if (rho.n_rows != rho.n_cols)
      throw Exception("qic::TrX",
                      Exception::type::MATRIX_NOT_SQUARE_OR_CVECTOR);

  if (dim.n_elem == 0 || arma::any(dim == 0))
    throw Exception("qic::TrX", Exception::type::INVALID_DIMS);

  //  if (arma::prod(dim) != p.n_rows)
  //  throw Exception("qic::TrX", Exception::type::DIMS_MISMATCH_MATRIX);

  if (dim.n_elem < subsys.n_elem || arma::any(subsys == 0) ||
      arma::any(subsys > dim.n_elem) ||
      subsys.n_elem != arma::unique(subsys).eval().n_elem)
    throw Exception("qic::TrX", Exception::type::INVALID_SUBSYS);

  if ((Sbasis.n_elem != rho.n_rows) || arma::any(Sbasis > arma::prod(dim)) ||
      Sbasis.n_elem != arma::unique(Sbasis).eval().n_elem)
    throw Exception("qic::TrX", "Invalid basis!");
#endif

  if (subsys.n_elem == dim.n_elem) {
    if (checkV)
      return {arma::trace(rho)};
    else
      return {rho.t() * rho};
  }

  if (subsys.n_elem == 0) {
    if (checkV)
      return rho;
    else
      return rho * rho.t();
  }

  _internal::dim_collapse_sys(dim, subsys);
  const arma::uword n = dim.n_elem;
  const arma::uword m = subsys.n_elem;

  arma::uvec keep(n - m);
  arma::uword keep_count(0);
  for (arma::uword run = 0; run < n; ++run) {
    if (!arma::any(subsys == run + 1)) {
      keep.at(keep_count) = run + 1;
      ++keep_count;
    }
  }

  arma::uword dimtrace = arma::prod(dim(subsys - 1));
  arma::uword dimkeep = arma::prod(dim) / dimtrace;

  arma::uword product[_internal::MAXQDIT];
  product[n - 1] = 1;
  for (arma::sword i = n - 2; i > -1; --i)
    product[i] = product[i + 1] * dim.at(i + 1);

  arma::uword productr[_internal::MAXQDIT];
  productr[n - m - 1] = 1;
  for (arma::sword i = n - m - 2; i > -1; --i)
    productr[i] = productr[i + 1] * dim.at(keep.at(i + 1) - 1);

  arma::Mat<trait::eT<T1> > tr_rho(dimkeep, dimkeep, arma::fill::zeros);

  const arma::uword loop_no = 2 * n;
  constexpr auto loop_no_buffer = 2 * _internal::MAXQDIT + 1;
  arma::uword loop_counter[loop_no_buffer] = {0};
  arma::uword MAX[loop_no_buffer];

  for (arma::uword i = 0; i < n; ++i) {
    MAX[i] = dim.at(i);
    if (arma::any(subsys == (i + 1)))
      MAX[i + n] = 1;
    else
      MAX[i + n] = dim.at(i);
  }
  MAX[loop_no] = 2;

  arma::uword p1 = 0;

  while (loop_counter[loop_no] == 0) {
    arma::uword I(0), J(0), K(0), L(0), n_to_k(0);
    arma::uword Icount(0), Jcount(0);

    for (arma::uword i = 0; i < n; ++i) {
      if (arma::any(subsys == i + 1)) {
        I += product[i] * loop_counter[i];
        J += product[i] * loop_counter[i];

      } else {
        I += product[i] * loop_counter[i];
        J += product[i] * loop_counter[i + n];
      }
    }

    while (Icount < Sbasis.n_elem && I != Sbasis.at(Icount)) ++Icount;

    while (Jcount < Sbasis.n_elem && J != Sbasis.at(Jcount)) ++Jcount;

    if (Icount < Sbasis.n_elem && Jcount < Sbasis.n_elem) {

      for (arma::uword i = 0; i < n; ++i) {
        if (arma::any(keep == i + 1)) {
          K += productr[n_to_k] * loop_counter[i];
          L += productr[n_to_k] * loop_counter[i + n];
          ++n_to_k;
        }
      }

      tr_rho.at(K, L) += checkV
                           ? rho.at(Icount, Jcount)
                           : rho.at(Icount) * _internal::conj2(rho.at(Jcount));
    }

    ++loop_counter[0];
    while (loop_counter[p1] == MAX[p1]) {
      loop_counter[p1] = 0;
      loop_counter[++p1]++;
      if (loop_counter[p1] != MAX[p1])
        p1 = 0;
    }
  }

  return tr_rho;
}

//******************************************************************************

}  // namespace experimental

//******************************************************************************

}  // namespace qic
