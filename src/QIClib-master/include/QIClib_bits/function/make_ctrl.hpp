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

//*******************************************************************************

#ifdef QICLIB_USE_SERIAL_MAKE_CTRL
// USE SERIAL ALGORITHM

//*******************************************************************************

template <typename T1, typename TR = typename std::enable_if<
                         is_floating_point_var<trait::pT<T1> >::value,
                         arma::Mat<trait::eT<T1> > >::type>

inline TR make_ctrl(const T1& A1, arma::uvec ctrl, arma::uvec subsys,
                    arma::uvec dim) {
  const auto& A = _internal::as_Mat(A1);

  arma::uword d = ctrl.n_elem > 0 ? dim.at(ctrl.at(0) - 1) : 1;
  arma::uword N = arma::prod(dim);

#ifndef QICLIB_NO_DEBUG
  if (A.n_elem == 0)
    throw Exception("qic::make_ctrl", Exception::type::ZERO_SIZE);

  if (A.n_rows != A.n_cols)
    throw Exception("qic::make_ctrl", Exception::type::MATRIX_NOT_SQUARE);

  for (arma::uword i = 1; i < ctrl.n_elem; ++i)
    if (dim.at(ctrl.at(i) - 1) != d)
      throw Exception("qic::make_ctrl", Exception::type::DIMS_NOT_EQUAL);

  if (dim.n_elem == 0 || arma::any(dim == 0))
    throw Exception("qic::make_ctrl", Exception::type::INVALID_DIMS);

  if (arma::prod(dim(subsys - 1)) != A.n_rows)
    throw Exception("qic::make_ctrl", Exception::type::DIMS_MISMATCH_MATRIX);

  const arma::uvec ctrlsubsys = arma::join_cols(subsys, ctrl);

  if (ctrlsubsys.n_elem > dim.n_elem ||
      arma::unique(ctrlsubsys).eval().n_elem != ctrlsubsys.n_elem ||
      arma::any(ctrlsubsys > dim.n_elem) || arma::any(ctrlsubsys == 0))
    throw Exception("qic::make_ctrl", Exception::type::INVALID_SUBSYS);
#endif

  _internal::dim_collapse_sys_ctrl(dim, subsys, ctrl);

  const arma::uword n = dim.n_elem;
  const arma::uword m = subsys.n_elem;
  const arma::uword o = ctrl.n_elem;

  arma::uvec keep(n - m);
  arma::uword keep_count(0);
  for (arma::uword run = 0; run < n; ++run) {
    if (!arma::any(subsys == run + 1)) {
      keep.at(keep_count) = run + 1;
      ++keep_count;
    }
  }

  arma::uword product[_internal::MAXQDIT];
  product[n - 1] = 1;
  for (arma::uword i = 1; i < n; ++i)
    product[n - 1 - i] = product[n - i] * dim.at(n - i);

  arma::uword productr[_internal::MAXQDIT];
  productr[m - 1] = 1;
  for (arma::uword i = 1; i < m; ++i)
    productr[m - 1 - i] = productr[m - i] * dim.at(subsys.at(m - i) - 1);

  arma::uword p_num = std::max(static_cast<arma::uword>(1), d - 1);

  arma::field<arma::Mat<trait::eT<T1> > > Ap(p_num + 1);
  for (arma::uword i = 0; i <= p_num; ++i)
    Ap.at(i) = _internal::POWM_GEN_INT(A, i);

  arma::Mat<trait::eT<T1> > U(N, N, arma::fill::zeros);

  const arma::uword loop_no = 2 * n;
  constexpr auto loop_no_buffer = 2 * _internal::MAXQDIT + 1;
  arma::uword loop_counter[loop_no_buffer] = {0};
  arma::uword MAX[loop_no_buffer];

  for (arma::uword i = 0; i < n; ++i) {
    MAX[i] = dim.at(i);
    if (arma::any(keep == i + 1))
      MAX[i + n] = 1;
    else
      MAX[i + n] = dim.at(i);
  }
  MAX[loop_no] = 2;

  arma::uword p1 = 0;

  while (loop_counter[loop_no] == 0) {
    arma::uword count1(0), count2(0);
    for (arma::uword i = 0; i < n; ++i) {
      count1 += (arma::any(ctrl == i + 1) && loop_counter[i] != 0) ? 1 : 0;
      count2 += loop_counter[i + n] == 0 ? 1 : 0;
    }

    if ((count1 != o) && (count2 == n)) {
      arma::uword I(0);
      for (arma::uword i = 0; i < n; ++i) I += product[i] * loop_counter[i];
      U.at(I, I) = static_cast<trait::eT<T1> >(1.0);

    } else if (count1 == o) {
      arma::uword I(0), J(0), K(0), L(0);
      int power = o == 0 ? 1 : 0;

      for (arma::uword i = 0; i < n; ++i) {
        if (arma::any(keep == i + 1)) {
          I += product[i] * loop_counter[i];
          J += product[i] * loop_counter[i];

        } else {
          I += product[i] * loop_counter[i];
          J += product[i] * loop_counter[i + n];
        }

        arma::uword counter(0);
        while (any(subsys == i + 1)) {
          if (subsys.at(counter) != i + 1) {
            ++counter;
          } else {
            K += productr[counter] * loop_counter[i];
            L += productr[counter] * loop_counter[i + n];
            break;
          }
        }
      }

      if (o != 0) {
        arma::uword counter_1(1);
        for (arma::uword j = 1; j < o; ++j)
          counter_1 +=
            loop_counter[ctrl.at(0) - 1] == loop_counter[ctrl.at(j) - 1] ? 1
                                                                         : 0;
        power = counter_1 == o ? loop_counter[ctrl.at(0) - 1] : 0;
      }
      U.at(I, J) = Ap.at(power).at(K, L);
    }

    ++loop_counter[0];
    while (loop_counter[p1] == MAX[p1]) {
      loop_counter[p1] = 0;
      loop_counter[++p1]++;
      if (loop_counter[p1] != MAX[p1])
        p1 = 0;
    }
  }
  return U;
}

//*******************************************************************************

#else
// USE PARALLEL ALGORITHM

//*******************************************************************************

template <typename T1, typename TR = typename std::enable_if<
                         is_floating_point_var<trait::pT<T1> >::value,
                         arma::Mat<trait::eT<T1> > >::type>

inline TR make_ctrl(const T1& A1, arma::uvec ctrl, arma::uvec subsys,
                    arma::uvec dim) {
  const auto& A = _internal::as_Mat(A1);

  arma::uword d = ctrl.n_elem > 0 ? dim.at(ctrl.at(0) - 1) : 1;
  arma::uword N = arma::prod(dim);

#ifndef QICLIB_NO_DEBUG
  if (A.n_elem == 0)
    throw Exception("qic::make_ctrl", Exception::type::ZERO_SIZE);

  if (A.n_rows != A.n_cols)
    throw Exception("qic::make_ctrl", Exception::type::MATRIX_NOT_SQUARE);

  for (arma::uword i = 1; i < ctrl.n_elem; ++i)
    if (dim.at(ctrl.at(i) - 1) != d)
      throw Exception("qic::make_ctrl", Exception::type::DIMS_NOT_EQUAL);

  if (dim.n_elem == 0 || arma::any(dim == 0))
    throw Exception("qic::make_ctrl", Exception::type::INVALID_DIMS);

  if (arma::prod(dim(subsys - 1)) != A.n_rows)
    throw Exception("qic::make_ctrl", Exception::type::DIMS_MISMATCH_MATRIX);

  const arma::uvec ctrlsubsys = arma::join_cols(subsys, ctrl);

  if (ctrlsubsys.n_elem > dim.n_elem ||
      arma::unique(ctrlsubsys).eval().n_elem != ctrlsubsys.n_elem ||
      arma::any(ctrlsubsys > dim.n_elem) || arma::any(ctrlsubsys == 0))
    throw Exception("qic::make_ctrl", Exception::type::INVALID_SUBSYS);
#endif

  _internal::dim_collapse_sys_ctrl(dim, subsys, ctrl);

  const arma::uword n = dim.n_elem;
  const arma::uword m = subsys.n_elem;
  const arma::uword o = ctrl.n_elem;

  arma::uvec keep(n - m);
  arma::uword keep_count(0);
  for (arma::uword run = 0; run < n; ++run) {
    if (!arma::any(subsys == run + 1)) {
      keep.at(keep_count) = run + 1;
      ++keep_count;
    }
  }

  arma::uword productr[_internal::MAXQDIT];
  productr[m - 1] = 1;
  for (arma::uword i = 1; i < m; ++i)
    productr[m - 1 - i] = productr[m - i] * dim.at(subsys.at(m - i) - 1);

  arma::uword p_num = std::max(static_cast<arma::uword>(1), d - 1);

  arma::field<arma::Mat<trait::eT<T1> > > Ap(p_num + 1);
  for (arma::uword i = 0; i <= p_num; ++i)
    Ap.at(i) = _internal::POWM_GEN_INT(A, i);

  arma::Mat<trait::eT<T1> > U(N, N);

  auto worker = [n, o, &dim, &subsys, &ctrl, &keep, &productr,
                 &Ap](arma::uword I, arma::uword J) noexcept -> trait::eT<T1> {

    bool equality_check = I == J;
    arma::uword Iindex[_internal::MAXQDIT];
    arma::uword Jindex[_internal::MAXQDIT];
    arma::uword count1(0);

    for (arma::uword i = 1; i < n; ++i) {
      Iindex[n - i] = I % dim.at(n - i);
      Jindex[n - i] = J % dim.at(n - i);
      I /= dim.at(n - i);
      J /= dim.at(n - i);

      if (arma::any(keep == n - i + 1) && (Iindex[n - i] != Jindex[n - i]))
        return static_cast<trait::eT<T1> >(0);

      count1 += (arma::any(ctrl == n - i + 1) && Iindex[n - i] != 0) ? 1 : 0;
    }

    Iindex[0] = I;
    Jindex[0] = J;

    if (arma::any(keep == 1) && (Iindex[0] != Jindex[0]))
      return static_cast<trait::eT<T1> >(0);

    count1 += (arma::any(ctrl == 1) && Iindex[0] != 0) ? 1 : 0;

    if (equality_check && count1 != o)
      return static_cast<trait::eT<T1> >(1);
    else if (count1 != o)
      return static_cast<trait::eT<T1> >(0);

    arma::uword K(0), L(0);
    for (arma::uword i = 0; i < n; ++i) {
      arma::uword count2(0);
      while (arma::any(subsys == i + 1)) {
        if (subsys.at(count2) != i + 1) {
          ++count2;
        } else {
          K += productr[count2] * Iindex[i];
          L += productr[count2] * Jindex[i];
          break;
        }
      }
    }

    arma::uword power = o == 0 ? 1 : 0;
    if (o != 0) {
      arma::uword count3(1);
      for (arma::uword j = 1; j < o; ++j)
        count3 += Iindex[ctrl.at(0) - 1] == Jindex[ctrl.at(j) - 1] ? 1 : 0;
      power = count3 == o ? Jindex[ctrl.at(0) - 1] : 0;
    }

    return Ap.at(power).at(K, L);
  };

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (arma::uword JJ = 0; JJ < N; ++JJ) {
    for (arma::uword II = 0; II < N; ++II) {
      U.at(II, JJ) = worker(II, JJ);
    }
  }
  return U;
}

//*******************************************************************************

#endif

//*******************************************************************************

template <typename T1, typename TR = typename std::enable_if<
                         is_floating_point_var<trait::pT<T1> >::value,
                         arma::Mat<trait::eT<T1> > >::type>

inline TR make_ctrl(const T1& A1, arma::uvec ctrl, arma::uvec subsys,
                    arma::uword n, arma::uword dim = 2) {
  const auto& A = _internal::as_Mat(A1);

#ifndef QICLIB_NO_DEBUG
  if (n == 0)
    throw Exception("qic::make_ctrl", Exception::type::OUT_OF_RANGE);
  if (dim == 0)
    throw Exception("qic::make_ctrl", Exception::type::INVALID_DIMS);
#endif

  arma::uvec dim2(n);
  dim2.fill(dim);
  return make_ctrl(A, std::move(ctrl), std::move(subsys), std::move(dim2));
}

//*******************************************************************************

}  // namespace qic
