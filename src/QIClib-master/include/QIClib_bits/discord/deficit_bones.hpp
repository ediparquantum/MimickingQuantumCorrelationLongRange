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
          typename Enable =
            typename std::enable_if<arma::is_Mat_only<T1>::value, void>::type>
class deficit_space;

template <typename T1> class deficit_space<T1> {

  // static_assert(
  // std::is_same<T1, arma::Mat<trait::eT<T1> > >::value,
  // "deficit_space requires Armadillo Mat object as template argument!");

 private:
  T1 _rho{};
  arma::uword _subsys{};
  arma::uword _n_cols{};
  arma::uword _n_rows{};

  arma::uword _party_no{};
  arma::uvec _dim{};

  trait::pT<T1> _S_A_B{};
  trait::pT<T1> _result{};
  arma::Col<trait::pT<T1> > _tp{};
  trait::pT<T1> _result_reg{};
  arma::Col<trait::pT<T1> > _result_reg_all{};

  bool _is_computed{false};
  bool _is_reg_computed{false};
  bool _is_sab_computed{false};
  bool _deficit2{false};
  bool _deficit3{false};

  nlopt::algorithm _deficit_global_opt{};
  double _deficit_global_xtol{};
  double _deficit_global_ftol{};
  bool _deficit_global{false};
  nlopt::algorithm _deficit_local_opt{};
  double _deficit_local_xtol{};
  double _deficit_local_ftol{};
  arma::vec _deficit_angle_range{};
  arma::vec _deficit_angle_ini{};

  inline void init(arma::uvec);
  inline void s_a_b();
  inline void default_setting();

  //****************************************************************************

 public:
  //****************************************************************************

  deficit_space() = delete;
  deficit_space(const deficit_space&) = delete;
  deficit_space& operator=(const deficit_space&) = delete;
  ~deficit_space() = default;

  //****************************************************************************

  inline deficit_space(const T1&, arma::uword, arma::uvec);
  inline deficit_space(T1&&, arma::uword, arma::uvec);
  inline deficit_space(const T1&, arma::uword, arma::uword = 2);
  inline deficit_space(T1&&, arma::uword, arma::uword = 2);

  //****************************************************************************

  inline deficit_space& global_algorithm(nlopt::algorithm) noexcept;
  inline deficit_space& global_xtol(double) noexcept;
  inline deficit_space& global_ftol(double) noexcept;
  inline deficit_space& use_global_opt(bool) noexcept;
  inline deficit_space& local_algorithm(nlopt::algorithm) noexcept;
  inline deficit_space& local_xtol(double) noexcept;
  inline deficit_space& local_ftol(double) noexcept;
  inline deficit_space& angle_range(const arma::vec&);
  inline deficit_space& initial_angle(const arma::vec&);

  //****************************************************************************

  inline deficit_space& compute();
  inline deficit_space& compute_reg();
  inline const arma::Col<trait::pT<T1> >& opt_angles();
  inline const trait::pT<T1>& result();
  inline const trait::pT<T1>& result_reg();
  inline const arma::Col<trait::pT<T1> >& result_reg_all();
  inline deficit_space& reset();
  inline deficit_space& reset(arma::uword);
  inline deficit_space& reset(const T1&, arma::uword, arma::uvec);
  inline deficit_space& reset(T1&&, arma::uword, arma::uvec);
  inline deficit_space& reset(const T1&, arma::uword, arma::uword = 2);
  inline deficit_space& reset(T1&&, arma::uword, arma::uword = 2);
};

//******************************************************************************

}  // namespace qic
