/** @file
 *****************************************************************************

 Implementation of interfaces for the "geometric sequence" evaluation domain.

 See geometric_sequence_domain.hpp .

 *****************************************************************************
 * @author     This file is part of libfqfft, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef GEOMETRIC_SEQUENCE_DOMAIN_TCC_
#define GEOMETRIC_SEQUENCE_DOMAIN_TCC_

#include <libfqfft/evaluation_domain/domains/basic_radix2_domain_aux.hpp>
#include <libfqfft/polynomial_arithmetic/basis_change.hpp>
#include <libff/common/profiling.hpp>
#include <fstream>

#ifdef MULTICORE
#include <omp.h>
#endif

namespace libfqfft {

template<typename FieldT>
geometric_sequence_domain<FieldT>::geometric_sequence_domain(const size_t m) : evaluation_domain<FieldT>(m + 1)
{
  if (m <= 1) throw InvalidSizeException("geometric(): expected m > 1");

  do_precomputation(m, FieldT::geometric_generator());
}

template<typename FieldT>
geometric_sequence_domain<FieldT>::geometric_sequence_domain(const size_t m, const FieldT &generator) : evaluation_domain<FieldT>(m + 1)
{
  if (m <= 1) throw InvalidSizeException("geometric(): expected m > 1");

  do_precomputation(m, generator);
}

template<typename FieldT>
void geometric_sequence_domain<FieldT>::FFT(std::vector<FieldT> &a)
{
  cosetFFT(a, FieldT::one());
}

template<typename FieldT>
void geometric_sequence_domain<FieldT>::iFFT(std::vector<FieldT> &a)
{
  icosetFFT(a, FieldT::one());
}

template<typename FieldT>
void geometric_sequence_domain<FieldT>::cosetFFT(std::vector<FieldT> &a, const FieldT &s)
{
  if (a.size() != this->m) throw DomainSizeException("geometric: expected a.size() == this->m");

  monomial_to_newton_basis_geometric(a, this->geometric_sequence, this->geometric_triangular_sequence, this->m, s);

  /* Newton to Evaluation */
  std::vector<FieldT> T(this->m);
  T[0] = FieldT::one();

  std::vector<FieldT> g(this->m);
  g[0] = a[0];

  FieldT scale = s;
  for (size_t i = 1; i < this->m; i++)
  {
    T[i] = T[i-1] * (this->geometric_sequence[i] - FieldT::one()).inverse();
    g[i] = this->geometric_triangular_sequence[i] * a[i] * scale;

    scale *= s;
  }

  _polynomial_multiplication(a, g, T);
  a.resize(this->m);

#ifdef MULTICORE
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < this->m; i++)
  {
    a[i] *= T[i].inverse();
  }
}

template<typename FieldT>
void geometric_sequence_domain<FieldT>::icosetFFT(std::vector<FieldT> &a, const FieldT &s)
{
  if (a.size() != this->m) throw DomainSizeException("geometric: expected a.size() == this->m");

  /* Interpolation to Newton */
  std::vector<FieldT> T(this->m);
  T[0] = FieldT::one();

  std::vector<FieldT> W(this->m);
  W[0] = a[0] * T[0];

  FieldT prev_T = T[0];
  for (size_t i = 1; i < this->m; i++)
  {
    prev_T *= (this->geometric_sequence[i] - FieldT::one()).inverse();

    W[i] = a[i] * prev_T;
    T[i] = this->geometric_triangular_sequence[i] * prev_T;
    if (i % 2 == 1) T[i] = -T[i];
  }

  _polynomial_multiplication(a, W, T);
  a.resize(this->m);

  FieldT s_inv = s.inverse();
  FieldT scale = FieldT::one();

#ifdef MULTICORE
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < this->m; i++)
  {
    a[i] *= this->geometric_triangular_sequence[i].inverse() * scale;
    scale *= s_inv;
  }

  newton_to_monomial_basis_geometric(a, this->geometric_sequence, this->geometric_triangular_sequence, this->m, s);
}

template<typename FieldT>
std::vector<FieldT> geometric_sequence_domain<FieldT>::evaluate_all_lagrange_polynomials(const FieldT &t)
{
  /* Compute Lagrange polynomial of size m, with m+1 points (x_0, y_0), ... ,(x_m, y_m) */
  /* Evaluate for x = t */
  /* Return coeffs for each l_j(x) = (l / l_i[j]) * w[j] */

  /* for all i: w[i] = (1 / r) * w[i-1] * (1 - a[i]^m-i+1) / (1 - a[i]^-i) */

  /**
   * If t equals one of the geometric progression values,
   * then output 1 at the right place, and 0 elsewhere.
   */
  for (size_t i = 0; i < this->m; ++i)
  {
    if (this->geometric_sequence[i] == t) // i.e., t equals a[i]
    {
      std::vector<FieldT> res(this->m, FieldT::zero());
      res[i] = FieldT::one();
      return res;
    }
  }

  /**
   * Otherwise, if t does not equal any of the geometric progression values,
   * then compute each Lagrange coefficient.
   */
  std::vector<FieldT> l(this->m);
  l[0] = t - this->geometric_sequence[0];

  std::vector<FieldT> g(this->m);
  g[0] = FieldT::zero();

  FieldT l_vanish = l[0];
  FieldT g_vanish = FieldT::one();
  for (size_t i = 1; i < this->m; i++)
  {
    l[i] = t - this->geometric_sequence[i];
    g[i] = FieldT::one() - this->geometric_sequence[i];

    l_vanish *= l[i];
    g_vanish *= g[i];
  }

  FieldT r = this->geometric_sequence[this->m-1].inverse();
  FieldT r_i = r;

  std::vector<FieldT> g_i(this->m);
  g_i[0] = g_vanish.inverse();

  l[0] = l_vanish * l[0].inverse() * g_i[0];
  for (size_t i = 1; i < this->m; i++)
  {
    g_i[i] = g_i[i-1] * g[this->m-i] * -g[i].inverse() * this->geometric_sequence[i];
    l[i] = l_vanish * r_i * l[i].inverse() * g_i[i];
    r_i *= r;
  }

  return l;
}

template<typename FieldT>
FieldT geometric_sequence_domain<FieldT>::get_domain_element(const size_t idx)
{
  return this->geometric_sequence[idx];
}

template<typename FieldT>
FieldT geometric_sequence_domain<FieldT>::compute_vanishing_polynomial(const FieldT &t)
{
  /* Notes: Z = prod_{i = 0 to m} (t - a[i]) */
  /* Better approach: Montgomery Trick + Divide&Conquer/FFT */
  FieldT Z = FieldT::one();
  for (size_t i = 0; i < this->m; i++)
  {
    Z *= (t - this->geometric_sequence[i]);
  }
  return Z;
}

template<typename FieldT>
void geometric_sequence_domain<FieldT>::add_poly_Z(const FieldT &coeff, std::vector<FieldT> &H)
{
  if (H.size() != this->m+1) throw DomainSizeException("geometric: expected H.size() == this->m+1");

  if (coeff == FieldT::zero()) {
      H.resize(this->m+1, FieldT::zero());
      return;
  }

#ifdef MULTICORE
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < this->m+1; i++)
  {
    H[i] += (this->vanishing_poly[i] * coeff);
  }
}

template<typename FieldT>
void geometric_sequence_domain<FieldT>::divide_by_Z_on_coset(std::vector<FieldT> &P)
{
  for (size_t i = 0; i < this->m; ++i)
  {
    P[i] *= this->vanishing_poly_eval_coset_inverse[i];
  }
}

template<typename FieldT>
FieldT geometric_sequence_domain<FieldT>::get_coset()
{
  return coset;
}

template<typename FieldT>
void geometric_sequence_domain<FieldT>::do_precomputation(const size_t m, const FieldT &generator)
{
  if (generator == FieldT::zero())
    throw InvalidSizeException("geometric(): non-zero element required to initialize geometric sequence");

  this->geometric_sequence = std::vector<FieldT>(m + 1, FieldT::zero());
  this->geometric_sequence[0] = FieldT::one();

  this->geometric_triangular_sequence = std::vector<FieldT>(m + 1, FieldT::zero());
  this->geometric_triangular_sequence[0] = FieldT::one();

  // Compute domain elements
  for (size_t i = 1; i < m + 1; i++) {
    this->geometric_sequence[i] = this->geometric_sequence[i-1] * generator;
    this->geometric_triangular_sequence[i] = this->geometric_triangular_sequence[i-1] * this->geometric_sequence[i-1];
  }

  // Coset multiplier
  this->coset = this->geometric_sequence[m];

  // Compute vanishing polynomial
  this->vanishing_poly = std::vector<FieldT>(m + 1, FieldT::zero());

  this->vanishing_poly[m] = FieldT::one();
  for (size_t i = 0; i < m; i++) {
    this->vanishing_poly[m] *= (this->geometric_sequence[m] - this->geometric_sequence[i]);
  }
  iFFT(this->vanishing_poly);

  // Evaluate vanishing polynomial on coset
  this->vanishing_poly_eval_coset_inverse = this->vanishing_poly;
  cosetFFT(this->vanishing_poly_eval_coset_inverse, this->coset);

  for (size_t i = 0; i < m; i++) {
    this->vanishing_poly_eval_coset_inverse[i] = this->vanishing_poly_eval_coset_inverse[i].inverse();
  }

  this->vanishing_poly_eval_coset_inverse.resize(m);

  this->m = m;

  this->geometric_sequence.resize(m);
  this->geometric_triangular_sequence.resize(m);
}

} // libfqfft

#endif // GEOMETRIC_SEQUENCE_DOMAIN_TCC_
