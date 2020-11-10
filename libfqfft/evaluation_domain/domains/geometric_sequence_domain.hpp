/** @file
 *****************************************************************************

 Declaration of interfaces for the "geometric sequence" evaluation domain.

 These functions use a geometric sequence of size m to perform evaluation.

 *****************************************************************************
 * @author     This file is part of libfqfft, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef GEOMETRIC_SEQUENCE_DOMAIN_HPP
#define GEOMETRIC_SEQUENCE_DOMAIN_HPP

#include <libfqfft/evaluation_domain/evaluation_domain.hpp>

namespace libfqfft {

  template<typename FieldT>
  class geometric_sequence_domain : public evaluation_domain<FieldT> {
  public:

    std::vector<FieldT> geometric_sequence;
    std::vector<FieldT> geometric_triangular_sequence;
    std::vector<FieldT> vanishing_poly;
    std::vector<FieldT> vanishing_poly_eval_coset_inverse;
    FieldT coset;

    geometric_sequence_domain(size_t m);
    geometric_sequence_domain(size_t m, const FieldT &generator);

    void FFT(std::vector<FieldT> &a);
    void iFFT(std::vector<FieldT> &a);
    void cosetFFT(std::vector<FieldT> &a, const FieldT &s);
    void icosetFFT(std::vector<FieldT> &a, const FieldT &s);
    std::vector<FieldT> evaluate_all_lagrange_polynomials(const FieldT &t);
    FieldT get_domain_element(size_t idx);
    FieldT compute_vanishing_polynomial(const FieldT &t);
    void add_poly_Z(const FieldT &coeff, std::vector<FieldT> &H);
    void divide_by_Z_on_coset(std::vector<FieldT> &P);
    FieldT get_coset() override;

    void do_precomputation(size_t m, const FieldT &generator);
  };

} // libfqfft

#include <libfqfft/evaluation_domain/domains/geometric_sequence_domain.tcc>

#endif // GEOMETRIC_SEQUENCE_DOMAIN_HPP
