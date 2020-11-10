/** @file
 *****************************************************************************

 Declaration of interfaces for the "extended radix-2 large" evaluation domain.

 The domain has size m = 2^{k} * \ell and consists of
 "the m-th roots of unity" union "\ell cosets of these roots".

 *****************************************************************************/

#ifndef EXTENDED_RADIX2_LARGE_DOMAIN_HPP_
#define EXTENDED_RADIX2_LARGE_DOMAIN_HPP_

#include <libfqfft/evaluation_domain/evaluation_domain.hpp>
#include <libfqfft/evaluation_domain/domains/geometric_sequence_domain.hpp>

namespace libfqfft {

template<typename FieldT>
class extended_radix2_large_domain : public evaluation_domain<FieldT> {
public:

    size_t nroots;
    size_t ncosets;
    FieldT shift;
    FieldT omega;

    geometric_sequence_domain<FieldT> *geo_domain;
    std::vector<FieldT> domain_elements;
    std::vector<FieldT> vanishing_polynomial;
    std::vector<FieldT> vanishing_coset_invs;

    extended_radix2_large_domain(size_t m);

    void FFT(std::vector<FieldT> &a);
    void iFFT(std::vector<FieldT> &a);
    void cosetFFT(std::vector<FieldT> &a, const FieldT &g);
    void icosetFFT(std::vector<FieldT> &a, const FieldT &g);
    std::vector<FieldT> evaluate_all_lagrange_polynomials(const FieldT &t);
    FieldT get_domain_element(size_t idx);
    FieldT compute_vanishing_polynomial(const FieldT &t);
    void add_poly_Z(const FieldT &coeff, std::vector<FieldT> &H);
    void divide_by_Z_on_coset(std::vector<FieldT> &P);

};

} // libfqfft

#include <libfqfft/evaluation_domain/domains/extended_radix2_large_domain.tcc>

#endif // EXTENDED_RADIX2_LARGE_DOMAIN_HPP_

