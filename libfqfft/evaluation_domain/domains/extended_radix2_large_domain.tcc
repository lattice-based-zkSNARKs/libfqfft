/** @file
 *****************************************************************************

 Implementation of interfaces for the "extended radix-2 large" evaluation domain.

 See extended_radix2_large_domain.hpp .

 *****************************************************************************/

#ifndef EXTENDED_RADIX2_LARGE_DOMAIN_TCC_
#define EXTENDED_RADIX2_LARGE_DOMAIN_TCC_

#include <libfqfft/evaluation_domain/domains/basic_radix2_domain_aux.hpp>

namespace libfqfft {

template<typename FieldT>
extended_radix2_large_domain<FieldT>::extended_radix2_large_domain(const size_t m) : evaluation_domain<FieldT>(m) {
    if (m <= 1) throw DomainSizeException("extended_radix2(): expected m > 1");
    
    nroots = 1 << FieldT::s;
    // Important: field must have size at least 2*nroots*ncosets (since we are using g^2 H, g^4 H, ... as the cosets
    //            of the subgroup of 2^s-roots of unity H. This is to additionally support evaluation on the cosets
    //            g H, g^3 H, ...)
    if (m % nroots != 0) throw DomainSizeException("extended_radix2_large(): expected m == 2^s * n_cosets");

    ncosets = m / nroots;
    if (ncosets == 1) throw DomainSizeException("extended_radix2_large(): expected more than 1 coset");

    shift = FieldT::multiplicative_generator * FieldT::multiplicative_generator;
    
    try { omega = libff::get_root_of_unity<FieldT>(nroots); }
    catch (const std::invalid_argument& e) { throw DomainSizeException(e.what()); }

    geo_domain = new geometric_sequence_domain<FieldT>(ncosets, shift^nroots);

    FieldT omega_i = FieldT::one();
    domain_elements.resize(m);
    for (size_t i = 0; i < nroots; i++) {
        domain_elements[i] = omega_i;
        for (size_t j = 1; j < ncosets; j++) {
            domain_elements[j*nroots + i] = domain_elements[(j - 1)*nroots + i] * shift;
        }
        omega_i *= omega;
    }

    vanishing_polynomial.resize(ncosets + 1);
    vanishing_polynomial[0] = -FieldT::one();
    vanishing_polynomial[1] = FieldT::one();

    for (size_t j = 1; j < ncosets; j++) {
        const FieldT mult = -geo_domain->get_domain_element(j);
        for (size_t k = j + 1; k >= 1; k--) {
            vanishing_polynomial[k] *= mult;
            vanishing_polynomial[k] += vanishing_polynomial[(k-1)];
        }
        vanishing_polynomial[0] *= mult;
    }

    vanishing_coset_invs.resize(ncosets);
    FieldT shift_j = FieldT::one();
    for (size_t j = 0; j < ncosets; j++) {
        vanishing_coset_invs[j] = compute_vanishing_polynomial(this->get_coset() * shift_j).inverse();
        shift_j *= shift;
    }
}

template<typename FieldT>
void extended_radix2_large_domain<FieldT>::FFT(std::vector<FieldT> &a) {
    if (a.size() != this->m) throw DomainSizeException("extended_radix2: expected a.size() == this->m");

    std::vector<std::vector<FieldT>> b(ncosets);
    for (size_t j = 0; j < ncosets; j++) {
        b[j].resize(nroots, FieldT::zero());
    }

    std::vector<FieldT> scale(ncosets, FieldT::one());
    for (size_t j = 1; j < ncosets; j++) {
        scale[j] = shift * scale[j - 1];
    }
    std::vector<FieldT> scale_j(ncosets, FieldT::one());

    for (size_t i = 0; i < nroots; i++) {
        std::vector<FieldT> p(ncosets);
        for (size_t j = 0; j < ncosets; j++) {
            p[j] = a[j*nroots + i];
        }

        this->geo_domain->FFT(p);

        for (size_t j = 0; j < ncosets; j++) {
            b[j][i] = p[j] * scale_j[j];
            scale_j[j] *= scale[j];
        }
    }

    for (size_t j = 0; j < ncosets; j++) {
        _basic_radix2_FFT(b[j], omega);
    }

    for (size_t i = 0; i < nroots; i++) {
        for (size_t j = 0; j < ncosets; j++) {
            a[j*nroots + i] = b[j][i];
        }
    }
}

template<typename FieldT>
void extended_radix2_large_domain<FieldT>::iFFT(std::vector<FieldT> &a) {
    if (a.size() != this->m) throw DomainSizeException("extended_radix2: expected a.size() == this->m");
    
    const FieldT omega_inverse = omega.inverse();

    std::vector<std::vector<FieldT>> b(ncosets);
    for (size_t j = 0; j < ncosets; j++) {
        b[j]  = std::vector<FieldT>(a.begin() + j*nroots, a.begin() + (j + 1)*nroots);
        _basic_radix2_FFT(b[j], omega_inverse);
    }

    const FieldT inv = FieldT(nroots).inverse();
    const FieldT shiftinv = shift.inverse();
    std::vector<FieldT> scale(ncosets, FieldT::one());
    for (size_t j = 1; j < ncosets; j++) {
        scale[j] = shiftinv * scale[j - 1];
    }
    std::vector<FieldT> scale_j(ncosets, inv);

    for (size_t i = 0; i < nroots; i++) {
        std::vector<FieldT> p(ncosets);
        for (size_t j = 0; j < ncosets; j++) {
            p[j] = b[j][i] * scale_j[j];
            scale_j[j] *= scale[j];
        }

        this->geo_domain->iFFT(p);

        for (size_t j = 0; j < ncosets; j++) {
            a[j*nroots + i] = p[j];
        }
    }
}

template<typename FieldT>
void extended_radix2_large_domain<FieldT>::cosetFFT(std::vector<FieldT> &a, const FieldT &g) {
    if (g != FieldT::multiplicative_generator) {
        throw "Coset multiplier should be FieldT::multiplicative_generator";
    }
    _multiply_by_coset(a, g);
    FFT(a);
}

template<typename FieldT>
void extended_radix2_large_domain<FieldT>::icosetFFT(std::vector<FieldT> &a, const FieldT &g) {
    if (g != FieldT::multiplicative_generator) {
        throw "Coset multiplier should be FieldT::multiplicative_generator";
    }
    iFFT(a);
    _multiply_by_coset(a, g.inverse());
}

template<typename FieldT>
std::vector<FieldT> extended_radix2_large_domain<FieldT>::evaluate_all_lagrange_polynomials(const FieldT &t) {
    std::vector<FieldT> scale(ncosets, FieldT::zero());
    const FieldT td = t^nroots;
    FieldT num = FieldT::one();
    for (size_t j = 0; j < ncosets; j++) {
         num *= (td - geo_domain->get_domain_element(j));
         
         if (num == FieldT::zero()) {
            scale[j] = FieldT::one();
            break;
         }
    }

    if (num != FieldT::zero()) {
        for (size_t j = 0; j < ncosets; j++) {
            scale[j] = num * (td - geo_domain->get_domain_element(j)).inverse();
            for (size_t k = 0; k < ncosets; k++) {
                if (k == j) {
                    continue;
                }
                scale[j] *= (geo_domain->get_domain_element(j) - geo_domain->get_domain_element(k)).inverse();
            }
        }
    }

    std::vector<FieldT> result(this->m);

    FieldT z = t;

    for(size_t j = 0; j < ncosets; j++) {
        const std::vector<FieldT> L = _basic_radix2_evaluate_all_lagrange_polynomials(nroots, z);
        for (size_t i = 0; i < nroots; i++) {
            result[j*nroots + i] = L[i] * scale[j];
        }
        z *= shift.inverse();
    }

    return result;
}

template<typename FieldT>
FieldT extended_radix2_large_domain<FieldT>::get_domain_element(const size_t idx) {
    return domain_elements[idx];
}

template<typename FieldT>
FieldT extended_radix2_large_domain<FieldT>::compute_vanishing_polynomial(const FieldT &t) {
    const FieldT td = t^nroots;
    FieldT val = FieldT::one();
    for (size_t j = 0; j < ncosets; j++) {
        val *= (td - geo_domain->get_domain_element(j));
    }
    return val;
}

template<typename FieldT>
void extended_radix2_large_domain<FieldT>::add_poly_Z(const FieldT &coeff, std::vector<FieldT> &H) {
    if (H.size() != this->m+1) throw DomainSizeException("extended_radix2_large: expected H.size() == this->m+1");

    for (size_t j = 0; j < ncosets + 1; j++) {
        H[j*nroots] += (vanishing_polynomial[j] * coeff);
    }
}

template<typename FieldT>
void extended_radix2_large_domain<FieldT>::divide_by_Z_on_coset(std::vector<FieldT> &P) {
    for (size_t i = 0; i < nroots; i++) {
        for (size_t j = 0; j < ncosets; j++) {
           P[j*nroots + i] *= vanishing_coset_invs[j];
        }
    }
}

} // libfqfft

#endif // EXTENDED_RADIX2_LARGE_DOMAIN_TCC_
