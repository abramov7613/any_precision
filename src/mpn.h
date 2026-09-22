#pragma once

#include <ostream>
#include <vector>
#include <string>

typedef unsigned int mpn_digit;

class mpn_manager {

public:
    int compare(mpn_digit const * a, unsigned lnga,
                mpn_digit const * b, unsigned lngb) const;

    void add(mpn_digit const * a, unsigned lnga,
             mpn_digit const * b, unsigned lngb,
             mpn_digit *c, unsigned lngc_alloc,
             unsigned * plngc) const;

    void sub(mpn_digit const * a, unsigned lnga,
             mpn_digit const * b, unsigned lngb,
             mpn_digit * c, mpn_digit * pborrow) const;

    void mul(mpn_digit const * a, unsigned lnga,
             mpn_digit const * b, unsigned lngb,
             mpn_digit * c) const;

    void div(mpn_digit const * numer, unsigned lnum,
             mpn_digit const * denom, unsigned lden,
             mpn_digit * quot,
             mpn_digit * rem);

    char * to_string(mpn_digit const * a, unsigned lng,
                     char * buf, unsigned lbuf) const;

    std::string to_string(mpn_digit const * a, unsigned lng) const;

private:
    using mpn_sbuffer = std::vector<mpn_digit>;

    unsigned div_normalize(mpn_digit const * numer, unsigned lnum,
                         mpn_digit const * denom, unsigned lden,
                         mpn_sbuffer & n_numer,
                         mpn_sbuffer & n_denom) const;

    void div_unnormalize(mpn_sbuffer & numer, mpn_sbuffer & denom,
                         unsigned d, mpn_digit * rem) const;

    void div_1(mpn_sbuffer & numer, mpn_digit denom,
               mpn_digit * quot) const;

    void div_n(mpn_sbuffer & numer, mpn_sbuffer const & denom,
               mpn_digit * quot, mpn_digit * rem,
               mpn_sbuffer & ms, mpn_sbuffer & ab) const;
};
