#include <cassert>
#include <cstdint>
#include "mpn.h"

typedef uint64_t mpn_double_digit;
static_assert(sizeof(mpn_double_digit) == 2 * sizeof(mpn_digit), "size alignment");

int mpn_manager::compare(mpn_digit const * a, unsigned lnga,
                         mpn_digit const * b, unsigned lngb) const {
    int res = 0;

    unsigned j = std::max(lnga, lngb);
    for (; j-- > 0 && res == 0;) {
        mpn_digit u_j = (j < lnga) ? a[j] : 0;
        mpn_digit v_j = (j < lngb) ? b[j] : 0;
        if (u_j > v_j)
            res = 1;
        else if (u_j < v_j)
            res = -1;
    }

    return res;
}

bool mpn_manager::add(mpn_digit const * a, unsigned lnga,
                      mpn_digit const * b, unsigned lngb,
                      mpn_digit * c, unsigned lngc_alloc,
                      unsigned * plngc) const {
    // Essentially Knuth's Algorithm A
    unsigned len = std::max(lnga, lngb);
    assert(lngc_alloc == len+1 && len > 0);
    mpn_digit k = 0;
    mpn_digit r;
    bool c1, c2;
    for (unsigned j = 0; j < len; ++j) {
        mpn_digit u_j = (j < lnga) ? a[j] : 0;
        mpn_digit v_j = (j < lngb) ? b[j] : 0;
        r = u_j + v_j; c1 = r < u_j;
        c[j] = r + k;  c2 = c[j] < r;
        k = c1 | c2;
    }
    c[len] = k;
    unsigned &os = *plngc;
    for (os = len+1; os > 1 && c[os-1] == 0; ) os--;
    assert(os > 0 && os <= len+1);
    return true; // return k != 0?
}

bool mpn_manager::sub(mpn_digit const * a, unsigned lnga,
                      mpn_digit const * b, unsigned lngb,
                      mpn_digit * c, mpn_digit * pborrow) const {
    // Essentially Knuth's Algorithm S
    unsigned len = std::max(lnga, lngb);
    mpn_digit & k = *pborrow; k = 0;
    mpn_digit r;
    bool c1, c2;
    for (unsigned j = 0; j < len; ++j) {
        mpn_digit u_j = (j < lnga) ? a[j] : 0;
        mpn_digit v_j = (j < lngb) ? b[j] : 0;
        r = u_j - v_j; c1 = r > u_j;
        c[j] = r - k;  c2 = c[j] > r;
        k = c1 | c2;
    }
    return true; // return k != 0?
}

bool mpn_manager::mul(mpn_digit const * a, unsigned lnga,
                      mpn_digit const * b, unsigned lngb,
                      mpn_digit * c) const {
    // Essentially Knuth's Algorithm M.
    // Perhaps implement a more efficient version, see e.g., Knuth, Section 4.3.3.
    unsigned i;
    mpn_digit k;

#define DIGIT_BITS (sizeof(mpn_digit)*8)
#define HALF_BITS (sizeof(mpn_digit)*4)

    for (unsigned i = 0; i < lnga; ++i)
        c[i] = 0;

    for (unsigned j = 0; j < lngb; ++j) {
        mpn_digit v_j = b[j];
        if (v_j == 0) { // This branch may be omitted according to Knuth.
            c[j+lnga] = 0;
        }
        else {
            k = 0;
            for (i = 0; i < lnga; ++i) {
                mpn_digit u_i = a[i];
                mpn_double_digit t;
                t = ((mpn_double_digit)u_i * (mpn_double_digit)v_j) +
                    (mpn_double_digit) c[i+j] +
                    (mpn_double_digit) k;

                c[i+j] = (t << DIGIT_BITS) >> DIGIT_BITS;
                k = t >> DIGIT_BITS;
            }
            c[j+lnga] = k;
        }
    }

    return true;
}

#define MASK_FIRST (~((mpn_digit)(-1) >> 1))
#define FIRST_BITS(N, X) ((X) >> (DIGIT_BITS-(N)))
#define LAST_BITS(N, X) (((X) << (DIGIT_BITS-(N))) >> (DIGIT_BITS-(N)))
#define BASE ((mpn_double_digit)0x01 << DIGIT_BITS)

bool mpn_manager::div(mpn_digit const * numer, unsigned lnum,
                      mpn_digit const * denom, unsigned lden,
                      mpn_digit * quot,
                      mpn_digit * rem) {
    // Проверка деления на ноль
    assert(lden > 0 && "Division by zero");
    assert(lnum > 0 && "Empty numerator");
    assert(denom != nullptr && numer != nullptr);
    assert(quot != nullptr && rem != nullptr);
    assert(denom[lden - 1] != 0 && "Leading digit of denominator must be non-zero");

    bool res = false;
    if (lnum == 1 && lden == 1) {
        // Однозначное деление — быстрый путь
        *quot = numer[0] / denom[0];
        *rem  = numer[0] % denom[0];
    }
    else if (lnum < lden ||
             (lnum == lden && numer[lnum - 1] < denom[lden - 1])) {
        // Числитель меньше знаменателя: частное = 0, остаток = числитель.
        // (включает случай lnum < lden — БЕЗ unsigned underflow)
        *quot = 0;
        for (unsigned i = 0; i < lden; ++i)
            rem[i] = (i < lnum) ? numer[i] : 0;
    }
    else {
        // Нормальный путь — алгоритм Кнута D
        mpn_sbuffer u, v, t_ms, t_ab;
        unsigned d = div_normalize(numer, lnum, denom, lden, u, v);
        if (lden == 1)
            res = div_1(u, v[0], quot);
        else
            res = div_n(u, v, quot, rem, t_ms, t_ab);
        div_unnormalize(u, v, d, rem);
    }
    return res;
}

unsigned mpn_manager::div_normalize(mpn_digit const * numer, unsigned lnum,
                                  mpn_digit const * denom, unsigned lden,
                                  mpn_sbuffer & n_numer,
                                  mpn_sbuffer & n_denom) const
{
    unsigned d = 0;
    while (lden > 0 && ((denom[lden-1] << d) & MASK_FIRST) == 0) d++;
    assert(d < DIGIT_BITS);

    n_numer.resize(lnum+1);
    n_denom.resize(lden);

    if (d == 0) {
        n_numer[lnum] = 0;
        for (unsigned i = 0; i < lnum; ++i)
            n_numer[i] = numer[i];
        for (unsigned i = 0; i < lden; ++i)
            n_denom[i] = denom[i];
    }
    else if (lnum != 0) {
        assert(lden > 0);
        mpn_digit q = FIRST_BITS(d, numer[lnum-1]);
        n_numer[lnum] = q;
        for (unsigned i = lnum-1; i > 0; i--)
            n_numer[i] = (numer[i] << d) | FIRST_BITS(d, numer[i-1]);
        n_numer[0] = numer[0] << d;
        for (unsigned i = lden-1; i > 0; i--)
            n_denom[i] = denom[i] << d | FIRST_BITS(d, denom[i-1]);
        n_denom[0] = denom[0] << d;
    }
    else {
        d = 0;
    }

    return d;
}

void mpn_manager::div_unnormalize(mpn_sbuffer & numer, mpn_sbuffer & denom,
                                  unsigned d, mpn_digit * rem) const {
    if (d == 0) {
        for (unsigned i = 0; i < denom.size(); ++i)
            rem[i] = numer[i];
    }
    else {
        for (unsigned i = 0; i < denom.size()-1; ++i)
            rem[i] = numer[i] >> d | (LAST_BITS(d, numer[i+1]) << (DIGIT_BITS-d));
        rem[denom.size()-1] = numer[denom.size()-1] >> d;
    }
}

bool mpn_manager::div_1(mpn_sbuffer & numer, mpn_digit const denom,
                        mpn_digit * quot) const {
    mpn_double_digit q_hat, temp, ms;
    mpn_digit borrow;

    for (unsigned j = numer.size()-1; j > 0; j--) {
        temp = (((mpn_double_digit)numer[j]) << DIGIT_BITS) | ((mpn_double_digit)numer[j-1]);
        q_hat = temp / (mpn_double_digit) denom;
        assert(q_hat < BASE);
        ms = temp - (q_hat * (mpn_double_digit) denom);
        borrow = ms > temp;
        numer[j-1] = (mpn_digit) ms;
        numer[j] = ms >> DIGIT_BITS;
        quot[j-1] = (mpn_digit) q_hat;
        if (borrow) {
            quot[j-1]--;
            numer[j] = numer[j-1] + denom;
        }
    }

    return true; // return rem != 0?
}

bool mpn_manager::div_n(mpn_sbuffer & numer, mpn_sbuffer const & denom,
                        mpn_digit * quot, mpn_digit * rem,
                        mpn_sbuffer & ms, mpn_sbuffer & ab) const {
    assert(denom.size() > 1);

    // This is essentially Knuth's Algorithm D.
    unsigned m = numer.size() - denom.size();
    unsigned n = denom.size();

    assert(numer.size() == m+n);

    ms.resize(n+1);

    mpn_double_digit q_hat, temp, r_hat;
    mpn_digit borrow;

    for (unsigned j = m; j-- > 0; ) {
        temp = (((mpn_double_digit)numer[j+n]) << DIGIT_BITS) | ((mpn_double_digit)numer[j+n-1]);
        q_hat = temp / (mpn_double_digit) denom[n-1];
        r_hat = temp % (mpn_double_digit) denom[n-1];
        recheck:
        if (q_hat >= BASE ||
            ((q_hat * denom[n-2]) > ((r_hat << DIGIT_BITS) + numer[j+n-2]))) {
                q_hat--;
                r_hat += denom[n-1];
                if (r_hat < BASE) goto recheck;
        }
        assert(q_hat < BASE);
        // Replace numer[j+n]...numer[j] with
        // numer[j+n]...numer[j] - q * (denom[n-1]...denom[0])
        mpn_digit q_hat_small = (mpn_digit)q_hat;
        mul(&q_hat_small, 1, denom.data(), n, ms.data());
        sub(&numer[j], n+1, ms.data(), n+1, &numer[j], &borrow);
        quot[j] = q_hat_small;
        if (borrow) {
            quot[j]--;
            ab.resize(n+2);
            unsigned real_size;
            add(denom.data(), n, &numer[j], n+1, ab.data(), n+2, &real_size);
            for (unsigned i = 0; i < n+1; ++i)
                numer[j+i] = ab[i];
        }
    }

    return true; // return rem != 0?
}

char * mpn_manager::to_string(mpn_digit const * a, unsigned lng, char * buf, unsigned lbuf) const {
    assert(buf && lbuf > 0);

    if (lng == 1) {
#ifdef _WINDOWS
        sprintf_s(buf, lbuf, "%u", *a);
#else
        snprintf(buf, lbuf, "%u", *a);
#endif
    }
    else {
        mpn_sbuffer temp(lng, 0), t_numer(lng+1, 0), t_denom(1, 0);
        for (unsigned i = 0; i < lng; ++i)
            temp[i] = a[i];

        unsigned j = 0;
        mpn_digit rem;
        mpn_digit ten = 10;
        while (!temp.empty() && (temp.size() > 1 || temp[0] != 0)) {
            unsigned d = div_normalize(&temp[0], temp.size(), &ten, 1, t_numer, t_denom);
            div_1(t_numer, t_denom[0], &temp[0]);
            div_unnormalize(t_numer, t_denom, d, &rem);
            buf[j++] = '0' + rem;
            while (!temp.empty() && temp.back() == 0)
                temp.pop_back();
        }
        buf[j] = 0;

        j--;
        unsigned mid = (j/2) + ((j % 2) ? 1 : 0);
        for (unsigned i = 0; i < mid; ++i)
            std::swap(buf[i], buf[j-i]);
    }
    return buf;
}
