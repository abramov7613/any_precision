#include <cassert>
#include <cstdint>
#include <cstdio>
#include <algorithm>
#include "mpn.h"

typedef uint64_t mpn_double_digit;
static_assert(sizeof(mpn_double_digit) == 2 * sizeof(mpn_digit), "size alignment");

constexpr unsigned DIGIT_BITS = sizeof(mpn_digit) * 8;
constexpr mpn_double_digit BASE = (mpn_double_digit)0x01 << DIGIT_BITS;
constexpr mpn_digit MASK_FIRST = mpn_digit(1) << (DIGIT_BITS - 1);

constexpr mpn_digit FIRST_BITS(unsigned n, mpn_digit x) {
    if (n >= DIGIT_BITS) return x;
    return x >> (DIGIT_BITS - n);
}

constexpr mpn_digit LAST_BITS(unsigned n, mpn_digit x) {
    if (n >= DIGIT_BITS) return x;
    return x << (DIGIT_BITS - n) >> (DIGIT_BITS - n);
}


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

void mpn_manager::add(mpn_digit const * a, unsigned lnga,
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
}

void mpn_manager::sub(mpn_digit const * a, unsigned lnga,
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
}

void mpn_manager::mul(mpn_digit const * a, unsigned lnga,
                      mpn_digit const * b, unsigned lngb,
                      mpn_digit * c) const {
    // Essentially Knuth's Algorithm M.
    // Perhaps implement a more efficient version, see e.g., Knuth, Section 4.3.3.
    unsigned i;
    mpn_digit k;

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
}

void mpn_manager::div(mpn_digit const * numer, unsigned lnum,
                      mpn_digit const * denom, unsigned lden,
                      mpn_digit * quot,
                      mpn_digit * rem) {
    // Проверка деления на ноль
    assert(lden > 0 && "Division by zero");
    assert(lnum > 0 && "Empty numerator");
    assert(denom != nullptr && numer != nullptr);
    assert(quot != nullptr && rem != nullptr);
    assert(denom[lden - 1] != 0 && "Leading digit of denominator must be non-zero");

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
            div_1(u, v[0], quot);
        else
            div_n(u, v, quot, rem, t_ms, t_ab);
        div_unnormalize(u, v, d, rem);
    }
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
    const auto denom_size = denom.size();
    if (denom_size == 0) return; // ...
    if (d == 0) {
        for (unsigned i = 0; i < denom_size; ++i)
            rem[i] = numer[i];
    } else {
        const unsigned limit = denom_size - 1;
        for (unsigned i = 0; i < limit; ++i) {
            rem[i] = numer[i] >> d | (LAST_BITS(d, numer[i+1]) << (DIGIT_BITS-d));
        }
        rem[limit] = numer[limit] >> d;
    }
}

void mpn_manager::div_1(mpn_sbuffer & numer, mpn_digit const denom,
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
}

void mpn_manager::div_n(mpn_sbuffer & numer, mpn_sbuffer const & denom,
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
        while (q_hat >= BASE ||
              (q_hat * denom[n-2]) > ((r_hat << DIGIT_BITS) + numer[j+n-2])) {
            q_hat--;
            r_hat += denom[n-1];
            if (r_hat >= BASE) break;
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
}

std::string mpn_manager::to_string(mpn_digit const * a, unsigned lng) const {
    if (lng == 0) return "0";
    if (lng == 1) return std::to_string(a[0]);
    // Проверка на нулевое число
    bool all_zero = true;
    for (unsigned i = 0; i < lng; ++i) {
        if (a[i] != 0) { all_zero = false; break; }
    }
    if (all_zero) return "0";
    // Деление на 10^9 за итерацию — 9 десятичных цифр сразу.
    // 10^9 помещается в uint32_t (макс ~4.29 * 10^9).
    constexpr mpn_digit ten_pow9 = 1000000000u;
    mpn_sbuffer temp(a, a + lng);
    mpn_sbuffer t_numer, t_denom;
    mpn_digit rem;
    std::vector<mpn_digit> groups;
    groups.reserve(lng); // грубая оценка сверху
    while (!temp.empty() && (temp.size() > 1 || temp[0] != 0)) {
        unsigned d = div_normalize(temp.data(), static_cast<unsigned>(temp.size()),
                                   &ten_pow9, 1, t_numer, t_denom);
        div_1(t_numer, t_denom[0], temp.data());
        div_unnormalize(t_numer, t_denom, d, &rem);
        groups.push_back(rem);

        while (!temp.empty() && temp.back() == 0)
            temp.pop_back();
    }
    // Сборка строки: groups хранит группы от младшей к старшей,
    // поэтому идём с конца. Первая (старшая) группа — без дополнения нулями,
    // остальные добиваются до 9 цифр.
    std::string result;
    result.reserve(groups.size() * 9);
    result += std::to_string(groups.back());
    for (size_t i = groups.size() - 1; i > 0; --i) {
        char buf[10];
        std::snprintf(buf, sizeof(buf), "%09u", groups[i - 1]);
        result += buf;
    }
    return result;
}

char * mpn_manager::to_string(mpn_digit const * a, unsigned lng, char * buf, unsigned lbuf) const {
    assert(buf && lbuf > 0);
    const auto str = to_string(a, lng);
    const unsigned strsz = str.size();
    const auto& count = std::min(strsz, lbuf);
    std::copy_n(str.data(), count, buf);
    return buf;
}
