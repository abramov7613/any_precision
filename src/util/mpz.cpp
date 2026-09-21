/*++
Copyright (c) 2006 Microsoft Corporation

Module Name:

    mpz.cpp

Abstract:

    <abstract>

Author:

    Leonardo de Moura (leonardo) 2010-06-17.

Revision History:

--*/
#include <cstring>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include "util/mpz.h"
#include "util/hash.h"
#include "util/bit_util.h"

constexpr auto QUOT_ONLY = 0;
constexpr auto REM_ONLY = 1;
constexpr auto QUOT_AND_REM = 2;

// Available GCD algorithms:
// #define EUCLID_GCD
// #define BINARY_GCD
// #define LS_BINARY_GCD
#define LEHMER_GCD

#ifdef __has_builtin
    #define HAS_BUILTIN(X) __has_builtin(X)
#else
    #define HAS_BUILTIN(X) 0
#endif
#if HAS_BUILTIN(__builtin_ctz)
#define _trailing_zeros32(X) __builtin_ctz(X)
#elif defined(_WINDOWS) && (defined(_M_X86) || (defined(_M_X64) && !defined(_M_ARM64EC))) && !defined(__clang__)
// This is needed for _tzcnt_u32 and friends.
#include <immintrin.h>
#define _trailing_zeros32(X) _tzcnt_u32(X)
#else
static uint32_t _trailing_zeros32(uint32_t x) {
    uint32_t r = 0;
    for (; 0 == (x & 1) && r < 32; ++r, x >>= 1);
    return r;
}
#endif

#if (defined(__LP64__) || defined(_WIN64)) && defined(_M_X64) && !defined(_M_ARM64EC)
#if HAS_BUILTIN(__builtin_ctzll)
#define _trailing_zeros64(X) __builtin_ctzll(X)
#elif !defined(__clang__)
#define _trailing_zeros64(X) _tzcnt_u64(X)
#endif
#else
static uint64_t _trailing_zeros64(uint64_t x) {
    uint64_t r = 0;
    for (; 0 == (x & 1) && r < 64; ++r, x >>= 1);
    return r;
}
#endif

#undef HAS_BUILTIN

unsigned trailing_zeros(uint32_t x) {
    return static_cast<unsigned>(_trailing_zeros32(x));
}

unsigned trailing_zeros(uint64_t x) {
    return static_cast<unsigned>(_trailing_zeros64(x));
}

#define _bit_min(x, y) (y + ((x - y) & ((int)(x - y) >> 31)))
#define _bit_max(x, y) (x - ((x - y) & ((int)(x - y) >> 31)))

template<typename Manager>
class _scoped_numeral {
public:
    typedef typename Manager::numeral numeral;
private:
    Manager & m_manager;
    numeral   m_num;
public:
    _scoped_numeral(Manager & m):m_manager(m) {}
    _scoped_numeral(_scoped_numeral const& n) :m_manager(n.m_manager) { m().set(m_num, n.m_num); }
    _scoped_numeral(_scoped_numeral && n) noexcept: m_manager(n.m_manager) { m().swap(m_num, n.m_num); }
    ~_scoped_numeral() { m_manager.del(m_num); }

    Manager & m() const { return m_manager; }

    operator numeral const &() const { return m_num; }
    operator numeral&() { return m_num; }
    numeral const & get() const { return m_num; }
    numeral & get() { return m_num; }

    _scoped_numeral & operator=(_scoped_numeral const & n) {
        if (this == &n)
            return *this;
        m().set(m_num, n.m_num);
        return *this;
    }

    _scoped_numeral & operator=(int n) {
        m().set(m_num, n);
        return *this;
    }

    _scoped_numeral & operator=(numeral const & n) {
        m().set(m_num, n);
        return *this;
    }

    void reset() {
        m().reset(m_num);
    }

    void swap(_scoped_numeral & n) noexcept {
        m().swap(m_num, n.m_num);
    }

    void swap(numeral & n) noexcept {
        m().swap(m_num, n);
    }

    _scoped_numeral & operator+=(numeral const & a) {
        m().add(m_num, a, m_num);
        return *this;
    }

    _scoped_numeral & operator-=(numeral const & a) {
        m().sub(m_num, a, m_num);
        return *this;
    }

    _scoped_numeral & operator*=(numeral const & a) {
        m().mul(m_num, a, m_num);
        return *this;
    }

    _scoped_numeral & operator/=(numeral const & a) {
        m().div(m_num, a, m_num);
        return *this;
    }

    _scoped_numeral & operator%=(numeral const & a) {
        m().rem(m_num, a, m_num);
        return *this;
    }

    friend bool operator==(_scoped_numeral const & a, numeral const & b) {
        return a.m().eq(a, b);
    }

    friend bool operator==(_scoped_numeral const & a, _scoped_numeral const & b) {
        return a.m().eq(a.m_num, b.m_num);
    }

    friend bool operator!=(_scoped_numeral const & a, numeral const & b) {
        return !a.m().eq(a, b);
    }

    friend bool operator!=(_scoped_numeral const & a, _scoped_numeral const & b) {
        return !(a == b);
    }


    friend bool operator<(_scoped_numeral const & a, numeral const & b) {
        return a.m().lt(a, b);
    }

    friend bool operator>(_scoped_numeral const & a, numeral const & b) {
        return a.m().gt(a, b);
    }

    friend bool operator<=(_scoped_numeral const & a, numeral const & b) {
        return a.m().le(a, b);
    }

    friend bool operator>=(_scoped_numeral const & a, numeral const & b) {
        return a.m().ge(a, b);
    }

    bool is_zero() const {
        return m().is_zero(*this);
    }

    bool is_pos() const {
        return m().is_pos(*this);
    }

    bool is_neg() const {
        return m().is_neg(*this);
    }

    bool is_nonpos() const {
        return m().is_nonpos(*this);
    }

    bool is_nonneg() const {
        return m().is_nonneg(*this);
    }

    friend bool is_zero(_scoped_numeral const & a) {
        return a.m().is_zero(a);
    }

    friend bool is_pos(_scoped_numeral const & a) {
        return a.m().is_pos(a);
    }

    friend bool is_neg(_scoped_numeral const & a) {
        return a.m().is_neg(a);
    }

    friend bool is_nonneg(_scoped_numeral const & a) {
        return a.m().is_nonneg(a);
    }

    friend bool is_nonpos(_scoped_numeral const & a) {
        return a.m().is_nonpos(a);
    }

    friend _scoped_numeral abs(_scoped_numeral const& a) {
        _scoped_numeral res(a);
        a.m().abs(res);
        return res;
    }

    void neg() {
        m().neg(m_num);
    }

    friend _scoped_numeral operator+(_scoped_numeral const & r1, numeral const & r2) {
        return _scoped_numeral(r1) += r2;
    }

    friend _scoped_numeral operator-(_scoped_numeral const & r1, numeral const & r2) {
        return _scoped_numeral(r1) -= r2;
    }

    friend _scoped_numeral operator*(_scoped_numeral const & r1, numeral const & r2) {
        return _scoped_numeral(r1) *= r2;
    }

    friend _scoped_numeral operator/(_scoped_numeral const & r1, numeral const & r2) {
        return _scoped_numeral(r1) /= r2;
    }

    friend std::ostream & operator<<(std::ostream & out, _scoped_numeral const & s) {
        s.m().display(out, s);
        return out;
    }
}; // class _scoped_numeral

unsigned u_gcd(unsigned u, unsigned v) {
    if (u == 0) return v;
    if (v == 0) return u;
    unsigned shift = _trailing_zeros32(u | v);
    u >>= _trailing_zeros32(u);
    if (u == 1 || v == 1) return 1 << shift;
    if (u == v) return u << shift;
    do {
        v >>= _trailing_zeros32(v);
        unsigned diff = u - v;
        unsigned mdiff = diff & (unsigned)((int)diff >> 31);
        u = v + mdiff; // min
        v = diff - 2 * mdiff;   // if v <= u: u - v, if v > u: v - u = u - v - 2 * (u - v)
    }
    while (v != 0);
    return u << shift;
}

uint64_t u64_gcd(uint64_t u, uint64_t v) {
    if (u == 0) return v;
    if (v == 0) return u;
    if (u == 1 || v == 1) return 1;
    auto shift = _trailing_zeros64(u | v);
    u >>= _trailing_zeros64(u);
    do {
        v >>= _trailing_zeros64(v);
        if (u > v) std::swap(u, v);
        v -= u;
    }
    while (v != 0);
    return u << shift;
}



template<bool SYNCH>
mpz_manager<SYNCH>::mpz_manager():
    m_allocator("mpz_manager") {
    set(m_int_min, -static_cast<int64_t>(INT_MIN));
    mpz one(1);
    set(m_two64, (uint64_t)UINT64_MAX);
    add(m_two64, one, m_two64);
}

template<bool SYNCH>
mpz_manager<SYNCH>::~mpz_manager() {
    del(m_two64);
    del(m_int_min);
}

template<bool SYNCH>
mpz_cell * mpz_manager<SYNCH>::allocate(unsigned capacity) {
    assert(capacity >= m_init_cell_capacity);
    mpz_cell * cell;
    if (SYNCH) {
        cell = reinterpret_cast<mpz_cell*>(memory::allocate(cell_size(capacity)));
    }
    else {
        cell = reinterpret_cast<mpz_cell*>(m_allocator.allocate(cell_size(capacity)));
    }
    cell->m_capacity = capacity;

    return cell;
}

template<bool SYNCH>
void mpz_manager<SYNCH>::deallocate(bool is_heap, mpz_cell * ptr) {
    if (is_heap) {
        if (SYNCH) {
            memory::deallocate(ptr);
        }
        else {
            m_allocator.deallocate(cell_size(ptr->m_capacity), ptr);
        }
    }
}

template<bool SYNCH>
mpz_manager<SYNCH>::sign_cell::sign_cell(mpz_manager& m, mpz const& a):
    m_local(reinterpret_cast<mpz_cell*>(m_bytes)), m_a(a) {
    m_local.m_ptr->m_capacity = capacity;
    m.get_sign_cell(a, m_sign, m_cell, m_local.m_ptr);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::del(mpz_manager<SYNCH>* m, mpz & a) {
    if (a.m_ptr) {
        assert(m);
        m->deallocate(a.m_owner == mpz_self, a.m_ptr);
        a.m_ptr = nullptr;
        a.m_kind = mpz_small;
        a.m_owner = mpz_self;
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::add(mpz const & a, mpz const & b, mpz & c) {
    if (is_small(a) && is_small(b)) {
        set_i64(c, i64(a) + i64(b));
    }
    else {
        big_add(a, b, c);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::sub(mpz const & a, mpz const & b, mpz & c) {
    if (is_small(a) && is_small(b)) {
        set_i64(c, i64(a) - i64(b));
    }
    else {
        big_sub(a, b, c);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::set_big_i64(mpz & c, int64_t v) {
    if (c.m_ptr == nullptr) {
        c.m_ptr = allocate(m_init_cell_capacity);
        c.m_owner = mpz_self;
    }
    c.m_kind = mpz_large;
    assert(capacity(c) >= m_init_cell_capacity);
    uint64_t _v;
    if (v == std::numeric_limits<int64_t>::min()) {
        // min-int is even
        _v = -(v/2);
        c.m_val = -1;
    }
    else if (v < 0) {
        _v = -v;
        c.m_val = -1;
    }
    else {
        _v = v;
        c.m_val = 1;
    }
    if (sizeof(digit_t) == sizeof(uint64_t)) {
        // 64-bit machine
        digits(c)[0] = static_cast<digit_t>(_v);
        c.m_ptr->m_size = 1;
    }
    else {
        // 32-bit machine
        digits(c)[0] = static_cast<unsigned>(_v);
        digits(c)[1] = static_cast<unsigned>(_v >> 32);
        c.m_ptr->m_size = digits(c)[1] == 0 ? 1 : 2;
    }
    if (v == std::numeric_limits<int64_t>::min()) {
        big_add(c, c, c);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::set_big_ui64(mpz & c, uint64_t v) {
    if (c.m_ptr == nullptr) {
        c.m_ptr = allocate(m_init_cell_capacity);
        c.m_owner = mpz_self;
    }
    c.m_kind = mpz_large;
    assert(capacity(c) >= m_init_cell_capacity);
    c.m_val = 1;
    if (sizeof(digit_t) == sizeof(uint64_t)) {
        // 64-bit machine
        digits(c)[0] = static_cast<digit_t>(v);
        c.m_ptr->m_size = 1;
    }
    else {
        // 32-bit machine
        digits(c)[0] = static_cast<unsigned>(v);
        digits(c)[1] = static_cast<unsigned>(v >> 32);
        c.m_ptr->m_size = digits(c)[1] == 0 ? 1 : 2;
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::set(mpz_cell& src, mpz & a, int sign, unsigned sz) {
    unsigned i = sz;
    for (; i > 0 && src.m_digits[i-1] == 0; --i) ;

    if (i == 0) {
        // src is zero
        set(a, 0);
        return;
    }

    unsigned d = src.m_digits[0];
    if (i == 1 && d <= INT_MAX) {
        // src fits is a fixnum
        a.m_val = sign < 0 ? -static_cast<int>(d) : static_cast<int>(d);
        a.m_kind = mpz_small;
        return;
    }

    set_digits(a, i, src.m_digits);
    a.m_val = sign;

    assert(a.m_kind == mpz_large);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::set(mpz & a, char const * val) {
    set(a, 0);
    mpz ten(10);
    mpz tmp;
    char const * str = val;
    bool sign = false;
    while (str[0] == ' ') ++str;
    if (str[0] == '-')
        sign = true;
    while (str[0]) {
        if ('0' <= str[0] && str[0] <= '9') {
            assert(str[0] - '0' <= 9);
            mul(a, ten, tmp);
            add(tmp, mk_z(str[0] - '0'), a);
        }
        ++str;
    }
    del(tmp);
    if (sign)
        neg(a);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::set_digits(mpz & target, unsigned sz, digit_t const * digits) {
    // remove zero digits
    while (sz > 0 && digits[sz - 1] == 0)
        sz--;
    if (sz == 0)
        set(target, 0);
    else if (sz == 1)
        set(target, digits[0]);
    else {
        target.m_val = 1; // number is positive.
        if (target.m_ptr == nullptr) {
            unsigned c = sz < m_init_cell_capacity ? m_init_cell_capacity : sz;
            target.m_ptr             = allocate(c);
            target.m_ptr->m_size     = sz;
            target.m_ptr->m_capacity = c;
            target.m_kind            = mpz_large;
            target.m_owner           = mpz_self;
            memcpy(target.m_ptr->m_digits, digits, sizeof(digit_t) * sz);
        }
        else if (capacity(target) < sz) {
            assert(sz > m_init_cell_capacity);
            mpz_cell* ptr = allocate(sz);
            memcpy(ptr->m_digits, digits, sizeof(digit_t) * sz);
            ptr->m_size = sz;
            ptr->m_capacity = sz;
            deallocate(target);
            target.m_val    = 1;
            target.m_ptr    = ptr;
            target.m_kind   = mpz_large;
            target.m_owner  = mpz_self;
        }
        else {
            target.m_ptr->m_size = sz;
            if (target.m_ptr->m_digits != digits)
                memcpy(target.m_ptr->m_digits, digits, sizeof(digit_t) * sz);
            target.m_kind        = mpz_large;
        }
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::mul(mpz const & a, mpz const & b, mpz & c) {
    if (is_small(a) && is_small(b)) {
        set_i64(c, i64(a) * i64(b));
    }
    else {
        big_mul(a, b, c);
    }
}

// d <- a + b*c
template<bool SYNCH>
void mpz_manager<SYNCH>::addmul(mpz const & a, mpz const & b, mpz const & c, mpz & d) {
    if (is_one(b)) {
        add(a, c, d);
    }
    else if (is_minus_one(b)) {
        sub(a, c, d);
    }
    else {
        mpz tmp;
        mul(b,c,tmp);
        add(a,tmp,d);
        del(tmp);
    }
}


// d <- a - b*c
template<bool SYNCH>
void mpz_manager<SYNCH>::submul(mpz const & a, mpz const & b, mpz const & c, mpz & d) {
    if (is_one(b)) {
        sub(a, c, d);
    }
    else if (is_minus_one(b)) {
        add(a, c, d);
    }
    else {
        mpz tmp;
        mul(b,c,tmp);
        sub(a,tmp,d);
        del(tmp);
    }
}


template<bool SYNCH>
void mpz_manager<SYNCH>::machine_div_rem(mpz const & a, mpz const & b, mpz & q, mpz & r) {
    if (is_small(a) && is_small(b)) {
        int64_t _a = i64(a);
        int64_t _b = i64(b);
        set_i64(q, _a / _b);
        set_i64(r, _a % _b);
    }
    else {
        big_div_rem(a, b, q, r);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::machine_div(mpz const & a, mpz const & b, mpz & c) {
    if (is_small(b) && i64(b) == 0)
        throw std::runtime_error("division by 0");

    if (is_small(a) && is_small(b))
        set_i64(c, i64(a) / i64(b));
    else
        big_div(a, b, c);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::reset(mpz & a) {
    deallocate(a);
    set(a, 0);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::rem(mpz const & a, mpz const & b, mpz & c) {
    if (is_small(a) && is_small(b)) {
        set_i64(c, i64(a) % i64(b));
    }
    else {
        big_rem(a, b, c);
    }
}


template<bool SYNCH>
void mpz_manager<SYNCH>::div_gcd(mpz const& a, mpz const& b, mpz & c) {
    if (is_one(b)) {
        set(c, a);
    }
    else {
        machine_div(a, b, c);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::div(mpz const & a, mpz const & b, mpz & c) {
    assert(!is_zero(b));
    if (is_one(b)) {
        set(c, a);
    }
    else if (is_neg(a)) {
        mpz tmp;
        machine_div_rem(a, b, c, tmp);
        if (!is_zero(tmp)) {
            if (is_neg(b))
                add(c, mk_z(1), c);
            else
                sub(c, mk_z(1), c);
        }
        del(tmp);
    }
    else {
        machine_div(a, b, c);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::mod(mpz const & a, mpz const & b, mpz & c) {
    rem(a, b, c);
    if (is_neg(c)) {
        if (is_pos(b))
            add(c, b, c);
        else
            sub(c, b, c);
    }
}

template<bool SYNCH>
mpz mpz_manager<SYNCH>::mod2k(mpz const & a, unsigned k) {
    if (is_zero(a))
        return 0;

    mpz result;

    if (is_small(a) && k < 64) {
        uint64_t mask = ((1ULL << k) - 1);
        uint64_t uval = static_cast<uint64_t>(i64(a));
        set_i64(result, static_cast<int64_t>(uval & mask));
        return result;
    }

    if (is_nonneg(a) && bitsize(a) <= k) {
        return dup(a);
    }

    sign_cell ca(*this, a);
    unsigned digit_size = sizeof(digit_t) * 8;
    unsigned digit_count = k / digit_size;
    unsigned rem_bits = k % digit_size;
    unsigned total_digits = digit_count + (rem_bits > 0);
    digit_t mask = (1ULL << rem_bits) - 1;
    bool is_zero = true;

    allocate_if_needed(result, total_digits);

    // compute |a| mod 2^k-
    for (unsigned i = 0, e = std::min(digit_count, ca.cell()->m_size); i < e; ++i) {
        is_zero &= (digits(result)[i] = ca.cell()->m_digits[i]) == 0;
    }
    for (unsigned i = ca.cell()->m_size; i < total_digits; ++i) {
        digits(result)[i] = 0;
    }

    if (rem_bits > 0 && digit_count < ca.cell()->m_size) {
        is_zero &= (digits(result)[digit_count] = ca.cell()->m_digits[digit_count] & mask) == 0;
    }
    result.m_ptr->m_size = total_digits;

    if (ca.sign() < 0 && !is_zero) {
        // Negative case: if non-zero, result = 2^k - (|a| mod 2^k)
        // which boils down to computing ~result + 1
        for (unsigned i = 0; i < total_digits; ++i) {
            digits(result)[i] = ~digits(result)[i];
        }

        // Increment result
        digit_t carry = 1;
        for (unsigned i = 0; i < total_digits && carry; ++i) {
            digit_t sum = digits(result)[i] + carry;
            carry = sum < digits(result)[i];
            digits(result)[i] = sum;
        }

        // Clamp to k bits
        if (rem_bits != 0) {
            digits(result)[digit_count] &= mask;
        }
    }
    normalize(result);
    return result;
}

template<bool SYNCH>
void mpz_manager<SYNCH>::neg(mpz & a) {
    if (is_small(a) && a.m_val == INT_MIN) {
        // neg(INT_MIN) is not a small int
        set_big_i64(a, - static_cast<long long>(INT_MIN));
        return;
    }
    a.m_val = -a.m_val;
}

template<bool SYNCH>
void mpz_manager<SYNCH>::abs(mpz & a) {
    if (is_small(a)) {
        if (a.m_val < 0) {
            if (a.m_val == INT_MIN) {
                // abs(INT_MIN) is not a small int
                set_big_i64(a, - static_cast<long long>(INT_MIN));
            }
            else
                a.m_val = -a.m_val;
        }
    }
    else {
        a.m_val = 1;
    }
}


// TBD: replace use of 'tmp' by 'c'.
template<bool SYNCH>
template<bool SUB>
void mpz_manager<SYNCH>::big_add_sub(mpz const & a, mpz const & b, mpz & c) {
    sign_cell ca(*this, a), cb(*this, b);
    int sign_b = cb.sign();
    mpz_stack tmp;
    if (SUB)
        sign_b = -sign_b;
    unsigned real_sz;
    if (ca.sign() == sign_b) {
        unsigned sz  = std::max(ca.cell()->m_size, cb.cell()->m_size)+1;
        allocate_if_needed(tmp, sz);
        m_mpn_manager.add(ca.cell()->m_digits, ca.cell()->m_size,
                          cb.cell()->m_digits, cb.cell()->m_size,
                          tmp.m_ptr->m_digits, sz, &real_sz);
        assert(real_sz <= sz);
        set(*tmp.m_ptr, c, ca.sign(), real_sz);
    }
    else {
        digit_t borrow;
        int r = m_mpn_manager.compare(ca.cell()->m_digits, ca.cell()->m_size,
                                      cb.cell()->m_digits, cb.cell()->m_size);
        if (r == 0) {
            set(c, 0);
        }
        else if (r < 0) {
            // a < b
            unsigned sz = cb.cell()->m_size;
            allocate_if_needed(tmp, sz);
            m_mpn_manager.sub(cb.cell()->m_digits,
                              cb.cell()->m_size,
                              ca.cell()->m_digits,
                              ca.cell()->m_size,
                              tmp.m_ptr->m_digits,
                              &borrow);
            assert(borrow == 0);
            set(*tmp.m_ptr, c, sign_b, sz);
        }
        else {
            // a > b
            unsigned sz = ca.cell()->m_size;
            allocate_if_needed(tmp, sz);
            m_mpn_manager.sub(ca.cell()->m_digits,
                              ca.cell()->m_size,
                              cb.cell()->m_digits,
                              cb.cell()->m_size,
                              tmp.m_ptr->m_digits,
                              &borrow);
            assert(borrow == 0);
            set(*tmp.m_ptr, c, ca.sign(), sz);
        }
    }
    del(tmp);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::big_add(mpz const & a, mpz const & b, mpz & c) {
    big_add_sub<false>(a, b, c);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::big_sub(mpz const & a, mpz const & b, mpz & c) {
    big_add_sub<true>(a, b, c);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::big_mul(mpz const & a, mpz const & b, mpz & c) {
    // TBD replace tmp by c.
    mpz_stack tmp;
    sign_cell ca(*this, a), cb(*this, b);
    unsigned sz  = ca.cell()->m_size + cb.cell()->m_size;
    allocate_if_needed(tmp, sz);
    m_mpn_manager.mul(ca.cell()->m_digits,
                      ca.cell()->m_size,
                      cb.cell()->m_digits,
                      cb.cell()->m_size,
                      tmp.m_ptr->m_digits);
    set(*tmp.m_ptr, c, ca.sign() == cb.sign() ? 1 : -1, sz);
    del(tmp);
}


template<bool SYNCH>
void mpz_manager<SYNCH>::big_div_rem(mpz const & a, mpz const & b, mpz & q, mpz & r) {
    quot_rem_core<QUOT_AND_REM>(a, b, q, r);
}

template<bool SYNCH>
template<qr_mode MODE>
void mpz_manager<SYNCH>::quot_rem_core(mpz const & a, mpz const & b, mpz & q, mpz & r)
{
    /*
      +26 / +7 = +3, remainder is +5
      -26 / +7 = -3, remainder is -5
      +26 / -7 = -3, remainder is +5
      -26 / -7 = +3, remainder is -5
    */
    mpz_stack q1, r1;
    sign_cell ca(*this, a), cb(*this, b);
    if (cb.cell()->m_size > ca.cell()->m_size) {
        if (MODE == REM_ONLY || MODE == QUOT_AND_REM)
            set(r, a);
        if (MODE == QUOT_ONLY || MODE == QUOT_AND_REM)
            set(q, 0);
        return;
    }
    unsigned q_sz = ca.cell()->m_size - cb.cell()->m_size + 1;
    unsigned r_sz = cb.cell()->m_size;
    allocate_if_needed(q1, q_sz);
    allocate_if_needed(r1, r_sz);
    m_mpn_manager.div(ca.cell()->m_digits, ca.cell()->m_size,
                      cb.cell()->m_digits, cb.cell()->m_size,
                      q1.m_ptr->m_digits,
                      r1.m_ptr->m_digits);
    if (MODE == QUOT_ONLY || MODE == QUOT_AND_REM)
        set(*q1.m_ptr, q, ca.sign() == cb.sign() ? 1 : -1, q_sz);
    if (MODE == REM_ONLY || MODE == QUOT_AND_REM)
        set(*r1.m_ptr, r, ca.sign(), r_sz);
    del(q1);
    del(r1);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::big_div(mpz const & a, mpz const & b, mpz & c) {
    mpz dummy;
    quot_rem_core<QUOT_ONLY>(a, b, c, dummy);
    assert(is_zero(dummy));
    del(dummy);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::big_rem(mpz const & a, mpz const & b, mpz & c) {
    mpz dummy;
    quot_rem_core<REM_ONLY>(a, b, dummy, c);
    assert(is_zero(dummy));
    del(dummy);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::gcd(mpz const & a, mpz const & b, mpz & c) {
    static_assert(sizeof(a.m_val) == sizeof(int), "size mismatch");
    static_assert(sizeof(mpz) <= 16, "mpz size overflow");
    if (is_small(a) && is_small(b) && a.m_val != INT_MIN && b.m_val != INT_MIN) {
        int _a = a.m_val;
        int _b = b.m_val;
        if (_a < 0) _a = -_a;
        if (_b < 0) _b = -_b;
        unsigned r = u_gcd(_a, _b);
        set(c, r);
    }
    else {
        if (is_zero(a)) {
            set(c, b);
            abs(c);
            return;
        }
        if (is_zero(b)) {
            set(c, a);
            abs(c);
            return;
        }
#ifdef BINARY_GCD
        // Binary GCD for big numbers
        // - It doesn't use division
        // - The initial experiments, don't show any performance improvement
        // - It only works with _MP_INTERNAL
        mpz u, v, diff;
        set(u, a);
        set(v, b);
        abs(u);
        abs(v);

        unsigned k_u = power_of_two_multiple(u);
        unsigned k_v = power_of_two_multiple(v);
        unsigned k   = k_u < k_v ? k_u : k_v;

        machine_div2k(u, k_u);

        while (true) {
            machine_div2k(v, k_v);

            if (lt(u, v)) {
                sub(v, u, v);
            }
            else {
                sub(u, v, diff);
                swap(u, v);
                swap(v, diff);
            }

            if (is_zero(v) || is_one(v))
                break;

            // reset least significant bit
            if (is_small(v))
                v.m_val &= ~1;
            else
                v.m_ptr->m_digits[0] &= ~static_cast<digit_t>(1);
            k_v = power_of_two_multiple(v);
        }

        mul2k(u, k, c);
        del(u); del(v); del(diff);
#endif // BINARY_GCD

#ifdef EUCLID_GCD
        mpz tmp1;
        mpz tmp2;
        mpz aux;
        set(tmp1, a);
        set(tmp2, b);
        abs(tmp1);
        abs(tmp2);
        if (lt(tmp1, tmp2))
            swap(tmp1, tmp2);
        if (is_zero(tmp2)) {
            swap(c, tmp1);
        }
        else {
            while (true) {
                if (is_uint64(tmp1) && is_uint64(tmp2)) {
                    set(c, u64_gcd(get_uint64(tmp1), get_uint64(tmp2)));
                    break;
                }
                rem(tmp1, tmp2, aux);
                if (is_zero(aux)) {
                    swap(c, tmp2);
                    break;
                }
                swap(tmp1, tmp2);
                swap(tmp2, aux);
            }
        }
        del(tmp1); del(tmp2); del(aux);
#endif // EUCLID_GCD

#ifdef LS_BINARY_GCD
        mpz u, v, t, u1, u2;
        set(u, a);
        set(v, b);
        abs(u);
        abs(v);
        if (lt(u, v))
            swap(u, v);
        while (!is_zero(v)) {
            // Basic idea:
            // compute t = 2^e*v  such that t <= u < 2t
            // u := min{u - t, 2t - u}
            //
            // The assignment u := min{u - t, 2t - u}
            // can be replaced with u := u - t
            //
            // Since u and v are positive, we have:
            //    2^{log2(u)}     <= u < 2^{(log2(u) + 1)}
            //    2^{log2(v)}     <= v < 2^{(log2(v) + 1)}
            //  -->
            //    2^{log2(v)}*2^{log2(u)-log2(v)} <= v*2^{log2(u)-log2(v)} < 2^{log2(v) + 1}*2^{log2(u)-log2(v)}
            //  -->
            //    2^{log2(u)} <= v*2^{log2(u)-log2(v)} < 2^{log2(u) + 1}
            //
            // Now, let t be v*2^{log2(u)-log2(v)}
            // If t <= u, then we found t
            // Otherwise t = t div 2
            unsigned k_u = log2(u);
            unsigned k_v = log2(v);
            assert(k_v <= k_u);
            unsigned e   = k_u - k_v;
            mul2k(v, e, t);
            sub(u, t, u1);
            if (is_neg(u1)) {
                // t is too big
                machine_div2k(t, 1);
                // Now, u1 contains u - 2t
                neg(u1);
                // Now, u1 contains 2t - u
                sub(u, t, u2); // u2 := u - t
            }
            else {
                // u1 contains u - t
                mul2k(t, 1);
                sub(t, u, u2);
                // u2 contains 2t - u
            }
            assert(is_nonneg(u1));
            assert(is_nonneg(u2));
            if (lt(u1, u2))
                swap(u, u1);
            else
                swap(u, u2);
            if (lt(u, v))
                swap(u,v);
        }
        swap(u, c);
        del(u); del(v); del(t); del(u1); del(u2);
#endif // LS_BINARY_GCD

#ifdef LEHMER_GCD
        // For now, it only works if sizeof(digit_t) == sizeof(unsigned)
        static_assert(sizeof(digit_t) == sizeof(unsigned), "");

        int64_t a_hat, b_hat, A, B, C, D, T, q, a_sz, b_sz;
        mpz a1, b1, t, r, tmp;
        set(a1, a);
        set(b1, b);
        abs(a1);
        abs(b1);
        if (lt(a1, b1))
            swap(a1, b1);
        while (true) {
            assert(ge(a1, b1));
            if (is_small(b1)) {
                if (is_small(a1)) {
                    unsigned r = u_gcd(a1.m_val, b1.m_val);
                    set(c, r);
                    break;
                }
                else {
                    while (!is_zero(b1)) {
                        assert(ge(a1, b1));
                        rem(a1, b1, tmp);
                        swap(a1, b1);
                        swap(b1, tmp);
                    }
                    swap(c, a1);
                    break;
                }
            }
            assert(!is_small(a1));
            assert(!is_small(b1));
            a_sz  = a1.m_ptr->m_size;
            b_sz  = b1.m_ptr->m_size;
            assert(b_sz <= a_sz);
            a_hat = a1.m_ptr->m_digits[a_sz - 1];
            b_hat = (b_sz == a_sz) ? b1.m_ptr->m_digits[b_sz - 1] : 0;
            A = 1;
            B = 0;
            C = 0;
            D = 1;
            while (true) {
                // Loop invariants
                assert(a_hat + A <= static_cast<int64_t>(UINT_MAX) + 1);
                assert(a_hat + B <  static_cast<int64_t>(UINT_MAX) + 1);
                assert(b_hat + C <  static_cast<int64_t>(UINT_MAX) + 1);
                assert(b_hat + D <= static_cast<int64_t>(UINT_MAX) + 1);
                // overflows can't happen since I'm using int64
                if (b_hat + C == 0 || b_hat + D == 0)
                    break;
                q  = (a_hat + A)/(b_hat + C);
                if (q != (a_hat + B)/(b_hat + D))
                    break;
                T = A - q*C;
                A = C;
                C = T;
                T = B - q*D;
                B = D;
                D = T;
                T = a_hat - q*b_hat;
                a_hat = b_hat;
                b_hat = T;
            }
            assert(ge(a1, b1));
            if (B == 0) {
                rem(a1, b1, t);
                swap(a1, b1);
                swap(b1, t);
                assert(ge(a1, b1));
            }
            else {
                // t <- A*a1
                set(tmp, A);
                mul(a1, tmp, t);
                // t <- t + B*b1
                set(tmp, B);
                addmul(t, tmp, b1, t);
                // r <- C*a1
                set(tmp, C);
                mul(a1, tmp, r);
                // r <- r + D*b1
                set(tmp, D);
                addmul(r, tmp, b1, r);
                // a <- t
                swap(a1, t);
                // b <- r
                swap(b1, r);
                assert(ge(a1, b1));
            }
        }
        del(a1); del(b1); del(r); del(t); del(tmp);
#endif // LEHMER_GCD
    }
}

template<bool SYNCH>
unsigned mpz_manager<SYNCH>::size_info(mpz const & a) {
    if (is_small(a))
        return 1;
    return a.m_ptr->m_size + 1;
}



template<bool SYNCH>
struct mpz_manager<SYNCH>::sz_lt {
    mpz const * m_as;
    bool operator()(unsigned p1, unsigned p2) {
        return size_info(m_as[p1]) < size_info(m_as[p2]);
    }
};


template<bool SYNCH>
void mpz_manager<SYNCH>::gcd(unsigned sz, mpz const * as, mpz & g) {
#if 0
    // Optimization: sort numbers by size. Motivation: compute the gcd of the small ones first.
    // The optimization did not really help.
    switch (sz) {
    case 0:
        set(g, 0);
        return;
    case 1:
        set(g, as[0]);
        abs(g);
        return;
    case 2:
        gcd(as[0], as[1], g);
        return;
    default:
        break;
    }
    unsigned i;
    for (i = 0; i < sz; ++i) {
        if (!is_small(as[i]))
            break;
    }
    if (i != sz) {
        // array has big numbers
        sbuffer<unsigned, 1024> p;
        for (i = 0; i < sz; ++i)
            p.push_back(i);
        std::sort(p.begin(), p.end(), sz_lt{as});
        gcd(as[p[0]], as[p[1]], g);
        for (i = 2; i < sz; ++i) {
            if (is_one(g))
                return;
            gcd(g, as[p[i]], g);
        }
        return;
    }
    else {
        gcd(as[0], as[1], g);
        for (unsigned i = 2; i < sz; ++i) {
            if (is_one(g))
                return;
            gcd(g, as[i], g);
        }
    }
#else
    // Vanilla implementation
    switch (sz) {
    case 0:
        set(g, 0);
        return;
    case 1:
        set(g, as[0]);
        abs(g);
        return;
    default:
        break;
    }
    gcd(as[0], as[1], g);
    for (unsigned i = 2; i < sz; ++i) {
        if (is_one(g))
            return;
        gcd(g, as[i], g);
    }
#endif
}

template<bool SYNCH>
void mpz_manager<SYNCH>::gcd(mpz const & r1, mpz const & r2, mpz & a, mpz & b, mpz & r) {
    mpz tmp1, tmp2;
    mpz aux, quot;
    set(tmp1, r1);
    set(tmp2, r2);
    set(a, 1);
    set(b, 0);
    mpz nexta, nextb;
    set(nexta, 0);
    set(nextb, 1);

    abs(tmp1);
    abs(tmp2);
    if (lt(tmp1, tmp2)) {
        swap(tmp1, tmp2);
        swap(nexta, nextb);
        swap(a, b);
    }

    // tmp1 >= tmp2 >= 0
    // quot_rem in one function would be faster.
    while (is_pos(tmp2)) {
        assert(ge(tmp1, tmp2));

        // aux = tmp2
        set(aux, tmp2);
        // quot = div(tmp1, tmp2);
        machine_div(tmp1, tmp2, quot);
        // tmp2 = tmp1 % tmp2
        rem(tmp1, tmp2, tmp2);
        // tmp1 = aux
        set(tmp1, aux);
        // aux = nexta
        set(aux, nexta);
        // nexta = a - (quot*nexta)
        mul(quot, nexta, nexta);
        sub(a, nexta, nexta);
        // a = axu
        set(a, aux);
        // aux = nextb
        set(aux, nextb);
        // nextb = b - (quot*nextb)
        mul(nextb, quot, nextb);
        sub(b, nextb, nextb);
        // b = aux
        set(b, aux);
    }

    if (is_neg(r1))
        neg(a);
    if (is_neg(r2))
        neg(b);

    set(r, tmp1);
    del(tmp1);
    del(tmp2);
    del(aux);
    del(quot);
    del(nexta);
    del(nextb);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::lcm(mpz const & a, mpz const & b, mpz & c) {
    if (is_one(b)) {
        set(c, a);
    }
    else if (is_one(a) || eq(a, b)) {
        set(c, b);
    }
    else {
        mpz r;
        gcd(a, b, r);
        if (eq(r, a)) {
            set(c, b);
        }
        else if (eq(r, b)) {
            set(c, a);
        }
        else {
            // c contains gcd(a, b)
            // so c divides a, and machine_div(a, c) is equal to div(a, c)
            machine_div(a, r, r);
            mul(r, b, c);
        }
        del(r);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::bitwise_or(mpz const & a, mpz const & b, mpz & c) {
    assert(is_nonneg(a));
    assert(is_nonneg(b));
    if (is_small(a) && is_small(b)) {
        c.m_val = a.m_val | b.m_val;
        c.m_kind = mpz_small;
    }
    else {
        mpz a1, b1, a2, b2, m, tmp;
        set(a1, a);
        set(b1, b);
        set(m, 1);
        set(c, 0);
        while (!is_zero(a1) && !is_zero(b1)) {
            mod(a1, m_two64, a2);
            mod(b1, m_two64, b2);
            uint64_t v = get_uint64(a2) | get_uint64(b2);
            set(tmp, v);
            mul(tmp, m, tmp);
            add(c, tmp, c); // c += m * v
            mul(m, m_two64, m);
            div(a1, m_two64, a1);
            div(b1, m_two64, b1);
        }
        if (!is_zero(a1)) {
            mul(a1, m, a1);
            add(c, a1, c);
        }
        if (!is_zero(b1)) {
            mul(b1, m, b1);
            add(c, b1, c);
        }
        del(a1); del(b1); del(a2); del(b2); del(m); del(tmp);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::bitwise_and(mpz const & a, mpz const & b, mpz & c) {
    if (is_small(a) && is_small(b)) {
        c.m_val = a.m_val & b.m_val;
        c.m_kind = mpz_small;
    }
    else {
        mpz a1, b1, a2, b2, m, tmp;
        set(a1, a);
        set(b1, b);
        set(m, 1);
        set(c, 0);
        while (!is_zero(a1) && !is_zero(b1)) {
            mod(a1, m_two64, a2);
            mod(b1, m_two64, b2);
            uint64_t v = get_uint64(a2) & get_uint64(b2);
            set(tmp, v);
            mul(tmp, m, tmp);
            add(c, tmp, c); // c += m * v
            mul(m, m_two64, m);
            div(a1, m_two64, a1);
            div(b1, m_two64, b1);
        }
        del(a1); del(b1); del(a2); del(b2); del(m); del(tmp);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::bitwise_xor(mpz const & a, mpz const & b, mpz & c) {
    assert(is_nonneg(a));
    assert(is_nonneg(b));
    if (is_small(a) && is_small(b)) {
        set_i64(c, i64(a) ^ i64(b));
    }
    else {
        mpz a1, b1, a2, b2, m, tmp;
        set(a1, a);
        set(b1, b);
        set(m, 1);
        set(c, 0);
        while (!is_zero(a1) && !is_zero(b1)) {
            mod(a1, m_two64, a2);
            mod(b1, m_two64, b2);
            uint64_t v = get_uint64(a2) ^ get_uint64(b2);
            set(tmp, v);
            mul(tmp, m, tmp);
            add(c, tmp, c); // c += m * v
            mul(m, m_two64, m);
            div(a1, m_two64, a1);
            div(b1, m_two64, b1);
        }
        if (!is_zero(a1)) {
            mul(a1, m, a1);
            add(c, a1, c);
        }
        if (!is_zero(b1)) {
            mul(b1, m, b1);
            add(c, b1, c);
        }
        del(a1); del(b1); del(a2); del(b2); del(m); del(tmp);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::bitwise_not(unsigned sz, mpz const & a, mpz & c) {
    assert(is_nonneg(a));
    if (is_small(a) && sz <= 64) {
        uint64_t v = ~get_uint64(a);
        unsigned zero_out = 64 - sz;
        v = (v << zero_out) >> zero_out;
        set(c, v);
    }
    else {
        mpz a1, a2, m, tmp;
        set(a1, a);
        set(m, 1);
        set(c, 0);
        while (sz > 0) {
            mod(a1, m_two64, a2);
            uint64_t n = get_uint64(a2);
            uint64_t v = ~n;
            assert(~v == n);
            if (sz < 64) {
                uint64_t mask = (1ull << static_cast<uint64_t>(sz)) - 1ull;
                v = mask & v;
            }
            set(tmp, v);
            assert(get_uint64(tmp) == v);
            mul(tmp, m, tmp);
            add(c, tmp, c); // c += m * v
            mul(m, m_two64, m);
            div(a1, m_two64, a1);
            sz -= (sz<64) ? sz : 64;
        }
        del(a1); del(a2); del(m); del(tmp);
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::big_set(mpz & target, mpz const & source) {
    if (&target == &source)
        return;
    target.m_val = source.m_val;
    if (target.m_ptr == nullptr) {
        target.m_ptr = allocate(capacity(source));
        target.m_ptr->m_size     = size(source);
        target.m_ptr->m_capacity = capacity(source);
        target.m_kind = mpz_large;
        target.m_owner = mpz_self;
        memcpy(target.m_ptr->m_digits, source.m_ptr->m_digits, sizeof(digit_t) * size(source));
    }
    else if (capacity(target) < size(source)) {
        deallocate(target);
        target.m_ptr = allocate(capacity(source));
        target.m_ptr->m_size     = size(source);
        target.m_ptr->m_capacity = capacity(source);
        target.m_kind = mpz_large;
        target.m_owner = mpz_self;
        memcpy(target.m_ptr->m_digits, source.m_ptr->m_digits, sizeof(digit_t) * size(source));
    }
    else {
        target.m_ptr->m_size = size(source);
        memcpy(target.m_ptr->m_digits, source.m_ptr->m_digits, sizeof(digit_t) * size(source));
        target.m_kind = mpz_large;
    }
}

template<bool SYNCH>
int mpz_manager<SYNCH>::big_compare(mpz const & a, mpz const & b) {
    if (sign(a) > 0) {
        // a is positive
        if (sign(b) > 0) {
            // a & b are positive
            sign_cell ca(*this, a), cb(*this, b);
            return m_mpn_manager.compare(ca.cell()->m_digits, ca.cell()->m_size,
                                         cb.cell()->m_digits, cb.cell()->m_size);
        }
        else {
            // b is negative
            return 1; // a > b
        }
    }
    else {
        // a is negative
        if (sign(b) > 0) {
            // b is positive
            return -1; // a < b
        }
        else {
            // a & b are negative
            sign_cell ca(*this, a), cb(*this, b);
            return m_mpn_manager.compare(cb.cell()->m_digits, cb.cell()->m_size,
                                         ca.cell()->m_digits, ca.cell()->m_size);
        }
    }
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::is_uint64(mpz const & a) const {
    if (a.m_val < 0)
        return false;
    if (is_small(a))
        return true;
    if (sizeof(digit_t) == sizeof(uint64_t)) {
        return size(a) <= 1;
    }
    else {
        return size(a) <= 2;
    }
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::is_int64(mpz const & a) const {
    if (is_small(a))
        return true;
    if (!is_abs_uint64(a))
        return false;
    uint64_t num = big_abs_to_uint64(a);
    uint64_t msb = static_cast<uint64_t>(1) << 63;
    uint64_t msb_val = msb & num;
    if (a.m_val >= 0) {
        // non-negative number.
        return (0 == msb_val);
    }
    else {
        // negative number.
        // either the high bit is 0, or
        // the number is 2^64 which can be represented.
        //
        return 0 == msb_val || (msb_val == num);
    }
}

template<bool SYNCH>
uint64_t mpz_manager<SYNCH>::get_uint64(mpz const & a) const {
    if (is_small(a))
        return static_cast<uint64_t>(a.m_val);
    assert(a.m_ptr->m_size > 0);
    return big_abs_to_uint64(a);
}

template<bool SYNCH>
int64_t mpz_manager<SYNCH>::get_int64(mpz const & a) const {
    if (is_small(a))
        return static_cast<int64_t>(a.m_val);
    assert(is_int64(a));
    uint64_t num = big_abs_to_uint64(a);
    if (a.m_val < 0) {
        if (num != 0 && (num << 1) == 0)
            return INT64_MIN;
        return -static_cast<int64_t>(num);
    }
    return static_cast<int64_t>(num);
}

template<bool SYNCH>
double mpz_manager<SYNCH>::get_double(mpz const & a) const {
    if (is_small(a))
        return static_cast<double>(a.m_val);
    double r = 0.0;
    double d = 1.0;
    unsigned sz = size(a);
    for (unsigned i = 0; i < sz; ++i) {
        r += d * static_cast<double>(digits(a)[i]);
        if (sizeof(digit_t) == sizeof(uint64_t))
            d *= (1.0 + static_cast<double>(UINT64_MAX)); // 64-bit version, multiply by 2^64
        else
            d *= (1.0 + static_cast<double>(UINT_MAX));   // 32-bit version, multiply by 2^32
    }
    if (!(r >= 0.0)) {
        r = static_cast<double>(UINT64_MAX); // some large number
    }
    return a.m_val < 0 ? -r : r;
}

template<bool SYNCH>
void mpz_manager<SYNCH>::display(std::ostream & out, mpz const & a) const {
    if (is_small(a)) {
        out << a.m_val;
    }
    else {
        if (a.m_val < 0)
            out << '-';

        auto sz = sizeof(digit_t) == 4 ? 11 : 21;
        std::vector<char> buffer(sz * size(a), 0);
        out << m_mpn_manager.to_string(digits(a), size(a), buffer.data(), buffer.size());
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::display_smt2(std::ostream & out, mpz const & a, bool decimal) const {
    if (is_neg(a)) {
        mpz_manager<SYNCH> * _this = const_cast<mpz_manager<SYNCH>*>(this);
        _scoped_numeral<mpz_manager<SYNCH> > tmp(*_this);
        _this->set(tmp, a);
        _this->neg(tmp);
        out << "(- ";
        display(out, tmp);
        if (decimal)
            out << ".0";
        out << ")";
    }
    else {
        display(out, a);
        if (decimal)
            out << ".0";
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::display_hex(std::ostream & out, mpz const & a, unsigned num_bits) const {
    assert(num_bits % 4 == 0);
    std::ios fmt(nullptr);
    fmt.copyfmt(out);
    out << std::hex;
    if (is_small(a)) {
        out << std::setw(num_bits/4) << std::setfill('0') << get_uint64(a);
    } else {
        digit_t *ds = digits(a);
        unsigned sz = size(a);
        unsigned bitSize = sz * sizeof(digit_t) * 8;
        unsigned firstDigitSize;
        if (num_bits >= bitSize) {
            firstDigitSize = sizeof(digit_t) * 2;

            for (unsigned i = 0; i < (num_bits - bitSize)/4; ++i) {
                out << "0";
            }
        } else {
            firstDigitSize = num_bits % (sizeof(digit_t) * 8) / 4;
        }

        out << std::setfill('0') << std::setw(firstDigitSize) << ds[sz-1] << std::setw(sizeof(digit_t)*2);
        for (unsigned i = 1; i < sz; ++i) {
            out << ds[sz-i-1];
        }
    }
    out.copyfmt(fmt);
}

static void display_binary_data(std::ostream &out, uint64_t val, uint64_t numBits) {
    for (uint64_t shift = numBits; shift-- > 64ull; ) out << "0";
    if (numBits > 64) numBits = 64;
    for (uint64_t shift = numBits; shift-- > 0; ) {
        if (val & (1ull << shift)) {
            out << "1";
        } else {
            out << "0";
        }
    }
 }

template<bool SYNCH>
void mpz_manager<SYNCH>::display_bin(std::ostream & out, mpz const & a, unsigned num_bits) const {
    if (is_small(a)) {
        display_binary_data(out, get_uint64(a), num_bits);
    }
    else {
        digit_t *ds = digits(a);
        unsigned sz = size(a);
        const unsigned digitBitSize = sizeof(digit_t) * 8;
        unsigned bitSize = sz * digitBitSize;
        unsigned firstDigitLength;
        if (num_bits > bitSize) {
            firstDigitLength = 0;
            for (unsigned i = 0; i < (num_bits - bitSize); ++i) {
                out << "0";
            }
        } else {
            firstDigitLength = num_bits % digitBitSize;
        }
        for (unsigned i = 0; i < sz; ++i) {
            if (i == 0 && firstDigitLength != 0) {
                display_binary_data(out, ds[sz-1], firstDigitLength);
            } else {
                display_binary_data(out, ds[sz-i-1], digitBitSize);
            }
        }
    }
}

template<bool SYNCH>
std::string mpz_manager<SYNCH>::to_string(mpz const & a) const {
    std::ostringstream buffer;
    display(buffer, a);
    return std::move(buffer).str();
}

template<bool SYNCH>
unsigned mpz_manager<SYNCH>::hash(mpz const & a) {
    if (is_small(a))
        return ::abs(a.m_val);
    unsigned sz = size(a);
    if (sz == 1)
        return static_cast<unsigned>(digits(a)[0]);
    return string_hash(std::string_view(reinterpret_cast<char*>(digits(a)), sz * sizeof(digit_t)), 17);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::power(mpz const & a, unsigned p, mpz & b) {
    if (is_small(a)) {
        if (a.m_val == 2) {
            if (p < 8 * sizeof(int) - 1) {
                b.m_val = 1 << p;
                b.m_kind = mpz_small;
            }
            else {
                unsigned sz    = p/(8 * sizeof(digit_t)) + 1;
                unsigned shift = p%(8 * sizeof(digit_t));
                assert(sz > 0);
                allocate_if_needed(b, sz);
                assert(b.m_ptr->m_capacity >= sz);
                b.m_ptr->m_size     = sz;
                for (unsigned i = 0; i < sz - 1; ++i)
                    b.m_ptr->m_digits[i] = 0;
                b.m_ptr->m_digits[sz-1] = 1 << shift;
                b.m_val = 1;
                b.m_kind = mpz_large;
            }
            return;
        }
        if (a.m_val == 0) {
            assert(p != 0);
            set(b, 0);
            return;
        }
        if (a.m_val == 1) {
            set(b, 1);
            return;
        }
    }
    // general purpose
    unsigned mask = 1;
    mpz power;
    set(power, a);
    set(b, 1);
    while (mask <= p) {
        if (mask & p)
            mul(b, power, b);
        mul(power, power, power);
        mask = mask << 1;
    }
    del(power);
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::is_power_of_two(mpz const & a) {
    unsigned shift;
    return is_power_of_two(a, shift);
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::is_power_of_two(mpz const & a, unsigned & shift) {
    if (is_nonpos(a))
        return false;
    if (is_small(a)) {
        if (::is_power_of_two(a.m_val)) {
            shift = ::log2((unsigned)a.m_val);
            return true;
        }
        else {
            return false;
        }
    }
    mpz_cell * c     = a.m_ptr;
    unsigned sz      = c->m_size;
    digit_t * ds     = c->m_digits;
    for (unsigned i = 0; i < sz - 1; ++i) {
        if (ds[i] != 0)
            return false;
    }
    digit_t v = ds[sz-1];
    if (!(v & (v - 1)) && v) {
        shift = log2(a);
        return true;
    }
    else {
        return false;
    }
}

// Expand capacity of a
template<bool SYNCH>
void mpz_manager<SYNCH>::ensure_capacity(mpz & a, unsigned capacity) {
    if (capacity <= 1)
        return;
    if (capacity < m_init_cell_capacity)
        capacity = m_init_cell_capacity;

    if (is_small(a)) {
        int val = a.m_val;
        allocate_if_needed(a, capacity);
        a.m_kind = mpz_large;
        assert(a.m_ptr->m_capacity >= capacity);
        if (val == INT_MIN) {
            unsigned intmin_sz = m_int_min.m_ptr->m_size;
            for (unsigned i = 0; i < intmin_sz; ++i)
                a.m_ptr->m_digits[i] = m_int_min.m_ptr->m_digits[i];
            a.m_val = -1;
            a.m_ptr->m_size = m_int_min.m_ptr->m_size;
        }
        else if (val < 0) {
            a.m_ptr->m_digits[0] = -val;
            a.m_val = -1;
            a.m_ptr->m_size = 1;
        }
        else {
            a.m_ptr->m_digits[0] =  val;
            a.m_val = 1;
            a.m_ptr->m_size = 1;
        }
    }
    else if (a.m_ptr->m_capacity < capacity) {
        mpz_cell * new_cell = allocate(capacity);
        assert(new_cell->m_capacity == capacity);
        unsigned old_sz  = a.m_ptr->m_size;
        new_cell->m_size = old_sz;
        for (unsigned i = 0; i < old_sz; ++i)
            new_cell->m_digits[i] = a.m_ptr->m_digits[i];
        deallocate(a);
        a.m_ptr = new_cell;
        a.m_owner = mpz_self;
        a.m_kind = mpz_large;
    }
}

template<bool SYNCH>
void mpz_manager<SYNCH>::normalize(mpz & a) {
    mpz_cell * c = a.m_ptr;
    digit_t * ds = c->m_digits;
    unsigned i = c->m_size;
    for (; i > 0; --i) {
        if (ds[i-1] != 0)
            break;
    }

    if (i == 0) {
        // a is zero...
        set(a, 0);
        return;
    }

    if (i == 1 && ds[0] <= INT_MAX) {
        // a is small
        int val = a.m_val < 0 ? -static_cast<int>(ds[0]) : static_cast<int>(ds[0]);
        a.m_val = val;
        a.m_kind = mpz_small;
        return;
    }
    // adjust size
    c->m_size = i;
}

template<bool SYNCH>
void mpz_manager<SYNCH>::machine_div2k(mpz & a, unsigned k) {
    if (k == 0 || is_zero(a))
        return;
    if (is_small(a)) {
        if (k < 32) {
            int64_t twok = 1ull << ((int64_t)k);
            int64_t val = a.m_val;
            a.m_val = (int)(val/twok);
        }
        else {
            a.m_val = 0;
        }
        return;
    }
    unsigned digit_shift = k / (8 * sizeof(digit_t));
    mpz_cell * c         = a.m_ptr;
    unsigned sz          = c->m_size;
    if (digit_shift >= sz) {
        set(a, 0);
        return;
    }
    unsigned bit_shift   = k % (8 * sizeof(digit_t));
    unsigned comp_shift  = (8 * sizeof(digit_t)) - bit_shift;
    unsigned new_sz      = sz - digit_shift;
    assert(new_sz >= 1);
    digit_t * ds = c->m_digits;
    if (new_sz < sz) {
        unsigned i       = 0;
        unsigned j       = digit_shift;
        if (bit_shift != 0) {
            for (; i < new_sz - 1; ++i, ++j) {
                ds[i] = ds[j];
                ds[i] >>= bit_shift;
                ds[i] |= (ds[j+1] << comp_shift);
            }
            ds[i] = ds[j];
            ds[i] >>= bit_shift;
        }
        else {
            for (; i < new_sz; ++i, ++j) {
                ds[i] = ds[j];
            }
        }
    }
    else {
        assert(new_sz == sz);
        assert(bit_shift != 0);
        unsigned i       = 0;
        for (; i < new_sz - 1; ++i) {
            ds[i] >>= bit_shift;
            ds[i] |= (ds[i+1] << comp_shift);
        }
        ds[i] >>= bit_shift;
    }

    c->m_size = new_sz;
    normalize(a);
}

template<bool SYNCH>
void mpz_manager<SYNCH>::mul2k(mpz & a, unsigned k) {
    if (k == 0 || is_zero(a))
        return;
    if (is_small(a) && k < 32) {
        set_i64(a, i64(a) * (static_cast<int64_t>(1) << k));
        return;
    }
    unsigned word_shift  = k / (8 * sizeof(digit_t));
    unsigned bit_shift   = k % (8 * sizeof(digit_t));
    unsigned old_sz      = is_small(a) ? 1 : a.m_ptr->m_size;
    unsigned new_sz      = old_sz + word_shift + 1;
    ensure_capacity(a, new_sz);
    assert(!is_small(a));
    mpz_cell * cell_a    = a.m_ptr;
    old_sz = cell_a->m_size;
    digit_t * ds         = cell_a->m_digits;
    for (unsigned i = old_sz; i < new_sz; ++i)
        ds[i] = 0;
    cell_a->m_size       = new_sz;

    if (word_shift > 0) {
        unsigned j = old_sz;
        unsigned i = old_sz + word_shift;
        while (j > 0) {
            --j; --i;
            ds[i] = ds[j];
        }
        while (i > 0) {
            --i;
            ds[i] = 0;
        }
    }
    if (bit_shift > 0) {
#ifdef NDEBUG
        for (unsigned i = 0; i < word_shift; ++i) {
            assert(ds[i] == 0);
        }
#endif
        unsigned comp_shift = (8 * sizeof(digit_t)) - bit_shift;
        digit_t prev = 0;
        for (unsigned i = word_shift; i < new_sz; ++i) {
            digit_t new_prev = (ds[i] >> comp_shift);
            ds[i] <<= bit_shift;
            ds[i] |= prev;
            prev = new_prev;
        }
    }
    normalize(a);
}

static_assert(sizeof(digit_t) == 4 || sizeof(digit_t) == 8, "");

template<bool SYNCH>
unsigned mpz_manager<SYNCH>::power_of_two_multiple(mpz const & a) {
    if (is_zero(a))
        return 0;
    if (is_small(a)) {
        unsigned r = 0;
        int v      = a.m_val;
#define COUNT_DIGIT_RIGHT_ZEROS()               \
        if (v % (1 << 16) == 0) {               \
            r += 16;                            \
            v /= (1 << 16);                     \
        }                                       \
        if (v % (1 << 8) == 0) {                \
            r += 8;                             \
            v /= (1 << 8);                      \
        }                                       \
        if (v % (1 << 4) == 0) {                \
            r += 4;                             \
            v /= (1 << 4);                      \
        }                                       \
        if (v % (1 << 2) == 0) {                \
            r += 2;                             \
            v /= (1 << 2);                      \
        }                                       \
        if (v % 2 == 0) {                       \
            r++;                                \
        }
        COUNT_DIGIT_RIGHT_ZEROS();
        return r;
    }
    mpz_cell * c        = a.m_ptr;
    unsigned sz         = c->m_size;
    unsigned r          = 0;
    digit_t * source    = c->m_digits;
    for (unsigned i = 0; i < sz; ++i) {
        if (source[i] != 0) {
            digit_t v = source[i];
            if (sizeof(digit_t) == 8) {
                // TODO: we can remove this if after we move to MPN
                // In MPN the digit_t is always an unsigned integer
                if (static_cast<uint64_t>(v) % (static_cast<uint64_t>(1) << 32) == 0) {
                    r += 32;
                    v = static_cast<digit_t>(static_cast<uint64_t>(v) / (static_cast<uint64_t>(1) << 32));
                }
            }
            COUNT_DIGIT_RIGHT_ZEROS();
            return r;
        }
        r += (8 * sizeof(digit_t));
    }
    return r;
}

template<bool SYNCH>
unsigned mpz_manager<SYNCH>::log2(mpz const & a) {
    if (is_nonpos(a))
        return 0;
    if (is_small(a))
        return ::log2((unsigned)a.m_val);
    static_assert(sizeof(digit_t) == 8 || sizeof(digit_t) == 4, "");
    mpz_cell * c     = a.m_ptr;
    unsigned sz      = c->m_size;
    digit_t * ds     = c->m_digits;
    if (sizeof(digit_t) == 8)
        return (sz - 1)*64 + uint64_log2(ds[sz-1]);
    else
        return (sz - 1)*32 + ::log2(static_cast<unsigned>(ds[sz-1]));
}

template<bool SYNCH>
unsigned mpz_manager<SYNCH>::mlog2(mpz const & a) {
    if (is_nonneg(a))
        return 0;
    if (is_small(a) && a.m_val == INT_MIN)
        return ::log2((unsigned)a.m_val);

    if (is_small(a))
        return ::log2((unsigned)-a.m_val);
    static_assert(sizeof(digit_t) == 8 || sizeof(digit_t) == 4, "");
    mpz_cell * c     = a.m_ptr;
    unsigned sz      = c->m_size;
    digit_t * ds     = c->m_digits;
    if (sizeof(digit_t) == 8)
        return (sz - 1)*64 + uint64_log2(ds[sz-1]);
    else
        return (sz - 1)*32 + ::log2(static_cast<unsigned>(ds[sz-1]));
}

template<bool SYNCH>
unsigned mpz_manager<SYNCH>::bitsize(mpz const & a) {
    if (is_nonneg(a))
        return log2(a) + 1;
    else
        return mlog2(a) + 1;
}

template<bool SYNCH>
unsigned mpz_manager<SYNCH>::next_power_of_two(mpz const & a) {
    if (is_nonpos(a))
        return 0;
    if (is_one(a))
        return 0;
    unsigned shift;
    if (is_power_of_two(a, shift))
        return shift;
    else
        return log2(a) + 1;
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::is_perfect_square(mpz const & a, mpz & root) {
    if (is_neg(a))
        return false;
    set(root, 0);
    if (is_zero(a)) {
        return true;
    }
    if (is_one(a)) {
        set(root, 1);
        return true;
    }
    // current contract is that root is set to an approximation within +1/-1 of actional root.
    // x^2 mod 16 in { 9, 1, 4, 0 }
    auto mod16 = get_least_significant(a) & 0xF;
    if (mod16 != 0 && mod16 != 1 && mod16 != 4 && mod16 != 9)
        return false;

    mpz lo, hi, mid, sq_lo, sq_mid;
    set(lo, 1);
    set(hi, a);
    set(sq_lo, 1);

    bool result = false;
    // lo*lo <= *this < hi*hi

    // first find small interval lo*lo <= a <<= hi*hi
    while (true) {
        assert(lt(lo, hi));

        if (eq(sq_lo, a)) {
            set(root, lo);
            result = true;
            break;
        }
        mpz& tmp = mid;
        mul(lo, mpz(2), tmp);
        if (gt(tmp, hi))
            break;
        mul(tmp, tmp, sq_mid);
        if (gt(sq_mid, a)) {
            set(hi, tmp);
            break;
        }
        set(lo, tmp);
        set(sq_lo, sq_mid);
    }

    while (!result) {
        assert(lt(lo, hi));

        if (eq(sq_lo, a)) {
            set(root, lo);
            result = true;
            break;
        }

        mpz & tmp = mid;

        add(lo, mpz(1), tmp);
        if (eq(tmp, hi)) {
            set(root, hi);
            result = false;
            break;
        }

        add(hi, lo, tmp);
        div(tmp, mpz(2), mid);

        assert(lt(lo, mid) && lt(mid, hi));

        mul(mid, mid, sq_mid);

        if (gt(sq_mid, a)) {
            set(hi, mid);
        }
        else {
            set(lo, mid);
            set(sq_lo, sq_mid);
        }
    }
    del(lo);
    del(hi);
    del(mid);
    del(sq_lo);
    del(sq_mid);
    return result;
}


static unsigned div_l(unsigned k, unsigned n) {
    return k/n;
}

static unsigned div_u(unsigned k, unsigned n) {
    return k%n == 0 ? k/n : k/n + 1;
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::root(mpz & a, unsigned n) {
    assert(n % 2 != 0 || is_nonneg(a));
    if (is_zero(a)) {
        return true; // precise
    }

    // Initial approximation
    //
    // We have that:
    // a >  0 ->    2^{log2(a)}     <= a <= 2^{(log2(a) + 1)}
    // a <  0 ->   -2^{log2(a) + 1} <= a <= -2^{log2(a)}
    //
    // Thus
    // a >  0 ->    2^{div_l(log2(a), n)}     <= a^{1/n} <=  2^{div_u(log2(a) + 1, n)}
    // a <  0 ->   -2^{div_u(log2(a) + 1, n)} <= a^{1/n} <= -2^{div_l(log2(a), n)}
    //
    mpz lower;
    mpz upper;
    mpz mid;
    mpz mid_n;

    if (is_pos(a)) {
        unsigned k = log2(a);
        power(mpz(2), div_l(k, n),     lower);
        power(mpz(2), div_u(k + 1, n), upper);
    }
    else {
        unsigned k = mlog2(a);
        power(mpz(2), div_u(k + 1, n), lower);
        power(mpz(2), div_l(k, n),     upper);
        neg(lower);
        neg(upper);
    }

    bool result;
    assert(le(lower, upper));
    if (eq(lower, upper)) {
        swap(a, lower);
        result = true;
    }
    else {
        // Refine using bisection. TODO: use Newton's method if this is a bottleneck
        while (true) {
            add(upper, lower, mid);
            machine_div2k(mid, 1);
            power(mid, n, mid_n);
            if (eq(mid_n, a)) {
                swap(a, mid);
                result = true;
                break;
            }
            if (eq(mid, lower) || eq(mid, upper)) {
                swap(a, upper);
                result = false;
                break;
            }
            if (lt(mid_n, a)) {
                // new lower bound
                swap(mid, lower);
            }
            else {
                assert(lt(a, mid_n));
                // new upper bound
                swap(mid, upper);
            }
        }
    }
    del(lower);
    del(upper);
    del(mid);
    del(mid_n);
    return result;
}

template<bool SYNCH>
digit_t mpz_manager<SYNCH>::get_least_significant(mpz const& a) {
    assert(!is_neg(a));
    if (is_small(a))
        return std::abs(a.m_val);
    mpz_cell* cell_a = a.m_ptr;
    unsigned sz = cell_a->m_size;
    if (sz == 0)
        return 0;
    return cell_a->m_digits[0];
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::decompose(mpz const & a, std::vector<digit_t> & digits) {
    digits.clear();
    if (is_small(a)) {
        if (a.m_val < 0) {
            digits.push_back(-a.m_val);
            return true;
        }
        else {
            digits.push_back(a.m_val);
            return false;
        }
    }
    else {
        mpz_cell * cell_a = a.m_ptr;
        unsigned sz = cell_a->m_size;
        for (unsigned i = 0; i < sz; ++i) {
            digits.push_back(cell_a->m_digits[i]);
        }
        return a.m_val < 0;
    }
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::get_bit(mpz const & a, unsigned index) {
    if (is_small(a)) {
        assert(a.m_val >= 0);
        if (index >= 8*sizeof(digit_t))
            return false;
        return 0 != (a.m_val & (1ull << (digit_t)index));
    }
    unsigned i = index / (sizeof(digit_t)*8);
    unsigned o = index % (sizeof(digit_t)*8);

    mpz_cell * cell_a = a.m_ptr;
    unsigned sz = cell_a->m_size;
    if (sz*sizeof(digit_t)*8 <= index)
        return false;
    return 0 != (cell_a->m_digits[i] & (1ull << (digit_t)o));
}

template<bool SYNCH>
bool mpz_manager<SYNCH>::divides(mpz const & a, mpz const & b) {
    _scoped_numeral<mpz_manager<SYNCH> > tmp(*this);
    bool r;
    if (is_zero(a)) {
        // I assume 0 | 0.
        // Remark a|b is a shorthand for (exists x. a x = b)
        // If b is zero, any x will do. If b != 0, then a does not divide b
        r = is_zero(b);
    }
    else {
        rem(b, a, tmp);
        r = is_zero(tmp);
    }
    return r;
}

#ifndef SINGLE_THREAD
template class mpz_manager<true>;
#endif
template class mpz_manager<false>;
