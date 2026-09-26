/**
 * @file mpn.h
 * @brief Multi-precision non-negative integer arithmetic library.
 *
 * This header declares the class, which provides
 * arbitrary-precision arithmetic operations on non-negative
 * integers stored as arrays of numbers.
 *
 * All operations follow algorithms described in Donald E. Knuth,
 * "The Art of Computer Programming", Vol. 2, Section 4.3.1–4.3.3.
 *
 * @par Digit representation
 *   Numbers are stored as arrays of unsigned int - 8, 16 or 32 bits
 *   in **little-endian order**: index 0 holds the least significant
 *   digit, index @c lng-1 holds the most significant digit.
 *   A number of length @c lng occupies @c lng consecutive elements
 *   of an @c MpnDigit array.  The most significant digit
 *   (at index @c lng-1) is assumed to be non-zero unless the number
 *   itself is zero.
 *
 * @par Conventions
 *   - Lengths are expressed in number of @c MpnDigit elements.
 *   - Caller is responsible for allocating sufficiently large output
 *     buffers.
 *   - Methods do not allocate memory for the primary operands;
 *     internal scratch space is managed via @c std::vector.
 */
#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace mpn_detail {

template<class T>
concept MpnDigit =
    std::is_integral_v<T> &&
    std::is_unsigned_v<T> &&
    !std::is_same_v<T, bool>;

// The arithmetic algorithms need a type that can hold two digits plus carry.
template<class T>
struct wide_type;

template<MpnDigit T>
requires (std::numeric_limits<T>::digits <= 8)
struct wide_type<T> { using type = std::uint16_t; };

template<MpnDigit T>
requires (std::numeric_limits<T>::digits > 8 && std::numeric_limits<T>::digits <= 16)
struct wide_type<T> { using type = std::uint32_t; };

template<MpnDigit T>
requires (std::numeric_limits<T>::digits > 16 && std::numeric_limits<T>::digits <= 32)
struct wide_type<T> { using type = std::uint64_t; };

template<MpnDigit T>
using wide_type_t = typename wide_type<T>::type;

template<MpnDigit T>
constexpr std::size_t digit_bits_v = std::numeric_limits<T>::digits;

template<MpnDigit T>
constexpr wide_type_t<T> base_v = wide_type_t<T>(1) << digit_bits_v<T>;

template<MpnDigit T>
constexpr T top_bit_v = T(1) << (digit_bits_v<T> - 1);

template<MpnDigit T>
constexpr T first_bits(std::size_t n, T x) noexcept {
    return n == 0 ? T(0) : (n >= digit_bits_v<T> ? x : T(x >> (digit_bits_v<T> - n)));
}

template<MpnDigit T>
constexpr T last_bits(std::size_t n, T x) noexcept {
    if (n == 0) return T(0);
    if (n >= digit_bits_v<T>) return x;
    return T((x << (digit_bits_v<T> - n)) >> (digit_bits_v<T> - n));
}

// Largest power of 10 that fits in Digit, and its decimal width.
template<MpnDigit T>
struct decimal_traits {
    static constexpr T base = [] {
        T value = 1;
        while (value <= (std::numeric_limits<T>::max() / T(10)))
            value = T(value * T(10));
        return value;
    }();

    static constexpr std::size_t digits = [] {
        std::size_t n = 0;
        T value = base;
        while (value > T(1)) {
            value = T(value / T(10));
            ++n;
        }
        return n;
    }();
};

} // namespace mpn_detail

/**
 * @brief Multi-precision unsigned integer arithmetic over configurable digits.
 * @tparam Digit Unsigned integral digit type. Supported widths are 8, 16, 32 bits.
 */
template<mpn_detail::MpnDigit Digit = std::uint32_t>
class mpn_manager {
public:
    using digit_type = Digit;
    using wide_type = mpn_detail::wide_type_t<Digit>;
    using buffer_type = std::vector<Digit>;

    /// Number of value bits stored in one digit.
    static constexpr std::size_t digit_bits = mpn_detail::digit_bits_v<Digit>;

    /// Numeric base represented by one digit position: 2^digit_bits.
    static constexpr wide_type base = mpn_detail::base_v<Digit>;

    /**
     * @brief Compares two non-negative multi-precision integers.
     *
     * Each operand is represented as a little-endian array: element zero is
     * the least-significant digit. Leading zero digits are allowed and do not
     * affect the result.
     *
     * @param a Pointer to the first operand, or `nullptr` when @p lnga is zero.
     * @param lnga Number of digits in the first operand.
     * @param b Pointer to the second operand, or `nullptr` when @p lngb is zero.
     * @param lngb Number of digits in the second operand.
     * @return `-1` if `a < b`, `0` if `a == b`, or `1` if `a > b`.
     * @throws std::invalid_argument If a non-empty operand has a null pointer.
     * @note The function does not modify either operand.
     */
    [[nodiscard]] int compare(Digit const* a, std::size_t lnga,
                              Digit const* b, std::size_t lngb) const;

    /**
     * @brief Adds two non-negative multi-precision integers.
     *
     * The result is written to @p c in little-endian form. The result length
     * excludes insignificant leading zero digits, except that zero itself is
     * represented by one digit. The output buffer must provide room for one
     * possible carry digit in addition to `max(lnga, lngb)` digits.
     *
     * @param a First operand.
     * @param lnga Number of digits in @p a.
     * @param b Second operand.
     * @param lngb Number of digits in @p b.
     * @param c Output buffer. It must not overlap either input buffer.
     * @param lngc_alloc Number of allocated digits in @p c.
     * @param plngc Receives the number of significant digits written to @p c.
     * @throws std::invalid_argument If an output length pointer is null, an
     *         operand pointer is null while its length is non-zero, or both
     *         operands are empty.
     * @throws std::out_of_range If @p c is too small.
     * @throws std::length_error If the required result size cannot be represented
     *         by `std::size_t`.
     */
    void add(Digit const* a, std::size_t lnga,
             Digit const* b, std::size_t lngb,
             Digit* c, std::size_t lngc_alloc,
             std::size_t* plngc) const;

    /**
     * @brief Subtracts one non-negative multi-precision integer from another.
     *
     * The operation computes `a - b` modulo `base^N`, where
     * `N = max(lnga, lngb)`. If `a < b`, @p pborrow receives one; otherwise it
     * receives zero. The output is always exactly `N` digits.
     *
     * @param a Minuend.
     * @param lnga Number of digits in @p a.
     * @param b Subtrahend.
     * @param lngb Number of digits in @p b.
     * @param c Output buffer with at least `max(lnga, lngb)` digits. It must not
     *        overlap either input buffer.
     * @param pborrow Receives the final borrow (`0` or `1`).
     * @throws std::invalid_argument If an output pointer is null, an operand
     *         pointer is null while its length is non-zero, or both operands are empty.
     */
    void sub(Digit const* a, std::size_t lnga,
             Digit const* b, std::size_t lngb,
             Digit* c, Digit* pborrow) const;

    /**
     * @brief Multiplies two non-negative multi-precision integers.
     *
     * The product is written in little-endian form to @p c. The caller must
     * provide at least `lnga + lngb` digits of storage. The output buffer must
     * not overlap either input buffer.
     *
     * @param a First factor.
     * @param lnga Number of digits in @p a.
     * @param b Second factor.
     * @param lngb Number of digits in @p b.
     * @param c Output buffer with at least `lnga + lngb` digits.
     * @throws std::invalid_argument If an operand is empty or an input/output
     *         pointer is null.
     * @throws std::length_error If `lnga + lngb` overflows `std::size_t`.
     */
    void mul(Digit const* a, std::size_t lnga,
             Digit const* b, std::size_t lngb,
             Digit* c) const;

    /**
     * @brief Divides a non-negative multi-precision integer by another.
     *
     * Computes `numer = quot * denom + rem`. Both operands are represented in
     * little-endian order. The denominator must be non-zero and must not have
     * leading zero digits. The quotient and remainder buffers must be large
     * enough for the result and must not overlap the input buffers.
     *
     * @param numer Numerator.
     * @param lnum Number of numerator digits.
     * @param denom Non-zero denominator.
     * @param lden Number of denominator digits.
     * @param quot Output quotient buffer. The caller must provide enough storage
     *        for the quotient; at most `max(1, lnum - lden + 1)` digits are used.
     * @param rem Output remainder buffer with at least `lden` digits.
     * @throws std::invalid_argument If an operand is empty, an input/output
     *         pointer is null, or the denominator has a leading zero digit.
     * @note This low-level overload does not receive output capacities, so the
     *       caller is responsible for providing sufficiently large output buffers.
     */
    void div(Digit const* numer, std::size_t lnum,
             Digit const* denom, std::size_t lden,
             Digit* quot, Digit* rem) const;

    /**
     * @brief Converts a multi-precision integer to its decimal representation.
     *
     * Digits are interpreted as a little-endian unsigned integer. An empty
     * input represents zero, as does an input consisting entirely of zero digits.
     *
     * @param a Pointer to the digit array, or `nullptr` when @p lng is zero.
     * @param lng Number of digits in @p a.
     * @return Decimal representation without a sign or leading zeroes.
     * @throws std::invalid_argument If @p a is null while @p lng is non-zero.
     */
    [[nodiscard]] std::string to_string(Digit const* a, std::size_t lng) const;

    /**
     * @brief Compares two multi-precision integers represented by spans.
     * @param a First little-endian operand.
     * @param b Second little-endian operand.
     * @return `-1`, `0`, or `1` according to the numerical ordering.
     * @throws std::invalid_argument If a non-empty span has a null data pointer.
     */
    [[nodiscard]] int compare(std::span<const Digit> a,
                              std::span<const Digit> b) const;

    /**
     * @brief Adds two integers represented by spans.
     * @param a First little-endian operand.
     * @param b Second little-endian operand.
     * @param c Output span; it must contain at least `max(a.size(), b.size()) + 1` digits.
     * @param result_size Receives the number of significant result digits.
     * @throws std::invalid_argument If both operands are empty.
     * @throws std::out_of_range If @p c is too small.
     * @throws std::length_error If the required output size cannot be represented.
     */
    void add(std::span<const Digit> a, std::span<const Digit> b,
             std::span<Digit> c, std::size_t& result_size) const;

    /**
     * @brief Subtracts the second span from the first span.
     * @param a Minuend in little-endian order.
     * @param b Subtrahend in little-endian order.
     * @param c Output span with at least `max(a.size(), b.size())` elements.
     * @param borrow Receives the final borrow (`0` or `1`).
     * @throws std::invalid_argument If both operands are empty.
     * @throws std::out_of_range If @p c is too small.
     */
    void sub(std::span<const Digit> a, std::span<const Digit> b,
             std::span<Digit> c, Digit& borrow) const;

    /**
     * @brief Multiplies two integers represented by spans.
     * @param a First factor.
     * @param b Second factor.
     * @param c Output span with at least `a.size() + b.size()` elements.
     * @throws std::invalid_argument If either operand is empty.
     * @throws std::out_of_range If @p c is too small.
     * @throws std::length_error If the required output size overflows `std::size_t`.
     */
    void mul(std::span<const Digit> a, std::span<const Digit> b,
             std::span<Digit> c) const;

    /**
     * @brief Divides two integers represented by spans.
     * @param numer Numerator in little-endian order.
     * @param denom Non-zero denominator in little-endian order.
     * @param quot Output quotient span. It must be large enough for the quotient.
     * @param rem Output remainder span with at least `denom.size()` elements.
     * @throws std::invalid_argument If either operand is empty or the denominator
     *         contains a leading zero digit.
     * @throws std::out_of_range If an output span is too small.
     */
    void div(std::span<const Digit> numer, std::span<const Digit> denom,
             std::span<Digit> quot, std::span<Digit> rem) const;

    /**
     * @brief Converts a span of little-endian digits to decimal text.
     * @param a Digit span; an empty span represents zero.
     * @return Decimal representation without a sign or leading zeroes.
     * @throws std::invalid_argument If the span is non-empty but does not contain valid data.
     */
    [[nodiscard]] std::string to_string(std::span<const Digit> a) const;

private:
    using mpn_sbuffer = std::vector<Digit>;

    /**
     * @brief Normalizes numerator and denominator for Knuth's division algorithm.
     * @param numer Input numerator in little-endian order.
     * @param lnum Number of numerator digits.
     * @param denom Input denominator in little-endian order.
     * @param lden Number of denominator digits.
     * @param n_numer Receives the normalized numerator and one extra high digit.
     * @param n_denom Receives the normalized denominator.
     * @return Number of left-shift positions applied during normalization.
     * @throws std::logic_error If the normalization shift would be outside the digit width.
     */
    std::size_t div_normalize(Digit const* numer, std::size_t lnum,
                              Digit const* denom, std::size_t lden,
                              mpn_sbuffer& n_numer,
                              mpn_sbuffer& n_denom) const;

    /**
     * @brief Restores the original scale of a remainder after division.
     * @param numer Normalized numerator/remainder workspace.
     * @param denom Normalized denominator workspace.
     * @param d Normalization shift returned by div_normalize().
     * @param rem Destination remainder buffer.
     */
    void div_unnormalize(mpn_sbuffer const& numer,
                         mpn_sbuffer const& denom,
                         std::size_t d, Digit* rem) const;

    /**
     * @brief Divides a mutable digit buffer by a single normalized digit.
     * @param numer Numerator workspace; it is replaced by the quotient/remainder workspace.
     * @param denom Single-digit divisor.
     * @param quot Destination quotient buffer.
     */
    void div_1(mpn_sbuffer& numer, Digit denom, Digit* quot) const;

    /**
     * @brief Performs the multi-digit step of Knuth's Algorithm D.
     * @param numer Mutable normalized numerator workspace.
     * @param denom Normalized multi-digit denominator.
     * @param quot Destination quotient buffer.
     * @param ms Scratch buffer for multiplication of the divisor by a trial quotient digit.
     * @param ab Scratch buffer used when a correction/add-back step is required.
     * @throws std::invalid_argument If the denominator contains fewer than two digits.
     */
    void div_n(mpn_sbuffer& numer, mpn_sbuffer const& denom,
               Digit* quot, mpn_sbuffer& ms,
               mpn_sbuffer& ab) const;
};

template<mpn_detail::MpnDigit Digit>
int mpn_manager<Digit>::compare(Digit const* a, std::size_t lnga,
                                Digit const* b, std::size_t lngb) const {
    if ((lnga != 0 && a == nullptr) || (lngb != 0 && b == nullptr))
        throw std::invalid_argument("compare: non-empty operand has a null pointer");
    std::size_t j = std::max(lnga, lngb);
    for (; j-- > 0;) {
        const Digit u = (j < lnga) ? a[j] : Digit(0);
        const Digit v = (j < lngb) ? b[j] : Digit(0);
        if (u > v) return 1;
        if (u < v) return -1;
    }
    return 0;
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::add(Digit const* a, std::size_t lnga,
                             Digit const* b, std::size_t lngb,
                             Digit* c, std::size_t lngc_alloc,
                             std::size_t* plngc) const {
    const std::size_t len = std::max(lnga, lngb);
    if (plngc == nullptr) throw std::invalid_argument("add: plngc must not be null");
    if (len == 0) throw std::invalid_argument("add: both operands are empty");
    if (len == std::numeric_limits<std::size_t>::max())
        throw std::length_error("add: result size overflows std::size_t");
    if (c == nullptr) throw std::invalid_argument("add: output buffer must not be null");
    if (lngc_alloc < len + 1) throw std::out_of_range("add: output buffer is too small");
    if (lnga > 0 && a == nullptr) throw std::invalid_argument("add: first operand is null");
    if (lngb > 0 && b == nullptr) throw std::invalid_argument("add: second operand is null");

    Digit carry = 0;
    for (std::size_t j = 0; j < len; ++j) {
        const Digit u = (j < lnga) ? a[j] : Digit(0);
        const Digit v = (j < lngb) ? b[j] : Digit(0);
        const wide_type t = static_cast<wide_type>(wide_type(u) + wide_type(v)) + wide_type(carry);
        c[j] = static_cast<Digit>(t);
        carry = static_cast<Digit>(t >> digit_bits);
    }
    c[len] = carry;

    std::size_t& out_len = *plngc;
    out_len = len + 1;
    while (out_len > 1 && c[out_len - 1] == 0) --out_len;
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::sub(Digit const* a, std::size_t lnga,
                             Digit const* b, std::size_t lngb,
                             Digit* c, Digit* pborrow) const {
    if (pborrow == nullptr) throw std::invalid_argument("sub: pborrow must not be null");
    if (std::max(lnga, lngb) == 0) throw std::invalid_argument("sub: both operands are empty");
    if (c == nullptr) throw std::invalid_argument("sub: output buffer must not be null");
    if (lnga > 0 && a == nullptr) throw std::invalid_argument("sub: first operand is null");
    if (lngb > 0 && b == nullptr) throw std::invalid_argument("sub: second operand is null");
    const std::size_t len = std::max(lnga, lngb);
    Digit borrow = 0;

    for (std::size_t j = 0; j < len; ++j) {
        const Digit u = (j < lnga) ? a[j] : Digit(0);
        const Digit v = (j < lngb) ? b[j] : Digit(0);
        const wide_type subtrahend = wide_type(v) + wide_type(borrow);
        c[j] = static_cast<Digit>(wide_type(u) - subtrahend);
        borrow = (wide_type(u) < subtrahend) ? Digit(1) : Digit(0);
    }
    *pborrow = borrow;
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::mul(Digit const* a, std::size_t lnga,
                             Digit const* b, std::size_t lngb,
                             Digit* c) const {
    if (lnga == 0 || lngb == 0) throw std::invalid_argument("mul: operands must not be empty");
    if (lnga > std::numeric_limits<std::size_t>::max() - lngb)
        throw std::length_error("mul: result size overflows std::size_t");
    if (c == nullptr) throw std::invalid_argument("mul: output buffer must not be null");
    if (a == nullptr || b == nullptr) throw std::invalid_argument("mul: operands must not be null");
    for (std::size_t i = 0; i < lnga + lngb; ++i) c[i] = 0;

    for (std::size_t j = 0; j < lngb; ++j) {
        const Digit v = b[j];
        if (v == 0) continue;

        Digit carry = 0;
        for (std::size_t i = 0; i < lnga; ++i) {
            const wide_type t = static_cast<wide_type>(wide_type(a[i]) * wide_type(v) +
                                wide_type(c[i + j])) + wide_type(carry);
            c[i + j] = static_cast<Digit>(t);
            carry = static_cast<Digit>(t >> digit_bits);
        }
        c[j + lnga] = carry;
    }
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::div(Digit const* numer, std::size_t lnum,
                             Digit const* denom, std::size_t lden,
                             Digit* quot, Digit* rem) const {
    if (lden == 0) throw std::invalid_argument("div: denominator must not be empty");
    if (lnum == 0) throw std::invalid_argument("div: numerator must not be empty");
    if (numer == nullptr || denom == nullptr) throw std::invalid_argument("div: input must not be null");
    if (quot == nullptr || rem == nullptr) throw std::invalid_argument("div: output must not be null");
    if (denom[lden - 1] == 0) throw std::invalid_argument("div: denominator has a leading zero");

    if (lnum == 1 && lden == 1) {
        quot[0] = static_cast<Digit>(numer[0] / denom[0]);
        rem[0] = static_cast<Digit>(numer[0] % denom[0]);
        return;
    }

    if (lnum < lden ||
        (lnum == lden && numer[lnum - 1] < denom[lden - 1])) {
        quot[0] = 0;
        for (std::size_t i = 0; i < lden; ++i)
            rem[i] = (i < lnum) ? numer[i] : Digit(0);
        return;
    }

    mpn_sbuffer u, v, t_ms, t_ab;
    const std::size_t d = div_normalize(numer, lnum, denom, lden, u, v);
    if (lden == 1)
        div_1(u, v[0], quot);
    else
        div_n(u, v, quot, t_ms, t_ab);
    div_unnormalize(u, v, d, rem);
}

template<mpn_detail::MpnDigit Digit>
std::size_t mpn_manager<Digit>::div_normalize(
    Digit const* numer, std::size_t lnum,
    Digit const* denom, std::size_t lden,
    mpn_sbuffer& n_numer,
    mpn_sbuffer& n_denom) const {
    std::size_t d = 0;
    while (((denom[lden - 1] << d) & mpn_detail::top_bit_v<Digit>) == 0) ++d;
    if (d >= digit_bits) throw std::logic_error("div_normalize: invalid normalization shift");

    n_numer.resize(lnum + 1);
    n_denom.resize(lden);

    if (d == 0) {
        n_numer[lnum] = 0;
        for (std::size_t i = 0; i < lnum; ++i) n_numer[i] = numer[i];
        for (std::size_t i = 0; i < lden; ++i) n_denom[i] = denom[i];
        return 0;
    }

    n_numer[lnum] = mpn_detail::first_bits<Digit>(d, numer[lnum - 1]);
    for (std::size_t i = lnum - 1; i > 0; --i)
        n_numer[i] = static_cast<Digit>((wide_type(numer[i]) << d) |
                                         wide_type(mpn_detail::first_bits<Digit>(d, numer[i - 1])));
    n_numer[0] = static_cast<Digit>(wide_type(numer[0]) << d);

    for (std::size_t i = lden - 1; i > 0; --i)
        n_denom[i] = static_cast<Digit>((wide_type(denom[i]) << d) |
                                         wide_type(mpn_detail::first_bits<Digit>(d, denom[i - 1])));
    n_denom[0] = static_cast<Digit>(wide_type(denom[0]) << d);

    return d;
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::div_unnormalize(
    mpn_sbuffer const& numer, mpn_sbuffer const& denom,
    std::size_t d, Digit* rem) const {
    const std::size_t denom_size = denom.size();
    if (denom_size == 0) return;

    if (d == 0) {
        for (std::size_t i = 0; i < denom_size; ++i) rem[i] = numer[i];
        return;
    }

    const std::size_t limit = denom_size - 1;
    for (std::size_t i = 0; i < limit; ++i) {
        rem[i] = static_cast<Digit>(
            (wide_type(numer[i]) >> d) |
            (wide_type(mpn_detail::last_bits<Digit>(d, numer[i + 1])) << (digit_bits - d)));
    }
    rem[limit] = static_cast<Digit>(numer[limit] >> d);
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::div_1(mpn_sbuffer& numer, Digit denom, Digit* quot) const {
    wide_type q_hat, temp, ms;

    for (std::size_t j = static_cast<std::size_t>(numer.size() - 1); j > 0; --j) {
        temp = static_cast<wide_type>(
            (wide_type(numer[j]) << digit_bits) | wide_type(numer[j - 1]));
        q_hat = static_cast<wide_type>(temp / wide_type(denom));
        ms = static_cast<wide_type>(temp - q_hat * wide_type(denom));

        numer[j - 1] = static_cast<Digit>(ms);
        numer[j] = static_cast<Digit>(ms >> digit_bits);
        quot[j - 1] = static_cast<Digit>(q_hat);

    }
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::div_n(
    mpn_sbuffer& numer, mpn_sbuffer const& denom,
    Digit* quot, mpn_sbuffer& ms, mpn_sbuffer& ab) const {
    if (denom.size() <= 1) throw std::invalid_argument("div_n: denominator must contain at least two digits");

    const std::size_t m = static_cast<std::size_t>(numer.size() - denom.size());
    const std::size_t n = denom.size();

    ms.resize(n + 1);

    wide_type q_hat, temp, r_hat;
    Digit borrow;

    for (std::size_t j = m; j-- > 0;) {
        temp = static_cast<wide_type>(
            (wide_type(numer[j + n]) << digit_bits) | wide_type(numer[j + n - 1]));
        q_hat = static_cast<wide_type>(temp / wide_type(denom[n - 1]));
        r_hat = static_cast<wide_type>(temp % wide_type(denom[n - 1]));

        while (q_hat >= base ||
               q_hat * wide_type(denom[n - 2]) >
                   (r_hat << digit_bits) + wide_type(numer[j + n - 2])) {
            --q_hat;
            r_hat += wide_type(denom[n - 1]);
            if (r_hat >= base) break;
        }

        const Digit q_small = static_cast<Digit>(q_hat);
        mul(&q_small, 1, denom.data(), n, ms.data());
        sub(&numer[j], n + 1, ms.data(), n + 1, &numer[j], &borrow);
        quot[j] = q_small;

        if (borrow) {
            --quot[j];
            ab.resize(n + 2);
            std::size_t real_size = 0;
            add(denom.data(), n, &numer[j], n + 1,
                ab.data(), n + 2, &real_size);
            for (std::size_t i = 0; i < n + 1; ++i)
                numer[j + i] = ab[i];
        }
    }
}

template<mpn_detail::MpnDigit Digit>
std::string mpn_manager<Digit>::to_string(Digit const* a, std::size_t lng) const {
    if (lng > 0 && a == nullptr) throw std::invalid_argument("to_string: input must not be null");
    if (lng == 0) return "0";

    bool all_zero = true;
    for (std::size_t i = 0; i < lng; ++i) {
        if (a[i] != 0) { all_zero = false; break; }
    }
    if (all_zero) return "0";

    constexpr Digit decimal_base = mpn_detail::decimal_traits<Digit>::base;
    constexpr std::size_t decimal_digits = mpn_detail::decimal_traits<Digit>::digits;

    mpn_sbuffer temp(a, a + lng);
    mpn_sbuffer t_numer, t_denom;
    Digit rem = 0;
    std::vector<Digit> groups;

    while (!temp.empty() && (temp.size() > 1 || temp[0] != 0)) {
        const std::size_t d = div_normalize(temp.data(), temp.size(),
                                         &decimal_base, 1, t_numer, t_denom);
        div_1(t_numer, t_denom[0], temp.data());
        div_unnormalize(t_numer, t_denom, d, &rem);
        groups.push_back(rem);

        while (!temp.empty() && temp.back() == 0) temp.pop_back();
    }

    std::string result;
    result.reserve(groups.size() * decimal_digits);
    result += std::to_string(groups.back());

    for (std::size_t i = groups.size() - 1; i > 0; --i) {
        const std::string group = std::to_string(groups[i - 1]);
        result.append(decimal_digits - group.size(), '0');
        result += group;
    }
    return result;
}

template<mpn_detail::MpnDigit Digit>
int mpn_manager<Digit>::compare(std::span<const Digit> a, std::span<const Digit> b) const {
    if ((a.size() != 0 && a.data() == nullptr) || (b.size() != 0 && b.data() == nullptr))
        throw std::invalid_argument("compare(span): non-empty span has a null data pointer");
    return compare(a.data(), a.size(),
                   b.data(), b.size());
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::add(std::span<const Digit> a, std::span<const Digit> b,
                             std::span<Digit> c, std::size_t& result_size) const {
    const std::size_t max_size = std::max(a.size(), b.size());
    if (max_size == std::numeric_limits<std::size_t>::max())
        throw std::length_error("add(span): result size overflows std::size_t");
    const std::size_t required = max_size + 1;
    if (c.size() < required) throw std::out_of_range("add(span): output buffer is too small");
    std::size_t out = 0;
    add(a.data(), a.size(), b.data(), b.size(),
        c.data(), c.size(), &out);
    result_size = out;
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::sub(std::span<const Digit> a, std::span<const Digit> b,
                             std::span<Digit> c, Digit& borrow) const {
    const std::size_t required = std::max(a.size(), b.size());
    if (c.size() < required) throw std::out_of_range("sub(span): output buffer is too small");
    sub(a.data(), a.size(), b.data(), b.size(),
        c.data(), &borrow);
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::mul(std::span<const Digit> a, std::span<const Digit> b,
                             std::span<Digit> c) const {
    if (a.size() > std::numeric_limits<std::size_t>::max() - b.size())
        throw std::length_error("mul(span): result size overflows std::size_t");
    if (c.size() < a.size() + b.size())
        throw std::out_of_range("mul(span): output buffer is too small");
    mul(a.data(), a.size(), b.data(), b.size(), c.data());
}

template<mpn_detail::MpnDigit Digit>
void mpn_manager<Digit>::div(std::span<const Digit> numer, std::span<const Digit> denom,
                             std::span<Digit> quot, std::span<Digit> rem) const {
    const std::size_t required_quot =
        numer.size() >= denom.size() ? numer.size() - denom.size() + 1 : 1;
    if (quot.size() < required_quot)
        throw std::out_of_range("div(span): quotient buffer is too small");
    if (rem.size() < denom.size())
        throw std::out_of_range("div(span): remainder buffer is too small");
    div(numer.data(), numer.size(), denom.data(), denom.size(), quot.data(), rem.data());
}

template<mpn_detail::MpnDigit Digit>
std::string mpn_manager<Digit>::to_string(std::span<const Digit> a) const {
    return to_string(a.data(), a.size());
}
