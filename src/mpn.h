/**
 * @file mpn.h
 * @brief Multi-precision non-negative integer arithmetic library.
 *
 * This header declares the @c mpn_manager class, which provides
 * arbitrary-precision (bignum) arithmetic operations on non-negative
 * integers stored as arrays of @c mpn_digit values.
 *
 * All operations follow algorithms described in Donald E. Knuth,
 * "The Art of Computer Programming", Vol. 2, Section 4.3.1–4.3.3.
 *
 * @par Digit representation
 *   Numbers are stored as arrays of @c mpn_digit (unsigned int, 32-bit)
 *   in **little-endian order**: index 0 holds the least significant
 *   digit, index @c lng-1 holds the most significant digit.
 *   A number of length @c lng occupies @c lng consecutive elements
 *   of an @c mpn_digit array.  The most significant digit
 *   (at index @c lng-1) is assumed to be non-zero unless the number
 *   itself is zero.
 *
 * @par Conventions
 *   - Lengths are expressed in number of @c mpn_digit elements.
 *   - Caller is responsible for allocating sufficiently large output
 *     buffers.
 *   - Methods do not allocate memory for the primary operands;
 *     internal scratch space is managed via @c std::vector.
 */

#pragma once
#include <vector>
#include <string>

/**
 * @brief A single digit of a multi-precision number.
 *
 * Each digit holds @c sizeof(unsigned int) * 8 = 32 bits of the number.
 * The full multi-precision integer is represented as an array of
 * these digits in little-endian order (least significant digit first).
 */
typedef unsigned int mpn_digit;

/**
 * @class mpn_manager
 * @brief Manager class for multi-precision non-negative integer operations.
 *
 * The @c mpn_manager class implements comparison, addition, subtraction,
 * multiplication, division, and decimal string conversion for
 * arbitrary-precision non-negative integers represented as arrays of
 * @c mpn_digit.
 *
 * @par Usage example
 * @code
 *   mpn_manager mpn;
 *
 *   // Multiply 0xFFFFFFFF by 0xFFFFFFFF
 *   mpn_digit a[] = {0xFFFFFFFF};
 *   mpn_digit b[] = {0xFFFFFFFF};
 *   mpn_digit c[2];
 *   mpn.mul(a, 1, b, 1, c);
 *
 *   // c now holds the 64-bit result in two 32-bit digits
 *   std::string s = mpn.to_string(c, 2);
 *   // s == "18446744065119617025"
 * @endcode
 *
 * @par Algorithms
 *   - Addition:   Knuth's Algorithm A (TAOCP Vol. 2, 4.3.1)
 *   - Subtraction: Knuth's Algorithm S (TAOCP Vol. 2, 4.3.1)
 *   - Multiplication: Knuth's Algorithm M (TAOCP Vol. 2, 4.3.1)
 *   - Division:   Knuth's Algorithm D (TAOCP Vol. 2, 4.3.1)
 *
 * @note All public methods are @c const except @c div(), which may
 *       use internal mutable state in future versions.
 * @note None of the methods are thread-safe when operating on
 *       shared mutable buffers.
 */
class mpn_manager {

public:
    /**
     * @brief Compare two multi-precision non-negative integers.
     *
     * Performs an element-wise comparison from the most significant
     * digit to the least significant digit, treating missing digits
     * (beyond the declared length) as zero.
     *
     * @param a    Pointer to the first operand (little-endian digits).
     * @param lnga Number of digits in @p a.
     * @param b    Pointer to the second operand (little-endian digits).
     * @param lngb Number of digits in @p b.
     * @return  1  if @p a > @p b,
     * @return -1  if @p a < @p b,
     * @return  0  if @p a == @p b.
     *
     * @pre @p a points to at least @p lnga valid @c mpn_digit elements.
     * @pre @p b points to at least @p lngb valid @c mpn_digit elements.
     */
    int compare(mpn_digit const * a, unsigned lnga,
                mpn_digit const * b, unsigned lngb) const;

    /**
     * @brief Add two multi-precision non-negative integers.
     *
     * Computes @p c = @p a + @p b using Knuth's Algorithm A.
     * The result buffer must be at least @c max(lnga, lngb) + 1 digits
     * long to accommodate a possible carry into the most significant
     * position.  Leading zero digits are stripped from the result
     * (the output length is adjusted so that the most significant
     * digit is non-zero, except when the result is zero itself,
     * in which case the length is 1).
     *
     * @param a         Pointer to the first addend (little-endian digits).
     * @param lnga       Number of digits in @p a.
     * @param b         Pointer to the second addend (little-endian digits).
     * @param lngb       Number of digits in @p b.
     * @param c         Output buffer for the sum.  Must have room for
     *                  at least @c max(lnga, lngb) + 1 digits.
     * @param lngc_alloc Allocated size of @p c (in digits).  Must equal
     *                   @c max(lnga, lngb) + 1.
     * @param plngc     On return, points to the actual number of digits
     *                  written into @p c (after stripping leading zeros).
     *
     * @pre @p lngc_alloc == @c max(lnga, lngb) + 1 and @c max(lnga, lngb) > 0.
     * @pre @p plngc is not @c nullptr.
     *
     * @par Algorithm
     *   Knuth's Algorithm A (TAOCP Vol. 2, 4.3.1, p. 265).
     *   Each digit pair is summed with carry propagation;
     *   the final carry digit is stored at position @c len.
     */
    void add(mpn_digit const * a, unsigned lnga,
             mpn_digit const * b, unsigned lngb,
             mpn_digit *c, unsigned lngc_alloc,
             unsigned * plngc) const;

    /**
     * @brief Subtract two multi-precision non-negative integers.
     *
     * Computes @p c = @p a - @p b using Knuth's Algorithm S.
     * If @p a < @p b, the result wraps (unsigned subtraction) and a
     * borrow flag of 1 is returned through @p pborrow; otherwise the
     * borrow is 0.
     *
     * @param a      Pointer to the minuend (little-endian digits).
     * @param lnga    Number of digits in @p a.
     * @param b      Pointer to the subtrahend (little-endian digits).
     * @param lngb    Number of digits in @p b.
     * @param c      Output buffer for the difference.  Must have room
     *               for at least @c max(lnga, lngb) digits.
     * @param pborrow On return, contains the final borrow value:
     *               0 if @p a >= @p b (no underflow),
     *               1 if @p a < @p b (underflow occurred).
     *
     * @pre @p c points to at least @c max(lnga, lngb) valid @c mpn_digit elements.
     * @pre @p pborrow is not @c nullptr.
     *
     * @par Algorithm
     *   Knuth's Algorithm S (TAOCP Vol. 2, 4.3.1, p. 272).
     *   Each digit pair is subtracted with borrow propagation
     *   analogous to carry in addition.
     */
    void sub(mpn_digit const * a, unsigned lnga,
             mpn_digit const * b, unsigned lngb,
             mpn_digit * c, mpn_digit * pborrow) const;

    /**
     * @brief Multiply two multi-precision non-negative integers.
     *
     * Computes @p c = @p a * @p b using Knuth's Algorithm M
     * (schoolbook multiplication).  The result buffer must be at
     * least @c lnga + lngb digits long.
     *
     * @param a    Pointer to the first factor (little-endian digits).
     * @param lnga  Number of digits in @p a.
     * @param b    Pointer to the second factor (little-endian digits).
     * @param lngb  Number of digits in @p b.
     * @param c    Output buffer for the product.  Must have room for
     *             at least @c lnga + lngb digits.
     *
     * @pre @p c points to at least @c lnga + lngb valid @c mpn_digit elements.
     * @pre @p a and @p b are not @c nullptr (unless their length is 0).
     *
     * @par Algorithm
     *   Knuth's Algorithm M (TAOCP Vol. 2, 4.3.1, p. 268).
     *   This is the O(n*m) schoolbook algorithm.  Each digit of @p b
     *   is multiplied by every digit of @p a with carry, accumulating
     *   the partial products in @p c.  The branch for zero digits
     *   of @p b is an optional optimisation.
     *
     * @note A more efficient algorithm (e.g., Karatsuba) could be
     *       implemented for large operands; see Knuth 4.3.3.
     */
    void mul(mpn_digit const * a, unsigned lnga,
             mpn_digit const * b, unsigned lngb,
             mpn_digit * c) const;

    /**
     * @brief Divide two multi-precision non-negative integers, producing
     *        quotient and remainder.
     *
     * Computes @p quot = @p numer / @p denom and @p rem = @p numer % @p denom
     * using Knuth's Algorithm D for multi-digit divisors and a simplified
     * routine for single-digit divisors.
     *
     * Three code paths are taken:
     *   - If both operands are single-digit, a native CPU division is used.
     *   - If the numerator is strictly smaller than the denominator,
     *     the quotient is 0 and the remainder equals the numerator.
     *   - Otherwise, the operands are normalised, division is performed
     *     (Algorithm D for multi-digit denominators, or @c div_1 for
     *     single-digit denominators), and the remainder is un-normalised.
     *
     * @param numer Pointer to the numerator (little-endian digits).
     * @param lnum   Number of digits in @p numer.
     * @param denom Pointer to the denominator (little-endian digits).
     * @param lden   Number of digits in @p denom.
     * @param quot  Output buffer for the quotient.  Must have room for
     *              at least @c max(lnum - lden + 1, 1) digits.
     * @param rem   Output buffer for the remainder.  Must have room for
     *              at least @c lden digits.
     *
     * @pre @p lden > 0 (no division by zero).
     * @pre @p lnum > 0 (numerator is non-empty).
     * @p numer, @p denom, @p quot, @p rem are not @c nullptr.
     * @pre The most significant digit of @p denom (@p denom[lden-1]) is non-zero.
     *
     * @par Algorithm
     *   Knuth's Algorithm D (TAOCP Vol. 2, 4.3.1, p. 272).
     *   The denominator is first normalised by left-shifting so that its
     *   most significant digit has its highest bit set.  The same shift
     *   is applied to the numerator.  This guarantees that each trial
     *   quotient digit @c q_hat satisfies @c q_hat <= BASE, reducing
     *   correction steps.  After division, the remainder is shifted back
     *   (un-normalised).
     *
     * @note This method is not @c const because future versions may
     *       cache internal state.
     */
    void div(mpn_digit const * numer, unsigned lnum,
             mpn_digit const * denom, unsigned lden,
             mpn_digit * quot,
             mpn_digit * rem);

    /**
     * @brief Convert a multi-precision integer to a decimal string
     *        using a caller-provided buffer.
     *
     * Calls @c to_string(const, unsigned) const and copies the result
     * into @p buf, truncating to @p lbuf characters if necessary.
     * The output is **not** null-terminated if truncated.
     *
     * @param a    Pointer to the number (little-endian digits).
     * @param lng   Number of digits in @p a.
     * @param buf  Caller-allocated output buffer.
     * @param lbuf Size of @p buf in bytes (including space for NUL).
     * @return Pointer to @p buf (the same pointer passed in).
     *
     * @pre @p buf is not @c nullptr and @p lbuf > 0.
     *
     * @warning If @p lbuf is too small, the result is silently truncated.
     *          No NUL terminator is written in that case.
     */
    char * to_string(mpn_digit const * a, unsigned lng,
                     char * buf, unsigned lbuf) const;

    /**
     * @brief Convert a multi-precision integer to a decimal string.
     *
     * Converts the number @p a of length @p lng into a base-10 string
     * representation.  The conversion is performed by repeatedly dividing
     * the number by 10^9 (which fits in a single @c mpn_digit), extracting
     * 9 decimal digits per iteration, and assembling the groups from
     * most significant to least significant.
     *
     * @param a   Pointer to the number (little-endian digits).
     * @param lng  Number of digits in @p a.
     * @return A @c std::string containing the decimal representation.
     *         Returns "0" if @p lng is 0 or if all digits are zero.
     *
     * @pre @p a is not @c nullptr (unless @p lng is 0).
     *
     * @par Algorithm
     *   The number is divided by 10^9 in each iteration using the
     *   internal @c div_normalize / @c div_1 / @c div_unnormalize
     *   pipeline.  Each remainder forms a group of up to 9 decimal
     *   digits.  Groups are collected from least significant to
     *   most significant, then concatenated in reverse order.
     *   The most significant group is printed without leading zeros;
     *   subsequent groups are zero-padded to exactly 9 digits.
     */
    std::string to_string(mpn_digit const * a, unsigned lng) const;

private:
    /** @brief Internal scratch buffer type for multi-precision digits. */
    using mpn_sbuffer = std::vector<mpn_digit>;

    /**
     * @brief Normalise the numerator and denominator for division.
     *
     * Left-shifts both the numerator and the denominator by @c d bits
     * so that the most significant digit of the normalised denominator
     * has its highest bit set.  This is the normalisation step
     * described in Knuth's Algorithm D.
     *
     * @param numer    Pointer to the original numerator.
     * @param lnum      Number of digits in @p numer.
     * @param denom    Pointer to the original denominator.
     * @param lden      Number of digits in @p denom.
     * @param n_numer  On return, the normalised numerator
     *                 (length @c lnum + 1, with a possible extra
     *                 zero digit at the top).
     * @param n_denom   On return, the normalised denominator
     *                 (length @c lden).
     * @return The shift amount @c d (0 <= d < 32).
     *
     * @pre @p lden > 0 and @p denom[lden-1] != 0.
     */
    unsigned div_normalize(mpn_digit const * numer, unsigned lnum,
                         mpn_digit const * denom, unsigned lden,
                         mpn_sbuffer & n_numer,
                         mpn_sbuffer & n_denom) const;

    /**
     * @brief Reverse the normalisation shift applied to the remainder.
     *
     * Right-shifts the remainder (stored in the low-order digits of
     * @p numer) by @c d bits to undo the effect of @c div_normalize.
     *
     * @param numer  The normalised numerator buffer; the low
     *               @c denom.size() digits hold the remainder.
     * @param denom  The normalised denominator (used only for its
     *               size to determine how many remainder digits to
     *               extract).
     * @param d      The shift amount returned by @c div_normalize.
     * @param rem    Output buffer for the un-normalised remainder.
     *               Must have room for at least @c denom.size() digits.
     *
     * @pre @p d < @c DIGIT_BITS (i.e. d < 32 for 32-bit digits).
     * @pre @p denom.size() > 0.
     */
    void div_unnormalize(mpn_sbuffer & numer, mpn_sbuffer & denom,
                         unsigned d, mpn_digit * rem) const;

    /**
     * @brief Divide a multi-precision number by a single digit.
     *
     * Performs long division of @p numer by the single-digit
     * @p denom, writing the quotient into @p quot.  The numerator
     * is consumed (modified in place) during the computation.
     *
     * @param numer  The normalised numerator (modified in place;
     *               the low digits are overwritten with intermediate
     *               remainders during the computation).
     * @param denom  A single normalised denominator digit.
     * @param quot   Output buffer for the quotient.  Must have room
     *               for at least @c numer.size() - 1 digits.
     *
     * @pre @p denom != 0.
     * @pre @p numer.size() >= 1.
     *
     * @par Algorithm
     *   Processes digits from most significant to least significant.
     *   At each step, a two-digit window is formed from the current
     *   and previous digits, divided by @p denom to obtain the
     *   trial quotient digit and remainder.  Borrow correction is
     *   applied when the trial quotient overshoots.
     */
    void div_1(mpn_sbuffer & numer, mpn_digit denom,
               mpn_digit * quot) const;

    /**
     * @brief Divide a multi-precision number by a multi-precision number
     *        (Knuth's Algorithm D core loop).
     *
     * Performs the main loop of Knuth's Algorithm D for divisors of
     * length >= 2.  For each quotient digit, a trial quotient @c q_hat
     * is estimated from the two most significant digits of the current
     * numerator window and the most significant digit of the
     * denominator, then refined using the next denominator digit.
     * The trial product @c q_hat * denom is subtracted from the
     * numerator window; if a borrow occurs, the quotient digit is
     * decremented and the denominator is added back.
     *
     * @param numer  The normalised numerator (modified in place).
     *               Must have length @c m + n where @c m is the
     *               number of quotient digits and @c n is the
     *               denominator length.
     * @param denom  The normalised denominator (length @c n >= 2).
     * @param quot   Output buffer for the quotient (@c m digits).
     * @param rem    Output buffer for the remainder (@c n digits).
     * @param ms     Scratch buffer for the trial product
     *               (@c q_hat * denom).  Resized to @c n + 1 internally.
     * @param ab     Scratch buffer for the add-back correction.
     *               Resized to @c n + 2 internally.
     *
     * @pre @p denom.size() > 1.
     * @pre @p numer.size() == (@p numer.size() - @p denom.size()) + @p denom.size().
     *
     * @par Algorithm
     *   Knuth's Algorithm D (TAOCP Vol. 2, 4.3.1, p. 272–276).
     *   The trial quotient digit @c q_hat is computed as
     *   @c (numer[j+n] * BASE + numer[j+n-1]) / denom[n-1],
     *   then corrected by testing whether
     *   @c q_hat * denom[n-2] > r_hat * BASE + numer[j+n-2].
     *   If so, @c q_hat is decremented.  The product @c q_hat * denom
     *   is then subtracted from the numerator window; a final
     *   correction adds back @p denom if the subtraction borrows.
     */
    void div_n(mpn_sbuffer & numer, mpn_sbuffer const & denom,
               mpn_digit * quot, mpn_digit * rem,
               mpn_sbuffer & ms, mpn_sbuffer & ab) const;
};
