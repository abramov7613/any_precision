# ez3i

**ez3i** (*extracted z3 integer*) is a standalone C++20 library for arbitrary-precision signed integer arithmetic. It was extracted from the [Z3 theorem prover](https://github.com/Z3Prover/z3) and packaged as an independent, dependency-free library.

The library provides exact, overflow-free arithmetic for integers of unlimited size. Small values (those fitting in a machine `int`) are handled inline without heap allocation; large values are stored as arrays of 32-bit digits in heap-allocated cells. This dual representation makes common-case operations fast while still supporting numbers of thousands of bits.

---

## Features

- **Signed big integers** with automatic small-number optimization — values that fit in a machine `int` are stored inline without heap allocation.
- **Full arithmetic**: addition, subtraction, multiplication, division, modulus, GCD, LCM, extended Euclid.
- **Bitwise operations** — AND, OR, XOR, NOT, bit extraction, power-of-two detection.
- **Utilities**: power-of-two, perfect-square test, integer logarithm, decimal string conversion.
- **Knuth's algorithms** — division follows Algorithm D from *The Art of Computer Programming*, Vol. 2, §4.3.
- **Modern C++20** — concepts, `std::span`, flexible digit widths (8/16/32-bit) at the low level.
- **Flexible value construction** — from `int`, `unsigned`, `int64_t`, `uint64_t`, C-strings, or raw digit arrays.
- **Conversion & display** — to `std::string`, to `double`, to SMT-LIB 2.0 format, to hex and binary streams.
- **Cross-platform** — builds on Windows (MSVC) and POSIX (GCC/Clang).
- **CMake** build system with static/shared library options and install targets.

---

## Architecture

The library is organised in two layers:

### `mpz_cell` — heap storage for large numbers

```cpp
struct mpz_cell {
    unsigned  m_size;      // number of significant digits currently stored
    unsigned  m_capacity;  // maximum number of digits the buffer can hold
    digit_t   m_digits[];  // flexible array member — little-endian digit storage
};
```

`mpz_cell` is a low-level POD structure that owns a variable-length array of `digit_t` (typedef for `unsigned int`, i.e. 32-bit on most platforms). It uses the C flexible array member idiom so that the digit buffer is allocated in the same block as the header — a single `operator new` call yields both the metadata and the storage. The `m_capacity` field records how many digits the buffer can hold; `m_size` records how many are currently significant (leading zero digits are trimmed by `normalize()`). The `mpz_manager` is solely responsible for allocating, resizing, and freeing `mpz_cell` objects.

**Why it exists.** When a number is too large to fit in a machine `int`, it must live on the heap. `mpz_cell` is that heap representation — a compact, self-contained block that avoids a separate pointer indirection for the digit array.

### `mpz` — the value handle

```cpp
class mpz {
    int        m_val;     // small value OR sign (+1 / -1) for large numbers
    unsigned   m_kind:1;  // mpz_small (0) or mpz_large (1)
    unsigned   m_owner:1; // mpz_self (0) or mpz_ext (1)
    mpz_cell*  m_ptr;     // pointer to cell (large numbers only)
    // ...
};
```

`mpz` is the user-facing value type. Every `mpz` object is one of two kinds:

| `m_kind` | `m_val` | `m_ptr` | Meaning |
|---|---|---|---|
| `mpz_small` | the actual integer value | `nullptr` | number fits in a machine `int` — no heap allocation |
| `mpz_large` | sign: +1 or -1 | pointer to `mpz_cell` | number is too large for `m_val` — stored on the heap |

The `m_owner` flag indicates whether the `mpz` object owns its `mpz_cell` (`mpz_self` — responsible for freeing it) or references an external buffer (`mpz_ext` — e.g. a stack-allocated cell in `mpz_stack`).

**Why it exists.** The dual representation avoids heap allocation for the common case where numbers are small (fit in `int`). This makes basic arithmetic on small values as fast as native `int` operations. Only when a result overflows does the library transparently promote to the large representation. `mpz` is move-only — it cannot be copied because copying a large `mpz` would require allocating a new cell, and the library is designed for explicit memory management through `mpz_manager`.

### `mpz_manager` — the public interface

```cpp
class mpz_manager {
    mpn_manager<digit_t> m_mpn_manager;  // low-level unsigned arithmetic engine
    // ...
public:
    // all public methods listed below
};
```

`mpz_manager` is the central facade. It owns the low-level unsigned arithmetic engine (`mpn_manager`), manages all heap allocation for `mpz_cell` objects, and dispatches every operation to either a fast inline path (when both operands are small) or a multi-precision path (when at least one operand is large). Users are expected to create one `mpz_manager` instance and use it for all operations — it tracks internal state such as the pre-computed value of 2⁶⁴.

> **Memory ownership rule.** Every `mpz` created or modified by `mpz_manager` must eventually be released with `del()` or `reset()`. The helper template `_scoped_numeral<mpz_manager>` provides RAII-style automatic cleanup.

---

## `mpz_manager` — Public API Reference

### Construction & lifecycle

| Method | Signature | Description |
|---|---|---|
| **Constructor** | `mpz_manager()` | Initialises the manager and pre-computes internal constants (2⁶⁴ and INT_MIN's large representation). |
| **Destructor** | `~mpz_manager()` | Releases all internally owned constants. Does **not** free user-created `mpz` objects — those must be released explicitly via `del()`. |
| `mk_z` | `static mpz mk_z(int val)` | Factory function: creates a small `mpz` from an `int`. The most concise way to create a temporary value for inline use. |
| `del` | `void del(mpz & a)` | Frees any heap memory owned by `a` and resets it to small zero. Safe to call on small `mpz` objects (no-op). |
| `del` (static) | `static void del(mpz_manager* m, mpz & a)` | Static variant — useful when only a pointer to the manager is available. |
| `reset` | `void reset(mpz & a)` | Deallocates `a`'s cell and sets `a` to zero. Equivalent to `del(a)` followed by `set(a, 0)`. |
| `dup` | `mpz dup(const mpz & source)` | Creates and returns a copy of `source`. The returned `mpz` owns its own cell (if any) and must be released. |

### Value assignment (`set` overloads)

| Method | Signature | Description |
|---|---|---|
| `set` | `void set(mpz & target, mpz const & source)` | Deep-copies `source` into `target`. Handles both small and large representations. |
| `set` | `void set(mpz & a, int val)` | Assigns a plain `int` value — `a` becomes a small `mpz`. |
| `set` | `void set(mpz & a, unsigned val)` | Assigns an `unsigned`. If `val > INT_MAX`, promotes to the large representation. |
| `set` | `void set(mpz & a, int64_t val)` | Assigns a 64-bit signed integer. Uses the fast path if `val` fits in `int`. |
| `set` | `void set(mpz & a, uint64_t val)` | Assigns a 64-bit unsigned integer. Promotes to large representation if needed. |
| `set` | `void set(mpz & a, char const * val)` | Parses a decimal C-string (e.g. `"12345"`, `"-42"`) and stores the result in `a`. Supports arbitrarily long strings. |
| `set_digits` | `void set_digits(mpz & target, unsigned sz, digit_t const * digits)` | Constructs a number from a raw little-endian digit array. Leading zero digits are stripped. The result is always positive. |
| `swap` | `void swap(mpz & a, mpz & b) noexcept` | Swaps two `mpz` values in O(1) by exchanging pointers and flags. |

### Arithmetic operations

| Method | Signature | Description |
|---|---|---|
| `add` | `void add(mpz const & a, mpz const & b, mpz & c)` | Computes `c = a + b`. |
| `sub` | `void sub(mpz const & a, mpz const & b, mpz & c)` | Computes `c = a - b`. |
| `mul` | `void mul(mpz const & a, mpz const & b, mpz & c)` | Computes `c = a * b`. |
| `addmul` | `void addmul(mpz const & a, mpz const & b, mpz const & c, mpz & d)` | Computes `d = a + b * c`. Optimised when `b` is ±1. |
| `submul` | `void submul(mpz const & a, mpz const & b, mpz const & c, mpz & d)` | Computes `d = a - b * c`. Optimised when `b` is ±1. |
| `inc` | `void inc(mpz & a)` | Increments `a` by 1 in place. |
| `dec` | `void dec(mpz & a)` | Decrements `a` by 1 in place. |
| `neg` | `void neg(mpz & a)` | Negates `a` in place. For small numbers, flips `m_val`; for large, flips the sign in `m_val`. |
| `abs` | `void abs(mpz & a)` | Replaces `a` with its absolute value in place. |

### Division & modular arithmetic

| Method | Signature | Description |
|---|---|---|
| `machine_div_rem` | `void machine_div_rem(mpz const & a, mpz const & b, mpz & q, mpz & r)` | Truncated division: `q = trunc(a / b)`, `r = a - q * b`. The remainder has the same sign as the dividend. |
| `machine_div` | `void machine_div(mpz const & a, mpz const & b, mpz & c)` | Truncated quotient only: `c = trunc(a / b)`. Throws `std::runtime_error` if `b` is zero. |
| `rem` | `void rem(mpz const & a, mpz const & b, mpz & c)` | Remainder of truncated division: `c = a - trunc(a/b) * b`. |
| `div` | `void div(mpz const & a, mpz const & b, mpz & c)` | **Floor division**: `c = floor(a / b)`. Adjusts the quotient when the dividend is negative and there is a non-zero remainder. |
| `mod` | `void mod(mpz const & a, mpz const & b, mpz & c)` | **Floor modulus**: `c = a - floor(a/b) * b`. The result always has the same sign as the divisor (or is zero). |
| `div_gcd` | `void div_gcd(mpz const & a, mpz const & b, mpz & c)` | Division shortcut for use inside GCD computations. If `b == 1`, copies `a`; otherwise calls `machine_div`. |
| `mod2k` | `mpz mod2k(mpz const & a, unsigned k)` | Returns `a mod 2^k`. Implemented via bit masking — no full division is performed. Returns a new `mpz` (caller must release it). |

### GCD, LCM & number theory

| Method | Signature | Description |
|---|---|---|
| `gcd` | `void gcd(mpz const & a, mpz const & b, mpz & c)` | Computes `c = gcd(a, b)` using the binary GCD algorithm. The result is always non-negative. |
| `gcd` (multi) | `void gcd(unsigned sz, mpz const * as, mpz & g)` | Computes the GCD of an array of `sz` numbers. |
| `gcd` (extended) | `void gcd(mpz const & r1, mpz const & r2, mpz & a, mpz & b, mpz & g)` | Extended Euclidean algorithm: finds `a`, `b`, `g` such that `r1 * a + r2 * b = g = gcd(r1, r2)`. |
| `lcm` | `void lcm(mpz const & a, mpz const & b, mpz & c)` | Computes `c = lcm(a, b)`. Handles edge cases: `lcm(a, 1) = a`, `lcm(a, a) = a`, `lcm(a, 0) = 0`. |
| `divides` | `bool divides(mpz const & a, mpz const & b)` | Returns `true` if `a` divides `b` evenly (i.e. `b % a == 0`). Defines `0 | 0` as true. |

### Bitwise operations

All bitwise operations require non-negative operands (asserted in debug builds). For large numbers, they process the value in 64-bit chunks.

| Method | Signature | Description |
|---|---|---|
| `bitwise_or` | `void bitwise_or(mpz const & a, mpz const & b, mpz & c)` | `c = a | b`. |
| `bitwise_and` | `void bitwise_and(mpz const & a, mpz const & b, mpz & c)` | `c = a & b`. |
| `bitwise_xor` | `void bitwise_xor(mpz const & a, mpz const & b, mpz & c)` | `c = a ^ b`. |
| `bitwise_not` | `void bitwise_not(unsigned sz, mpz const & a, mpz & c)` | `c = ~a`, masked to `sz` bits. The width parameter `sz` is required because the bitwise NOT of an arbitrary-precision integer is conceptually infinite. |
| `get_bit` | `bool get_bit(mpz const & a, unsigned index)` | Returns the value of the bit at position `index` (0 = LSB). Returns `false` if `index` exceeds the number's bit width. |

### Comparison & predicates

| Method | Signature | Description |
|---|---|---|
| `eq` | `bool eq(mpz const & a, mpz const & b)` | Equality test. Fast path for two small values. |
| `neq` | `bool neq(mpz const & a, mpz const & b)` | Inequality test. |
| `lt` | `bool lt(mpz const & a, mpz const & b)` | Less-than. Also has an overload `lt(mpz const & a, int b)`. |
| `gt` | `bool gt(mpz const & a, mpz const & b)` | Greater-than. |
| `le` | `bool le(mpz const & a, mpz const & b)` | Less-than-or-equal. |
| `ge` | `bool ge(mpz const & a, mpz const & b)` | Greater-than-or-equal. |
| `is_pos` | `static bool is_pos(mpz const & a)` | True if `a > 0`. |
| `is_neg` | `static bool is_neg(mpz const & a)` | True if `a < 0`. |
| `is_zero` | `static bool is_zero(mpz const & a)` | True if `a == 0`. |
| `is_nonpos` | `static bool is_nonpos(mpz const & a)` | True if `a <= 0`. |
| `is_nonneg` | `static bool is_nonneg(mpz const & a)` | True if `a >= 0`. |
| `sign` | `static int sign(mpz const & a)` | Returns -1, 0, or +1. |
| `is_small` | `static bool is_small(mpz const & a)` | True if `a` uses the inline (non-heap) representation. |

### Power, logarithm & bit-size

| Method | Signature | Description |
|---|---|---|
| `power` | `void power(mpz const & a, unsigned p, mpz & b)` | Computes `b = a^p`. Special-cased for `a = 0, 1, 2`; uses exponentiation by squaring for general bases. |
| `log2` | `unsigned log2(mpz const & a)` | Returns `floor(log2(a))`. Requires `a > 0`. |
| `mlog2` | `unsigned mlog2(mpz const & a)` | Returns `floor(log2(abs(a)))` for negative `a`. |
| `bitsize` | `unsigned bitsize(mpz const & a)` | Returns the minimum number of bits needed to represent `abs(a)`. |
| `is_power_of_two` | `bool is_power_of_two(mpz const & a)` | True if `a` is a positive power of two. |
| `is_power_of_two` (with shift) | `bool is_power_of_two(mpz const & a, unsigned & shift)` | Same, but also outputs the exponent. |
| `next_power_of_two` | `unsigned next_power_of_two(mpz const & a)` | Returns the exponent `k` such that `2^k` is the smallest power of two greater than `a`. Returns 0 for `a <= 0` or `a == 1`. |

### Scaling by powers of two

| Method | Signature | Description |
|---|---|---|
| `mul2k` | `void mul2k(mpz & a, unsigned k)` | Multiplies `a` by `2^k` in place. Implemented via bit shifts. |
| `machine_div2k` | `void machine_div2k(mpz & a, unsigned k)` | Divides `a` by `2^k` in place (truncated toward zero). Implemented via bit shifts. |

### Type conversion & extraction

| Method | Signature | Description |
|---|---|---|
| `is_uint64` | `bool is_uint64(mpz const & a) const` | True if `a` is non-negative and fits in `uint64_t`. |
| `is_int64` | `bool is_int64(mpz const & a) const` | True if `a` fits in `int64_t`. |
| `get_uint64` | `uint64_t get_uint64(mpz const & a) const` | Extracts the value as `uint64_t`. Caller must ensure `is_uint64(a)`. |
| `get_int64` | `int64_t get_int64(mpz const & a) const` | Extracts the value as `int64_t`. Caller must ensure `is_int64(a)`. |
| `is_uint` | `bool is_uint(mpz const & a) const` | True if `a` fits in `unsigned`. |
| `get_uint` | `unsigned get_uint(mpz const & a) const` | Extracts the value as `unsigned`. |
| `is_int` | `bool is_int(mpz const & a) const` | True if `a` fits in `int` (exclusive bounds). |
| `get_int` | `int get_int(mpz const & a) const` | Extracts the value as `int`. |
| `get_double` | `double get_double(mpz const & a) const` | Converts to `double` (may lose precision for very large numbers). |
| `get_least_significant` | `digit_t get_least_significant(mpz const & a)` | Returns the least significant digit of `abs(a)`. |
| `decompose` | `bool decompose(mpz const & a, std::vector<digit_t> & digits)` | Writes the little-endian digit representation into `digits`. Returns `true` if `a` is negative. |

### String & stream output

| Method | Signature | Description |
|---|---|---|
| `to_string` | `std::string to_string(mpz const & a) const` | Returns the decimal string representation (e.g. `"12345"`, `"-42"`). |
| `display` | `void display(std::ostream & out, mpz const & a) const` | Writes the decimal representation to an output stream. |
| `display_smt2` | `void display_smt2(std::ostream & out, mpz const & a, bool decimal) const` | Writes the value in [SMT-LIB 2.0](http://smtlib.cs.uiowa.edu/) format. If `decimal` is true, appends `.0` for real-number contexts. Negative values are wrapped as `(- ...)`. |
| `display_hex` | `void display_hex(std::ostream & out, mpz const & a, unsigned num_bits) const` | Writes the value as a zero-padded hexadecimal string of `num_bits / 4` characters. |
| `display_bin` | `void display_bin(std::ostream & out, mpz const & a, unsigned num_bits) const` | Writes the value as a zero-padded binary string of `num_bits` characters. |

---

## Building

### Requirements

- **C++20** compiler (GCC 10+, Clang 12+, MSVC 2019 16.11+)
- **CMake 3.21+**

### Standalone build

```bash
git clone https://github.com/<your-org>/ez3i.git
cd ez3i
cmake -B build -DEZ3I_BUILD_TESTS=ON
cmake --build build
ctest --test-directory build
```

### CMake options

| Option | Default | Description |
|---|---|---|
| `BUILD_SHARED_LIBS` | `OFF` | Build as shared (`ON`) or static (`OFF`) library. |
| `EZ3I_BUILD_TESTS` | `ON` (top-level) | Build unit tests. |
| `EZ3I_INSTALL` | `ON` (top-level) | Generate install targets. |

### Installation

```bash
cmake --install build --prefix /usr/local
```

---

## Integration

### FetchContent (recommended)

```cmake
include(FetchContent)
FetchContent_Declare(
    ez3i
    GIT_REPOSITORY https://github.com/abramov7613/ez3i.git
    GIT_TAG        v1.0.0
)
FetchContent_MakeAvailable(ez3i)

target_link_libraries(my_app PRIVATE ez3i)
```

### add_subdirectory

```cmake
add_subdirectory(third_party/ez3i)
target_link_libraries(my_app PRIVATE ez3i)
```

---

## Usage Examples

### Example 1: Basic arithmetic with large numbers

```cpp
#include <ez3i/mpz.h>
#include <iostream>

int main() {
    mpz_manager mgr;

    // Create two large numbers from strings
    mpz a, b;
    mgr.set(a, "123456789012345678901234567890");
    mgr.set(b, "987654321098765432109876543210");

    // Addition
    mpz sum;
    mgr.add(a, b, sum);
    std::cout << "a + b = " << mgr.to_string(sum) << "\n";

    // Multiplication
    mpz product;
    mgr.mul(a, b, product);
    std::cout << "a * b = " << mgr.to_string(product) << "\n";

    // Division with remainder
    mpz q, r;
    mgr.machine_div_rem(product, a, q, r);
    std::cout << "product / a = " << mgr.to_string(q) << "\n";
    std::cout << "product % a = " << mgr.to_string(r) << "\n";

    // Cleanup
    mgr.del(a);
    mgr.del(b);
    mgr.del(sum);
    mgr.del(product);
    mgr.del(q);
    mgr.del(r);

    return 0;
}
```

### Example 2: Using `mk_z` and GCD / LCM

This example demonstrates the `mk_z` factory for creating temporary small values inline, and computes the GCD and LCM of several numbers.

```cpp
#include <ez3i/mpz.h>
#include <iostream>

int main() {
    mpz_manager mgr;

    // mk_z creates a small mpz from an int — perfect for temporaries
    mpz a = mpz_manager::mk_z(360);
    mpz b = mpz_manager::mk_z(84);

    // GCD
    mpz g;
    mgr.gcd(a, b, g);
    std::cout << "gcd(" << mgr.to_string(a)
              << ", " << mgr.to_string(b)
              << ") = " << mgr.to_string(g) << "\n";

    // LCM
    mpz l;
    mgr.lcm(a, b, l);
    std::cout << "lcm(" << mgr.to_string(a)
              << ", " << mgr.to_string(b)
              << ") = " << mgr.to_string(l) << "\n";

    // Extended Euclidean: find x, y such that 360*x + 84*y = gcd(360, 84)
    mpz x, y, g2;
    mgr.gcd(a, b, x, y, g2);
    std::cout << "360 * (" << mgr.to_string(x)
              << ") + 84 * (" << mgr.to_string(y)
              << ") = " << mgr.to_string(g2) << "\n";

    // Use mk_z for quick comparisons
    if (mgr.lt(g, mpz_manager::mk_z(20))) {
        std::cout << "GCD is less than 20\n";
    }

    // Cleanup
    mgr.del(a);
    mgr.del(b);
    mgr.del(g);
    mgr.del(l);
    mgr.del(x);
    mgr.del(y);
    mgr.del(g2);

    return 0;
}
```

### Example 3: Bitwise operations and power-of-two detection

```cpp
#include <ez3i/mpz.h>
#include <iostream>

int main() {
    mpz_manager mgr;

    mpz a, b;
    mgr.set(a, "1024");  // 2^10
    mgr.set(b, "768");

    // Bitwise AND / OR / XOR
    mpz and_r, or_r, xor_r;
    mgr.bitwise_and(a, b, and_r);
    mgr.bitwise_or(a, b, or_r);
    mgr.bitwise_xor(a, b, xor_r);

    std::cout << "1024 & 768 = " << mgr.to_string(and_r) << "\n";
    std::cout << "1024 | 768 = " << mgr.to_string(or_r) << "\n";
    std::cout << "1024 ^ 768 = " << mgr.to_string(xor_r) << "\n";

    // Power-of-two check
    unsigned shift;
    if (mgr.is_power_of_two(a, shift)) {
        std::cout << "1024 is 2^" << shift << "\n";
    }

    // Exponentiation
    mpz pow_r;
    mgr.power(a, 3, pow_r);  // 1024^3 = 2^30
    std::cout << "1024^3 = " << mgr.to_string(pow_r) << "\n";

    // log2
    std::cout << "log2(1024^3) = " << mgr.log2(pow_r) << "\n";

    // Cleanup
    mgr.del(a); mgr.del(b);
    mgr.del(and_r); mgr.del(or_r); mgr.del(xor_r);
    mgr.del(pow_r);

    return 0;
}
```

---

## Project Structure

```
ez3i/
├── CMakeLists.txt
├── include/
│   └── ez3i/
│       ├── mpz.h          # Public API: mpz, mpz_cell, mpz_manager
│       └── mpn.hpp        # Internal: unsigned multi-precision arithmetic
├── src/
│   └── mpz.cpp           # Implementation of mpz_manager
└── tests/
    └── CMakeLists.txt     # Unit tests
```

---

## License

This project is derived from Z3, which is licensed under the MIT License.

## Acknowledgements

- The arithmetic algorithms are based on Donald E. Knuth, *The Art of Computer Programming*, Vol. 2, §4.3.
- The original integer implementation comes from the [Z3 theorem prover](https://github.com/Z3Prover/z3) by Microsoft Research.
