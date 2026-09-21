#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <boost/multiprecision/cpp_int.hpp>

#include "util/mpz.h"

using boost::multiprecision::cpp_int;
using manager_t = mpz_manager<false>;

namespace {

std::string text(const cpp_int& value) {
    return value.str();
}

cpp_int as_cpp(const mpz& value, manager_t& manager) {
    return cpp_int(manager.to_string(value));
}

cpp_int gcd_cpp(cpp_int a, cpp_int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        cpp_int r = a % b;
        a = b;
        b = r;
    }
    return a;
}

void set_value(manager_t& manager, mpz& result, const cpp_int& value) {
    const std::string value_text = text(value);
    manager.set(result, value_text.c_str());
}

cpp_int random_cpp_int(std::mt19937_64& generator, unsigned bits) {
    cpp_int result = 0;
    const unsigned words = (bits + 63) / 64;
    for (unsigned i = 0; i < words; ++i) {
        result <<= 64;
        result += generator();
    }
    if ((generator() & 1) != 0)
        result = -result;
    return result;
}

TEST(MpzTest1, BoundaryValuesRoundTrip) {
    manager_t manager;
    const std::vector<std::string> values = {
        "-2147483649", "-2147483648", "-2147483647",
        "-1", "0", "1",
        "2147483646", "2147483647", "2147483648", "2147483649",
        "-9223372036854775808", "9223372036854775807",
        "18446744073709551615"
    };

    for (const std::string& expected : values) {
        mpz value;
        manager.set(value, expected.c_str());
        EXPECT_EQ(manager.to_string(value), expected) << "value=" << expected;
    }
}

TEST(MpzTest1, BoundaryArithmetic) {
    manager_t manager;
    const cpp_int int_max = std::numeric_limits<int>::max();
    const cpp_int int_min = std::numeric_limits<int>::min();
    const cpp_int int64_max = std::numeric_limits<std::int64_t>::max();
    const cpp_int int64_min = std::numeric_limits<std::int64_t>::min();

    for (const cpp_int& value : {int_min, int_max, int64_min, int64_max}) {
        mpz input, one, result;
        set_value(manager, input, value);
        manager.set(one, 1);
        manager.add(input, one, result);
        EXPECT_EQ(as_cpp(result, manager), value + 1);
        manager.sub(input, one, result);
        EXPECT_EQ(as_cpp(result, manager), value - 1);
        manager.neg(input);
        EXPECT_EQ(as_cpp(input, manager), -value);
    }
}

TEST(MpzTest1, ArithmeticMatchesReference) {
    manager_t manager;
    std::mt19937_64 generator(0xA11CE1234ULL);

    for (int i = 0; i < 250; ++i) {
        const cpp_int a = random_cpp_int(generator, 512);
        const cpp_int b = random_cpp_int(generator, 512);
        mpz x, y, result;
        set_value(manager, x, a);
        set_value(manager, y, b);

        manager.add(x, y, result);
        EXPECT_EQ(as_cpp(result, manager), a + b);
        manager.sub(x, y, result);
        EXPECT_EQ(as_cpp(result, manager), a - b);
        manager.mul(x, y, result);
        EXPECT_EQ(as_cpp(result, manager), a * b);

        manager.addmul(x, y, x, result);
        EXPECT_EQ(as_cpp(result, manager), a + b * a);
        manager.submul(x, y, x, result);
        EXPECT_EQ(as_cpp(result, manager), a - b * a);
    }
}

TEST(MpzTest1, DivisionRemainderAndModuloMatchReference) {
    manager_t manager;
    std::mt19937_64 generator(0xD1V1234ULL);

    for (int i = 0; i < 250; ++i) {
        const cpp_int a = random_cpp_int(generator, 512);
        cpp_int b = random_cpp_int(generator, 256);
        if (b == 0)
            b = 1;

        mpz x, y, quotient, remainder, result;
        set_value(manager, x, a);
        set_value(manager, y, b);
        manager.machine_div_rem(x, y, quotient, remainder);

        EXPECT_EQ(as_cpp(quotient, manager), a / b);
        EXPECT_EQ(as_cpp(remainder, manager), a % b);
        EXPECT_EQ(as_cpp(quotient, manager) * b + as_cpp(remainder, manager), a);

        manager.mod(x, y, result);
        cpp_int expected_mod = a % b;
        if (expected_mod < 0)
            expected_mod += b < 0 ? -b : b;
        EXPECT_EQ(as_cpp(result, manager), expected_mod);
    }
}

TEST(MpzTest1, GcdLcmDividesAndExtendedGcd) {
    manager_t manager;
    const std::vector<std::pair<cpp_int, cpp_int>> cases = {
        {0, 0}, {0, -42}, {-48, 18}, {48, 18},
        {"12345678901234567890"_cpp, "9876543210"_cpp}
    };

    for (const auto& [a, b] : cases) {
        mpz x, y, gcd, lcm, coefficient_a, coefficient_b, bezout_gcd;
        set_value(manager, x, a);
        set_value(manager, y, b);
        manager.gcd(x, y, gcd);
        EXPECT_EQ(as_cpp(gcd, manager), gcd_cpp(a, b));

        manager.lcm(x, y, lcm);
        const cpp_int expected_lcm = (a == 0 || b == 0) ? 0 : (a * b) / gcd_cpp(a, b);
        EXPECT_EQ(as_cpp(lcm, manager), expected_lcm < 0 ? -expected_lcm : expected_lcm);

        EXPECT_EQ(manager.divides(x, y), a == 0 ? b == 0 : (b % a == 0));

        manager.gcd(x, y, coefficient_a, coefficient_b, bezout_gcd);
        EXPECT_EQ(as_cpp(bezout_gcd, manager), gcd_cpp(a, b));
        EXPECT_EQ(as_cpp(coefficient_a, manager) * a +
                  as_cpp(coefficient_b, manager) * b,
                  as_cpp(bezout_gcd, manager));
    }
}

TEST(MpzTest1, PowersShiftsAndModuloPowersOfTwo) {
    manager_t manager;
    std::mt19937_64 generator(99);
    const std::vector<unsigned> shifts = {0, 1, 31, 32, 33, 63, 64, 65, 127, 128, 129, 256, 1000};

    for (int i = 0; i < 100; ++i) {
        const cpp_int original = random_cpp_int(generator, 768);
        for (unsigned shift : shifts) {
            mpz value, shifted, divided;
            set_value(manager, value, original);
            manager.mul2k(value, shift);
            EXPECT_EQ(as_cpp(value, manager), original * (cpp_int(1) << shift));

            set_value(manager, value, original);
            manager.machine_div2k(value, shift);
            EXPECT_EQ(as_cpp(value, manager), original / (cpp_int(1) << shift));

            set_value(manager, value, original);
            shifted = manager.mod2k(value, shift);
            cpp_int expected = original % (cpp_int(1) << shift);
            if (expected < 0)
                expected += cpp_int(1) << shift;
            EXPECT_EQ(as_cpp(shifted, manager), expected);

            set_value(manager, value, 2);
            manager.power(value, shift, divided);
            EXPECT_EQ(as_cpp(divided, manager), cpp_int(1) << shift);
        }
    }
}

TEST(MpzTest1, BitOperationsAndBitQueries) {
    manager_t manager;
    std::mt19937_64 generator(123);

    for (int i = 0; i < 100; ++i) {
        cpp_int a = random_cpp_int(generator, 256);
        cpp_int b = random_cpp_int(generator, 256);
        if (a < 0) a = -a;
        if (b < 0) b = -b;

        mpz x, y, result;
        set_value(manager, x, a);
        set_value(manager, y, b);

        manager.bitwise_and(x, y, result);
        EXPECT_EQ(as_cpp(result, manager), a & b);
        manager.bitwise_or(x, y, result);
        EXPECT_EQ(as_cpp(result, manager), a | b);
        manager.bitwise_xor(x, y, result);
        EXPECT_EQ(as_cpp(result, manager), a ^ b);

        for (unsigned bit = 0; bit < 300; ++bit)
            EXPECT_EQ(manager.get_bit(x, bit), boost::multiprecision::bit_test(a, bit));
    }
}

TEST(MpzTest1, RootsAndNumberProperties) {
    manager_t manager;
    for (unsigned exponent = 0; exponent < 20; ++exponent) {
        mpz value, root;
        manager.power(mpz(3), exponent, value);
        EXPECT_TRUE(manager.is_perfect_square(value, root) == (exponent % 2 == 0));
        if (exponent % 2 == 0)
            EXPECT_EQ(manager.to_string(root), text(cpp_int(3) << 0));

        EXPECT_TRUE(manager.is_power_of_two(manager.power(mpz(2), exponent, value)));
    }

    mpz value;
    manager.set(value, 65);
    EXPECT_FALSE(manager.is_power_of_two(value));
    EXPECT_EQ(manager.log2(value), 6u);
    EXPECT_EQ(manager.bitsize(value), 7u);
}

TEST(MpzTest1, FormattingAndConversions) {
    manager_t manager;
    mpz value;
    manager.set(value, "255");

    std::ostringstream stream;
    manager.display_smt2(stream, value, true);
    EXPECT_EQ(stream.str(), "255.0");

    stream.str({});
    stream.clear();
    manager.display_hex(stream, value, 16);
    EXPECT_EQ(stream.str(), "00ff");

    stream.str({});
    stream.clear();
    manager.display_bin(stream, value, 8);
    EXPECT_EQ(stream.str(), "11111111");

    EXPECT_TRUE(manager.is_uint64(value));
    EXPECT_TRUE(manager.is_int64(value));
    EXPECT_EQ(manager.get_uint64(value), 255u);
    EXPECT_EQ(manager.get_int64(value), 255);
    EXPECT_EQ(manager.get_least_significant(value), 255u);
}

TEST(MpzTest1, DivisionByZeroThrows) {
    manager_t manager;
    mpz value(1), zero(0), result;
    EXPECT_ANY_THROW(manager.machine_div(value, zero, result));
}

} // namespace
