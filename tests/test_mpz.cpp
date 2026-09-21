#include <limits>
#include <cstdint>
#include <string>
#include <random>
#include <gtest/gtest.h>
#include <boost/multiprecision/cpp_int.hpp>
#include "mpz.h"

using boost::multiprecision::cpp_int;
using manager_t = mpz_manager<false>;

constexpr auto I_MIN = std::numeric_limits<int>::min();
constexpr auto I_MAX = std::numeric_limits<int>::max();
constexpr auto I64_MIN = std::numeric_limits<std::int64_t>::min();
constexpr auto I64_MAX = std::numeric_limits<std::int64_t>::max();

namespace {

//=============================  HELPERS =====================================

cpp_int to_cpp_int(std::int64_t v) { return cpp_int(v); }
cpp_int to_cpp_int(const std::string& v) { return cpp_int(v); }

std::int64_t rand_by_minmax(std::int64_t min, std::int64_t max) {
  static std::random_device rd;
  static std::mt19937_64 gen(rd());
  if (min > max) std::swap(min, max);
  return std::uniform_int_distribution<std::int64_t>(min, max)(gen);
}

std::string rand_by_digits(unsigned digits = 1) {
  static std::random_device rd{};
  static std::mt19937 gen(rd());
  std::string result;
  if (!digits) return result;
  if (!std::uniform_int_distribution<int>(0,1)(gen)) {
    result += '-';
    digits++;
  }
  result += static_cast<char>(std::uniform_int_distribution<int>(49, 57)(gen));
  digits--;
  while (result.size() < digits) {
    result += static_cast<char>(std::uniform_int_distribution<int>(48, 57)(gen));
  }
  return result;
}

cpp_int random_cpp_int(std::mt19937_64& rng, unsigned bits)
{
    cpp_int value = 0;

    unsigned words = (bits + 63) / 64;
    for (unsigned i = 0; i < words; ++i)
    {
        value <<= 64;
        value += rng();
    }

    if (rng() & 1)
        value = -value;

    return value;
}

std::string cpp_to_string(const cpp_int& v)
{
    return v.convert_to<std::string>();
}

void tst_div2k(manager_t & m, mpz const & v, unsigned k) {
    mpz x, y, two(2), pw;
    m.machine_div2k(v, k, x);
    m.power(two, k, pw);
    m.machine_div(v, pw, y);
    bool is_eq = m.eq(x, y);
    (void)is_eq;
    EXPECT_TRUE(is_eq);
}

void tst_div2k(manager_t & m, int v, unsigned k) {
    mpz x;
    m.set(x, v);
    tst_div2k(m, x, k);
}

void tst_div2k(manager_t & m, char const * v, unsigned k) {
    mpz x;
    m.set(x, v);
    tst_div2k(m, x, k);
}

void tst_mul2k(manager_t & m, mpz const & v, unsigned k) {
    mpz x, y, two(2), pw;
    m.mul2k(v, k, x);
    m.power(two, k, pw);
    m.mul(v, pw, y);
    bool is_eq = m.eq(x, y);
    (void)is_eq;
    EXPECT_TRUE(is_eq);
}

void tst_mul2k(manager_t & m, int v, unsigned k) {
    mpz x;
    m.set(x, v);
    tst_mul2k(m, x, k);
}

void tst_mul2k(manager_t & m, char const * v, unsigned k) {
    mpz x;
    m.set(x, v);
    tst_mul2k(m, x, k);
}

//=============================  TESTS =====================================

TEST(MpzTest, OriginalTst1)
{
    manager_t m;
    char const * str = "1002034040050606089383838288182";
    mpz v;
    m.set(v, str);
    mpz v2, v3;
    m.mul(v, m.mk_z(-2), v2);
    m.add(v, v2, v3);
    m.neg(v3);
    EXPECT_TRUE(m.eq(v, v3));
    EXPECT_TRUE(m.le(v, v3));
    EXPECT_TRUE(m.ge(v, v3));
    EXPECT_TRUE(m.lt(v2, v));
    EXPECT_TRUE(m.le(v2, v));
    EXPECT_TRUE(m.gt(v, v2));
    EXPECT_TRUE(m.ge(v, v2));
    EXPECT_TRUE(m.neq(v, v2));
    EXPECT_TRUE(!m.neq(v, v3));
}

TEST(MpzTest, OriginalBug1)
{
    manager_t m;
    mpz v1;
    m.set(v1, "1002043949858757875676767675747473");
    mpz v2;
    m.sub(v1, v1, v2);
    EXPECT_TRUE(m.is_zero(v2));
}

TEST(MpzTest, OriginalBug3)
{
    manager_t m;
    mpz v1, v2;
    m.set(v1, I_MIN);
    m.set(v2, I_MAX);
    m.add(v2, m.mk_z(1), v2);
    m.neg(v1);
    EXPECT_TRUE(m.eq(v1, v2));
}

TEST(MpzTest, OriginalBug4)
{
    manager_t m;
    mpz x, y;
    m.set(y, static_cast<uint64_t>(4294967295ull));
    m.set(x, static_cast<uint64_t>(4026531839ull));
    mpz result1;
    m.bitwise_or(x, y, result1);

    mpz result2;
    m.set(result2, x);
    m.bitwise_or(result2, y, result2);

    EXPECT_TRUE(m.eq(result1, result2));
}

TEST(MpzTest, OriginalTst2k)
{
    manager_t m;
    tst_mul2k(m, 120, 32);
    tst_mul2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 22);
    tst_div2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 22);

    tst_div2k(m, 3, 1);
    tst_div2k(m, 3, 2);
    tst_div2k(m, 3, 0);
    tst_div2k(m, 120, 32);
    tst_div2k(m, 120, 0);
    tst_div2k(m, 81, 2);
    tst_div2k(m, -3, 1);
    tst_div2k(m, -3, 2);
    tst_div2k(m, -3, 0);
    tst_div2k(m, -102, 4);
    tst_div2k(m, 0, 3);
    tst_div2k(m, 0, 1000);
    tst_div2k(m, 7, 10000);
    tst_div2k(m, -7, 1000);
    tst_div2k(m, -7, 2);
    tst_div2k(m, "1029384848584832828327176162636436484", 4);
    tst_div2k(m, "1029384848584832828327176162636436484", 100);
    tst_div2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 1);
    tst_div2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 0);
    tst_div2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 4);
    tst_div2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 7);
    tst_div2k(m, "102938484858483282832717616263643648433838737661626264364583983298239291919", 100);
    tst_div2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 100);
    tst_div2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 177);
    tst_div2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 77);
    tst_div2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 32);
    tst_div2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 64);
    tst_div2k(m, "-11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 64);
    tst_div2k(m, "-11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 128);
    tst_div2k(m, "-1092983874757371817626399990000000", 100);
    tst_div2k(m, "-109298387475737181762639999000000231", 8);
    tst_div2k(m, "-109298387475737181762639999000000231", 16);
    tst_div2k(m, "-109298387475737181762639999000000231", 17);
    tst_div2k(m, "-109298387475737181762639999000000231", 32);

    tst_mul2k(m, 3, 1);
    tst_mul2k(m, 3, 2);
    tst_mul2k(m, 3, 0);
    tst_mul2k(m, 120, 32);
    tst_mul2k(m, 120, 0);
    tst_mul2k(m, 81, 2);
    tst_mul2k(m, -3, 1);
    tst_mul2k(m, -3, 2);
    tst_mul2k(m, -3, 0);
    tst_mul2k(m, -102, 4);
    tst_mul2k(m, 0, 3);
    tst_mul2k(m, 0, 1000);
    tst_mul2k(m, 7, 10000);
    tst_mul2k(m, -7, 1000000);
    tst_mul2k(m, -7, 2);
    tst_mul2k(m, "1029384848584832828327176162636436484", 4);
    tst_mul2k(m, "1029384848584832828327176162636436484", 100);
    tst_mul2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 1);
    tst_mul2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 0);
    tst_mul2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 4);
    tst_mul2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 22);
    tst_mul2k(m, "102938484858483282832717616263643648481827437292943727163646457588332211", 7);
    tst_mul2k(m, "102938484858483282832717616263643648433838737661626264364583983298239291919", 100);
    tst_mul2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 100);
    tst_mul2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 177);
    tst_mul2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 77);
    tst_mul2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 32);
    tst_mul2k(m, "11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 64);
    tst_mul2k(m, "-11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 64);
    tst_mul2k(m, "-11579208923731619542357098500868790785326998466564056403945758400793872761761638458372762", 128);
    tst_mul2k(m, "-1092983874757371817626399990000000", 100);
    tst_mul2k(m, "-109298387475737181762639999000000231", 8);
    tst_mul2k(m, "-109298387475737181762639999000000231", 16);
    tst_mul2k(m, "-109298387475737181762639999000000231", 17);
    tst_mul2k(m, "-109298387475737181762639999000000231", 32);
}

TEST(MpzTest, OriginalTstintminbug)
{
    manager_t m;
    mpz intmin(INT_MIN);
    mpz big;
    mpz expected;
    mpz r;
    m.set(big, static_cast<uint64_t>(UINT64_MAX));
    m.set(expected, "18446744075857035263");
    m.sub(big, intmin, r);
    std::cout << "r: " << m.to_string(r) << "\nexpected: " << m.to_string(expected) << "\n";
    EXPECT_TRUE(m.eq(r, expected));
}

TEST(MpzTest, OriginalTstint64minbug)
{
    manager_t m;
    mpz intmin;
    mpz test;
    m.set(test, "-9223372036854775808");
    m.set(intmin, std::numeric_limits<int64_t>::min());
    std::cout << "minint: " << m.to_string(intmin) << "\n";
    EXPECT_TRUE(m.eq(test, intmin));
}

TEST(MpzTest, BoundaryInt64)
{
    manager_t m;

    mpz a;
    m.set(a, I64_MAX);

    EXPECT_EQ(m.get_int64(a), I64_MAX);
}

TEST(MpzTest, NegativeArithmetic)
{
    manager_t m;

    mpz a(-100);
    mpz b(7);
    mpz r;
    m.machine_div_rem(a, b, r, r);
    EXPECT_TRUE(manager_t::is_neg(r));

    mpz v1, v2;
    m.set(v1, I_MIN);
    m.set(v2, I_MAX);
    m.add(v2, m.mk_z(1), v2);
    m.neg(v1);
    EXPECT_TRUE(m.eq(v1, v2));
}

TEST(MpzTest, AddAcrossSmallBoundary)
{
    manager_t m;

    mpz a(I_MAX);
    mpz b(1);
    cpp_int c(I_MAX);
    c++;
    mpz r;
    m.add(a, b, r);

    EXPECT_EQ(m.to_string(r), c.str());
}

TEST(MpzTest, Multiplication)
{
    manager_t m;

    mpz a(1000000);
    mpz b(1000000);
    mpz r;

    m.mul(a, b, r);

    EXPECT_EQ(m.get_int64(r), 1000000000000LL);
}

TEST(MpzTest, DivisionSigns)
{
    manager_t m;

    mpz a(-100);
    mpz b(-7);
    mpz q, r;

    m.machine_div_rem(a, b, q, r);

    EXPECT_GT(m.get_int64(q), 0);
}

TEST(MpzTest, BitOperationsNegativeNumbers)
{
    manager_t m;

    mpz a(-1);
    mpz b(0xFF);
    mpz r;

    m.bitwise_and(a, b, r);

    EXPECT_EQ(m.get_int64(r), 255);
}

TEST(MpzTest, SmallArithmeticAgainstInt64)
{
    manager_t m;

    for (int64_t a = -100; a <= 100; ++a)
    {
        for (int64_t b = -100; b <= 100; ++b)
        {
            mpz x(a), y(b), r;

            m.add(x, y, r);

            EXPECT_EQ(m.get_int64(r), a + b);
        }
    }
}

TEST(MpzTest, GcdProperties)
{
    manager_t m;

    mpz a(48);
    mpz b(18);
    mpz g;

    m.gcd(a, b, g);

    EXPECT_EQ(m.get_int64(g), 6);
}

TEST(MpzTest, AdditionSmallRandom)
{
    manager_t m;

    for (int i = 0; i < 100; ++i)
    {
        int a = rand_by_minmax(0, I_MAX/2) ;
        int b = rand_by_minmax(0, I_MAX/2) ;

        mpz x(a), y(b), r;
        m.add(x, y, r);

        EXPECT_EQ(m.get_int64(r), a+b) ;
    }
}

TEST(MpzTest, AdditionBigRandom)
{
    manager_t m;

    for (int i = 0; i < 100; ++i)
    {
        std::string a = rand_by_digits(36) ;
        std::string b = rand_by_digits(32) ;

        mpz x, y, r, r2;
        m.set(x, a.c_str());
        m.set(y, b.c_str());
        m.add(x, y, r);
        std::string r2_str = (cpp_int(a) + cpp_int(b)).str();
        m.set(r2, r2_str.c_str());

        EXPECT_EQ(m.to_string(r), r2_str) ;
        EXPECT_TRUE(m.eq(r, r2));
    }
}

TEST(MpzTest, SubtractionBigRandom)
{
    manager_t m;

    for (int i = 0; i < 100; ++i)
    {
        std::string a = rand_by_digits(42) ;
        std::string b = rand_by_digits(42) ;

        mpz x, y, r, r2;
        m.set(x, a.c_str());
        m.set(y, b.c_str());
        m.sub(x, y, r);
        std::string r2_str = (cpp_int(a) - cpp_int(b)).str();
        m.set(r2, r2_str.c_str());

        EXPECT_EQ(m.to_string(r), r2_str) ;
        EXPECT_TRUE(m.eq(r, r2));
    }
}

TEST(MpzTest, MultiplicationBigRandom)
{
    manager_t m;

    for (int i = 0; i < 100; ++i)
    {
        std::string a = rand_by_digits(32) ;
        std::string b = rand_by_digits(6) ;

        mpz x, y, r, r2;
        m.set(x, a.c_str());
        m.set(y, b.c_str());
        m.mul(x, y, r);
        std::string r2_str = (cpp_int(a) * cpp_int(b)).str();
        m.set(r2, r2_str.c_str());

        EXPECT_EQ(m.to_string(r), r2_str) ;
        EXPECT_TRUE(m.eq(r, r2));
    }
}

TEST(MpzTest, Division64)
{
    manager_t m;
    std::mt19937_64 rng(999);

    for (int i = 0; i < 500; ++i)
    {
        int64_t a = static_cast<int64_t>(rng() % 1000000);
        int64_t b = static_cast<int64_t>(rng() % 1000) + 1;

        mpz x(a), y(b), q, r;

        m.machine_div_rem(x, y, q, r);

        EXPECT_EQ(m.get_int64(q), a / b);
        EXPECT_EQ(m.get_int64(r), a % b);
    }
}

TEST(MpzTest, DivisionBigRandom)
{
    manager_t m;

    for (int i = 0; i < 100; ++i)
    {
        std::string a = rand_by_digits(68) ;
        std::string b = rand_by_digits(42) ;

        mpz x, y, q, r, q2;
        m.set(x, a.c_str());
        m.set(y, b.c_str());
        m.machine_div_rem(x, y, q, r);
        std::string q2_str = (cpp_int(a) / cpp_int(b)).str();
        m.set(q2, q2_str.c_str());

        EXPECT_EQ(m.to_string(q), q2_str) ;
        EXPECT_TRUE(m.eq(q, q2));
    }
}

TEST(MpzTest, GcdRandom)
{
    manager_t m;
    std::mt19937_64 rng(3);

    for (int i = 0; i < 200; ++i)
    {
        int64_t a = static_cast<int64_t>(rng() % 1000000);
        int64_t b = static_cast<int64_t>(rng() % 1000000);

        mpz x(a), y(b), g;

        m.gcd(x, y, g);

        EXPECT_GE(m.get_int64(g), 0);
    }
}

TEST(MpzTest, ShiftRegression)
{
    manager_t m;

    for (int i = 0; i < 64; ++i)
    {
        mpz a(1);
        m.mul2k(a, i);

        EXPECT_EQ(m.to_string(a),
                  (cpp_int(1) << i).convert_to<std::string>());
    }
}

TEST(MpzTest, StringRoundTrip)
{
    manager_t m;

    for (long long i = -100000; i <= 100000; i += 997)
    {
        mpz a(i);

        auto s = m.to_string(a);

        mpz b;
        m.set(b, s.c_str());

        EXPECT_TRUE(m.eq(a, b));
    }
}


}

