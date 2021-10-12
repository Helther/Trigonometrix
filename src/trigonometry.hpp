/*
// FIle contains
*/
#pragma once


#include <cmath>
#include <type_traits>
#include <cassert>
#include "float_table.hpp"
#include "double_table.hpp"
#include "polynomials_coeffs.hpp"
#include <iostream>
#include <limits>


inline constexpr auto DEG_TO_RAD = 1.7453292519943295769236907684886E-2;
inline constexpr auto HALF_PI = 1.5707963267948966192;
inline constexpr auto QUARTER_PI = 7.853981633974483096E-1;
inline constexpr auto INV_QUARTER_PI = 1 / QUARTER_PI;
inline constexpr auto INV_HALF_PI = 1 / HALF_PI;
inline constexpr auto TWO_PI = 6.2831853071795864769;
// value that represents the moment when argument becomes too large standart reduction
inline constexpr auto RANGE_REDUCTION_SWITCH = std::numeric_limits<float>::max();// TODO devise fitting change point

namespace Trigonometrix
{

enum Octets
{
    Zero_Pi4     = 0,
    Pi4_Pi2      = 1,
    Pi2_Pi3by4   = 2,
    Pi3by4_Pi    = 3,
    PiZero_Pi4   = 4,
    PiPi4_Pi2    = 5,
    PiPi2_Pi3by4 = 6,
    PiPi3by4_2Pi = 7
};

enum Quads
{
    Zero_Pi2   = 0,
    Pi2_Pi     = 1,
    Pi_Pi3by2  = 2,
    Pi3by2_2Pi = 3
};

//=================================== INTERNAL ===============================//
    namespace _Internal
    {

    struct alignas (alignof (int)) ReductionRes
    {
        int quad;
        bool noReduciton;
    };

    inline constexpr uint32_t two_over_pi[] = { 0x0, 0x28be60db, 0x24e44152, 0x27f09d5f, 0x11f534dd, 0x3036d8a5, 0x1993c439, 0x107f945, 0x23abdebb, 0x31586dc9,
    0x6e3a424, 0x374b8019, 0x92eea09, 0x3464873f, 0x21deb1cb, 0x4a69cfb, 0x288235f5, 0xbaed121, 0xe99c702, 0x1ad17df9,
    0x13991d6, 0xe60d4ce, 0x1f49c845, 0x3e2ef7e4, 0x283b1ff8, 0x25fff781, 0x1980fef2, 0x3c462d68, 0xa6d1f6d, 0xd9fb3c9,
    0x3cb09b74, 0x3d18fd9a, 0x1e5fea2d, 0x1d49eeb1, 0x3ebe5f17, 0x2cf41ce7, 0x378a5292, 0x3a9afed7, 0x3b11f8d5, 0x3421580c,
    0x3046fc7b, 0x1aeafc33, 0x3bc209af, 0x10d876a7, 0x2391615e, 0x3986c219, 0x199855f1, 0x1281a102, 0xdffd880, 0x135cc9cc,
    0x10606155
    };

    inline constexpr uint32_t pi_over_two[] = { 0x1, 0x2487ed51, 0x42d1846, 0x26263314, 0x1701b839, 0x28948127 };

    union DoubleUInt
    {
        uint64_t u;
        double   d;
    };

    // radix or base of representation
    #define RADIX (30)
    #define DIGITS 6

    DoubleUInt two_pow_pradix = { (uint64_t) (1023 + RADIX) << 52 };
    DoubleUInt two_pow_mradix = { (uint64_t) (1023 - RADIX) << 52 };
    DoubleUInt two_pow_two_mradix = { (uint64_t) (1023-2*RADIX) << 52 };

    #define tp_pradix two_pow_pradix.d
    #define tp_mradix two_pow_mradix.d

    // extended fixed point representation of double precision
    // floating point number.

    // extended fixed point representation of double precision
    // floating point number.
    // x = sign * [ sum_{i = 0 to 2} ( X[i] * 2^(index - i)*RADIX ) ]
    struct Eprep
    {
        uint32_t X[3];		// three 32 bit integers are sufficient to represnt double in base_30
        int index;			// exponent bias
        int sign;			// sign of double
    };

    double eprep_to_double(Eprep epx)
    {
        double res = 0.0;

        res += std::ldexp((double) epx.X[0], (epx.index - 0)*RADIX);
        res += std::ldexp((double) epx.X[1], (epx.index - 1)*RADIX);
        res += std::ldexp((double) epx.X[2], (epx.index - 2)*RADIX);

        return copysign(res, epx.sign);
    }

    Eprep double_to_eprep(double x)
    {
        Eprep result;

        result.sign = (std::signbit( (float) x ) == 0) ? 1 : -1;
        x = fabs( x );

        int index = 0;
        while( x > tp_pradix ) {
            index++;
            x *= tp_mradix;
        }
        while( x < 1 ) {
            index--;
            x *= tp_pradix;
        }

        result.index = index;
        int i = 0;
        result.X[0] = result.X[1] = result.X[2] = 0;
        while( x != 0.0 ) {
            result.X[i] = (uint32_t) x;
            x = (x - (double) result.X[i]) * tp_pradix;
            i++;
        }
        return result;
    }
    // TODO make constexpr
    // Payne-Hayek Argument Reduction for Huge Arguments: Good to the Last Bit
    // KronosGroup SYCL implementation
    template<typename T> constexpr ReductionRes payneHayekRangeReduce(T& arg) noexcept
    {
                if ((arg >= 0 ? arg : -arg) <= HALF_PI)
            return {0, true};

        double x = arg;

        // After computation result[0] contains integer part while result[1]....result[DIGITS-1]
        // contain fractional part. So we are doing computation with (DIGITS-1)*RADIX precision.
        // Default DIGITS=6 and RADIX=30 so default precision is 150 bits. Kahan-McDonald algorithm
        // shows that a double precision x, closest to pi/2 is 6381956970095103 x 2^797 which can
        // cause 61 digits of cancellation in computation of f = x*2/pi - floor(x*2/pi) ... thus we need
        // at least 114 bits (61 leading zeros + 53 bits of mentissa of f) of precision to accurately compute
        // f in double precision. Since we are using 150 bits (still an overkill), we should be safe. Extra
        // bits can act as guard bits for correct rounding.
        uint64_t result[DIGITS+2];

        // compute extended precision representation of x
        Eprep epx = double_to_eprep(x);
        int index = epx.index;
        int i, j;
        // extended precision multiplication of 2/pi*x .... we will loose at max two RADIX=30 bit digits in
        // the worst case
        for(i = 0; i < (DIGITS+2); i++) {
            result[i] = 0;
            result[i] += ((index + i - 0) >= 0) ? ((uint64_t) two_over_pi[index + i - 0] * (uint64_t) epx.X[0]) : 0;
            result[i] += ((index + i - 1) >= 0) ? ((uint64_t) two_over_pi[index + i - 1] * (uint64_t) epx.X[1]) : 0;
            result[i] += ((index + i - 2) >= 0) ? ((uint64_t) two_over_pi[index + i - 2] * (uint64_t) epx.X[2]) : 0;
        }

        // Carry propagation.
        uint64_t tmp;
        for(i = DIGITS+2-1; i > 0; i--) {
            tmp = result[i] >> RADIX;
            result[i - 1] += tmp;
            result[i] -= (tmp << RADIX);
        }

        // we dont ned to normalize the integer part since only last two bits of this will be used
        // subsequently algorithm which remain unaltered by this normalization.
        // tmp = result[0] >> RADIX;
        // result[0] -= (tmp << RADIX);
        unsigned int N = (unsigned int) result[0];

        // if the result is > pi/4, bring it to (-pi/4, pi/4] range. Note that testing if the final
        // x_star = pi/2*(x*2/pi - k) > pi/4 is equivalent to testing, at this stage, if r[1] (the first fractional
        // digit) is greater than (2^RADIX)/2 and substracting pi/4 from x_star to bring it to mentioned
        // range is equivalent to substracting fractional part at this stage from one and changing the sign.
        int sign = 1;
        if(result[1] > (uint64_t)(1 << (RADIX - 1))) {
            for(i = 1; i < (DIGITS + 2); i++)
                result[i] = (~((unsigned int)result[i]) & 0x3fffffff);
            N += 1;
            sign = -1;
        }

        // Again as per Kahan-McDonald algorithim there may be 61 leading zeros in the worst case
        // (when x is multiple of 2/pi very close to an integer) so we need to get rid of these zeros
        // and adjust the index of final result. So in the worst case, precision of comupted result is
        // 90 bits (150 bits original bits - 60 lost in cancellation).
        int ind = 1;
        for(i = 1; i < (DIGITS+2); i++) {
            if(result[i] != 0)
                break;
            else
                ind++;
        }

        uint64_t r[DIGITS-1];
        for(i = 0; i < (DIGITS-1); i++) {
            r[i] = 0;
            for(j = 0; j <= i; j++) {
                r[i] += (result[ind+i-j] * (uint64_t) pi_over_two[j]);
            }
        }
        for(i = (DIGITS-2); i > 0; i--) {
            tmp = r[i] >> RADIX;
            r[i - 1] += tmp;
            r[i] -= (tmp << RADIX);
        }
        tmp = r[0] >> RADIX;
        r[0] -= (tmp << RADIX);

        Eprep epr;
        epr.sign = epx.sign*sign;
        if(tmp != 0) {
            epr.index = -ind + 1;
            epr.X[0] = (uint32_t) tmp;
            epr.X[1] = (uint32_t) r[0];
            epr.X[2] = (uint32_t) r[1];
        }
        else {
            epr.index = -ind;
            epr.X[0] = (uint32_t) r[0];
            epr.X[1] = (uint32_t) r[1];
            epr.X[2] = (uint32_t) r[2];
        }

        arg = eprep_to_double( epr );
        return {int(epx.sign*N), false};
    }

    // simple fast additive reduction with possible loss of accuracy for large arguments
    template<typename T> constexpr ReductionRes addRangeReduce(T& arg, double maxRange, double invMaxRange) noexcept
    {
        if ((arg >= 0 ? arg : -arg) <= maxRange)
            return {0, true};

        const int quad = (int)(arg * invMaxRange);
        arg = arg - (T)(quad * maxRange);

        return {quad, false};
    }

    template<typename T> constexpr auto getNearestInt(T x)
    {
        static_assert(std::is_floating_point<T>(), "Invalid arg type");
        constexpr auto h = T(0.5) - std::numeric_limits<T>::epsilon();
        const T sign = x > 0 ? 1 : -1;
        return (int64_t)(x + sign*h);
		}
    //======== LU table implementation with range folding to 0...Pi/4 range ==//
    constexpr float sin_inner_table(float x) noexcept
    {
        if (x != QUARTER_PI)
            x = x / QUARTER_PI * TABLE_COUNT_F;
        else
            x = TABLE_COUNT_F - 1;
        const auto index = getNearestInt(x);
        const float diff = x - index;
        // if int is higher that arg get gradient from previous range
        const int gradIndex = (diff < 0 && index > 0) ? index-1 : index;
        //return fmaf(diff, SIN_GRAD_F[gradIndex], SIN_TABLE_F[index]);// TODO test fma version
        return SIN_TABLE_F[index] + diff * SIN_GRAD_F[gradIndex];
    }

    constexpr double sin_inner_table(double x) noexcept
    {
        if (x != QUARTER_PI)
            x = x / QUARTER_PI * TABLE_COUNT_D;
        else
            x = TABLE_COUNT_D - 1;
        const auto index = getNearestInt(x);
        const float diff = x - index;

        const int gradIndex = (diff < 0 && index > 0) ? index-1 : index;

        return SIN_TABLE_D[index] + diff * SIN_GRAD_D[gradIndex];
    }

    constexpr float cos_inner_table(float x) noexcept
    {
        if (x != QUARTER_PI)
            x = x / QUARTER_PI * TABLE_COUNT_F;
        else
            x = TABLE_COUNT_F - 1;
        const auto index = getNearestInt(x);
        const float diff = x - index;

        const int gradIndex = (diff < 0 && index > 0) ? index-1 : index;

        return COS_TABLE_F[index] + diff * COS_GRAD_F[gradIndex];
    }

    constexpr double cos_inner_table(double x) noexcept
    {
        if (x != QUARTER_PI)
            x = x / QUARTER_PI * TABLE_COUNT_D;
        else
            x = TABLE_COUNT_D - 1;
        const auto index = getNearestInt(x);
        const float diff = x - index;

        const int gradIndex = (diff < 0 && index > 0) ? index-1 : index;

        return COS_TABLE_D[index] + diff * COS_GRAD_D[gradIndex];
    }
//============================= Polynomial implementation ====================//
    template <typename T, unsigned accuracyDegree>
    constexpr T sin_inner_polinomial(T x) noexcept
    {
        if (x == 0)
            return 0;
        const T x2 = x * x;
        constexpr auto polySize = std::get<PolyIndex>(SIN_POLIES[accuracyDegree]);
        T res = std::get<PolyData>(SIN_POLIES[accuracyDegree])[polySize-1];
        for (int i = polySize-1; i >= 0; --i)
            res = res * x2 + std::get<PolyData>(SIN_POLIES[accuracyDegree])[i];

        return res *= x;
    }

    template <typename T, unsigned accuracyDegree>
    constexpr T cos_inner_polinomial(T x) noexcept
    {
        if (x == 0)
            return 1;
        const T x2 = x * x;
        constexpr auto polySize = std::get<PolyIndex>(COS_POLIES[accuracyDegree]);
        T res = std::get<PolyData>(COS_POLIES[accuracyDegree])[polySize-1];
        for (int i = polySize-1; i >= 0; --i)
            res = res * x2 + std::get<PolyData>(COS_POLIES[accuracyDegree])[i];

        return res;
    }

}
//============================= END INTERNAL =================================//

template <typename T> constexpr T radToDeg(T value) noexcept
{
    return value / DEG_TO_RAD;
}

template <typename T> constexpr T degToRad(T value) noexcept
{
    return value * DEG_TO_RAD;
}

//================================== Interface ===============================//
/* polyApprox - defines implementation
* TODO better interface is needed
* accuracyDegree - integer value in range 0..7 that scales up the accuracy
* of the approximation (index for number of terms for polinomial)
*
* accuracy map(relation of index to number of accurate digits in fractional part of the result):
* 0 - 3 (0.00049 average absolute error)
* 1 - 5
* 2 - 6
* 3 - 7
* 4 - 9
* 5 - 11
* 6 - 13
* 7 - 15
*/
template <typename T, bool polyApprox = true, std::size_t accuracyDegree = accuracy<T>>
constexpr T cos(T x) noexcept requires(std::is_floating_point<T>::value)
{
    //static_assert (accuracyDegree >, );
    const auto range = polyApprox ? HALF_PI : QUARTER_PI;
    const auto invRange = polyApprox ? INV_HALF_PI : INV_QUARTER_PI;
    const _Internal::ReductionRes res = x > RANGE_REDUCTION_SWITCH ? _Internal::payneHayekRangeReduce(x) : _Internal::addRangeReduce(x, range, invRange);
    if constexpr (polyApprox)
    {
        if (res.noReduciton)
            return _Internal::cos_inner_polinomial<T,accuracyDegree>(x);
        const int sign = res.quad >= 0 ? 1 : -1;
        x *= sign;
        // split function period into 4 equal parts shifted by Pi/2
        switch ((res.quad*sign) & Pi3by2_2Pi)
        {
        case Zero_Pi2:
            return _Internal::cos_inner_polinomial<T,accuracyDegree>(x);
        case Pi2_Pi:
            return -_Internal::sin_inner_polinomial<T,accuracyDegree>(x);
        case Pi_Pi3by2:
            return -_Internal::cos_inner_polinomial<T,accuracyDegree>(x);
        case Pi3by2_2Pi:
            return _Internal::sin_inner_polinomial<T,accuracyDegree>(x);
        }
    }
    else
    {
        const int sign = x >= 0 ? 1 : -1; // get sign and negate it, so table won't out of range
        x *= sign;
        if (res.noReduciton)
            return _Internal::cos_inner_table(x);
        const int quad = res.quad >= 0 ? res.quad : -res.quad;
        // split function period into 8 equal parts shifted by Pi/4
        switch (quad & PiPi3by4_2Pi)
        {
        case Zero_Pi4:
             return _Internal::cos_inner_table(x);
        case Pi4_Pi2:
            return _Internal::sin_inner_table(QUARTER_PI - x);
        case Pi2_Pi3by4:
            return -_Internal::sin_inner_table(x);
        case Pi3by4_Pi:
            return -_Internal::cos_inner_table(QUARTER_PI - x);
        case PiZero_Pi4:
            return -_Internal::cos_inner_table(x);
        case PiPi4_Pi2:
            return -_Internal::sin_inner_table(QUARTER_PI - x);
        case PiPi2_Pi3by4:
            return _Internal::sin_inner_table(x);
        case PiPi3by4_2Pi:
            return _Internal::cos_inner_table(QUARTER_PI - x);
        }
    }
    assert(false && "invalid range");
}

template <typename T, bool polyApprox = true, std::size_t accuracyDegree = accuracy<T>>
constexpr T sin(T x) noexcept requires(std::is_floating_point<T>::value)
{
    if constexpr (polyApprox)
    {
        const auto range = polyApprox ? HALF_PI : QUARTER_PI;
        const auto invRange = polyApprox ? INV_HALF_PI : INV_QUARTER_PI;
        const _Internal::ReductionRes res = x > RANGE_REDUCTION_SWITCH ? _Internal::payneHayekRangeReduce(x) : _Internal::addRangeReduce(x, range, invRange);
        if (res.noReduciton)
            return _Internal::sin_inner_polinomial<T,accuracyDegree>(x);
        const int sign = res.quad >= 0 ? 1 : -1;
        x *= sign;
        switch ((res.quad*sign) & Pi3by2_2Pi)
        {
        case Zero_Pi2:
            return sign*_Internal::sin_inner_polinomial<T,accuracyDegree>(x);
        case Pi2_Pi:
            return sign*_Internal::cos_inner_polinomial<T,accuracyDegree>(x);
        case Pi_Pi3by2:
            return -sign*_Internal::sin_inner_polinomial<T,accuracyDegree>(x);
        case Pi3by2_2Pi:
            return -sign*_Internal::cos_inner_polinomial<T,accuracyDegree>(x);
        }
    }
    else
    {
        return cos(HALF_PI - x);
    }
    assert(false && "invalid range");
}


template <typename T, bool polyApprox = true, std::size_t accuracyDegree = accuracy<T>>
constexpr T cosDeg(T degrees) noexcept requires(std::is_floating_point<T>::value)
{
    return cos(degToRad(degrees));
}

template <typename T, bool polyApprox = true, std::size_t accuracyDegree = accuracy<T>>
constexpr T sinDeg(T degrees) noexcept requires(std::is_floating_point<T>::value)
{
    return sin(degToRad(degrees));
}

}

