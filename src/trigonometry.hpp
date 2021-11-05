/*
 * File contains interface and internal workings of such functions as: sine,
 * cosine, tangent, arc-sin/cos/tan. Function interface allows for accuracy/speed
 * tradeoff.
*/
#pragma once


#include <cmath>
#include <type_traits>
#include <cassert>
#include "polynomials_coeffs.hpp"
#include "lut_generator.hpp"
#include <limits>
#include <inttypes.h>


inline constexpr auto DEG_TO_RAD = 1.7453292519943295769236907684886E-2;
inline constexpr auto HALF_PI = 1.5707963267948966192;
inline constexpr auto QUARTER_PI = 7.853981633974483096E-1;
inline constexpr auto INV_QUARTER_PI = 1 / QUARTER_PI;
inline constexpr auto INV_HALF_PI = 1 / HALF_PI;
inline constexpr auto TWO_PI = 6.2831853071795864769;
// value that represents the moment when argument becomes too large for standart reduction
inline constexpr auto RANGE_REDUCTION_SWITCH = std::numeric_limits<float>::max();// TODO devise fitting change point

namespace Trigonometrix
{

enum Octants
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

/*
 * Approximation of acos(x) with a rational function f such that the
 * worst absolute error is minimal. That is, pick the function that performs best
 * in the worst case, and with following restrictions: acos(0) = Pi/2, acos(1) = 0, acos(-1) =
 * Source: https://github.com/ruuda/convector/blob/master/tools/approx_acos.py
 * Accuracy: av abs error 3e-11, max abs error 0.0167, accuracy get worse when aproaching the limits.
 * Up to 2 times faster than std, depending on compiler.
 * Returns the arccosine of a in the range [0,pi], expecting a to be in the range [-1,+1].
*/
template<typename T> constexpr T acos(T x) noexcept requires(std::is_floating_point_v<T>)
{
    assert(x >= -1 && x <= 1 && "invalid argument range");
    constexpr T c1 = -0.939115566365855;
    constexpr T c2 =  0.9217841528914573;
    constexpr T c3 = -1.2845906244690837;
    constexpr T c4 =  0.295624144969963174;
    const T x2 = x * x;
    const T x3 = x * x * x;
    const T x4 = x * x * x * x;

    return HALF_PI + (c1*x + c2*x3) / (1 + c3*x2 + c4*x4);

}

// wrapper function to handle interger arguments
template<typename T> constexpr auto acos(T x) noexcept requires(std::is_integral_v<T>)
{
    return acos<double>(x);
}

// just a shifted version of acos implementation
// Returns the arc cosine of a in the range [-pi/2,pi/2], expecting a to be in the range [-1,+1].
template<typename T> constexpr T asin(T x) noexcept requires(std::is_floating_point_v<T>)
{
    return HALF_PI - acos<T>(x);
}

// wrapper function to handle interger arguments
template<typename T> constexpr auto asin(T x) noexcept requires(std::is_integral_v<T>)
{
    return HALF_PI - acos<double>(x);
}


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

        double temp = arg;// use double for less loss of significance
        const int quad = (int)(temp * invMaxRange);
        temp = temp - (double)(quad * maxRange);
        arg = temp;
        return {quad, false};
    }

    template<typename T> constexpr auto getNearestInt(T x)
    {
        static_assert(std::is_floating_point<T>(), "Invalid arg type");
        constexpr auto h = T(0.5) - std::numeric_limits<T>::epsilon();
        const T sign = x > 0 ? 1 : -1;
        return (int64_t)(x + sign*h);
    }

//======================= Generic LU table implementation ====================//

    // loose approximation of what size of the table whould be for a given relative error
    constexpr int constLUTSizeFromAcc(double relError, int ratio)
    {
        return int(M_PI / Trigonometrix::acos(1 - relError) / ratio) + 1;
    }
    // maps number of accurate significant needed to the max error value for constLUTSizeFromAcc
    inline constexpr std::array<double,SIN_COS_ACC_MAP_COUNT> SC_LUT_ACC_MAP =
    { 0.1,0.01,0.001,0.0001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001 };

    // Information for periodic function LUT generation
    template <typename T, T(*Func)(T), std::size_t foldingRatio, std::size_t acc>
    struct LUTInfo
    {
        static constexpr std::size_t size = constLUTSizeFromAcc(SC_LUT_ACC_MAP[acc],foldingRatio);
        static constexpr double startValue = 0;
        static constexpr double endValue = 2*M_PI / foldingRatio;
        static constexpr double step = endValue / size;
        static constexpr auto table = getLUT<T,Func,LUTInfo>;
    };
    /*
     * calculates value with gradient approximation betweeen value points, based
     * on a given TableInfo static struct descibed above
     */
    template <typename T, class TableInfo>
    constexpr T generic_inner_table(T x) noexcept
    {
        x = x / TableInfo::endValue * TableInfo::size;
        auto index = getNearestInt(x);
        if (index == TableInfo::size)
            --index;
        const T diff = x - index;

        const int gradIndex = (diff < 0 && index > 0) ? index-1 : index;
        // TODO test FMA version
        return std::get<0>(TableInfo::table[index]) + diff * std::get<1>(TableInfo::table[gradIndex]);
    }

//============================= Polynomial implementation ====================//
    template <typename T, std::size_t accuracy>
    constexpr T sin_inner_polinomial(T x) noexcept
    {
        if (x == 0)
            return x;
        const T x2 = x * x;
        constexpr auto polySize = std::get<PolyIndex>(SIN_POLIES[accuracy]);
        T res = std::get<PolyData>(SIN_POLIES[accuracy])[polySize-1];
        for (int i = polySize-2; i >= 0; --i)
            res = res * x2 + std::get<PolyData>(SIN_POLIES[accuracy])[i];

        return res * x;
    }

    template <typename T, std::size_t accuracy>
    constexpr T cos_inner_polinomial(T x) noexcept
    {
        if (x == 0)
            return 1;
        const T x2 = x * x;
        constexpr auto polySize = std::get<PolyIndex>(COS_POLIES[accuracy]);
        T res = std::get<PolyData>(COS_POLIES[accuracy])[polySize-1];
        for (int i = polySize-2; i >= 0; --i)
            res = res * x2 + std::get<PolyData>(COS_POLIES[accuracy])[i];

        return res;
    }

    template <typename T, bool fast>
    constexpr T tan_inner_polynomial(T x) noexcept
    {
        x *= INV_QUARTER_PI;
        const T x2 = x * x;
        if constexpr(fast)
        {
            return x * TAN_DEGREE_2[0] / (TAN_DEGREE_2[1] + x2);
        }
        else
        {
            return x * (TAN_DEGREE_4[0] + TAN_DEGREE_4[1] * x2) /
                    (TAN_DEGREE_4[2] + x2 * (TAN_DEGREE_4[3] + x2));
        }
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
// forward declaractions for LUT generation
template <typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr T sin(T x) noexcept requires(std::is_floating_point_v<T>);
template <typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr auto sin(T x) noexcept requires (std::is_integral_v<T>);
//================================== Interface ===============================//
/*
 * Approximations of Sine/Cosine, where:
 * polyApprox - defines implementation, poly for polynomial, else LUT
 * accuracy - integer value in range 0..10 that scales up the accuracy
 * of the approximation, where value stands for number of digits of accurary required
 * in fractional part of the result (i.e. to get 0.0xx accuracy choose accuracy=1,
 * so it guarantees that maximum absolute error will be lower than 0.1 on the argument range).
 * Better accuracy, leads to slower runtime (obviously)
*/
template <typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr T cos(T x) noexcept requires(std::is_floating_point_v<T>)
{
    static_assert (accuracy < SIN_COS_ACC_MAP_COUNT, "invalid accuracy");
    if (x == std::numeric_limits<T>::infinity()) // don't try to compute inf and signal a nan
        return std::numeric_limits<T>::signaling_NaN();
    const auto range = polyApprox ? HALF_PI : QUARTER_PI;
    const auto invRange = polyApprox ? INV_HALF_PI : INV_QUARTER_PI;
    const _Internal::ReductionRes res = x > RANGE_REDUCTION_SWITCH ? _Internal::payneHayekRangeReduce(x) : _Internal::addRangeReduce(x, range, invRange);
    if constexpr (polyApprox)
    {
        if (res.noReduciton)
            return _Internal::cos_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        const int sign = res.quad >= 0 ? 1 : -1;
        x *= sign;
        // split function period into 4 equal parts shifted by Pi/2
        switch ((res.quad*sign) & Pi3by2_2Pi)
        {
        case Zero_Pi2:
            return _Internal::cos_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        case Pi2_Pi:
            return -_Internal::sin_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        case Pi_Pi3by2:
            return -_Internal::cos_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        case Pi3by2_2Pi:
            return _Internal::sin_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        }
    }
    else
    {
        const int sign = x >= 0 ? 1 : -1; // get sign and negate it, so table won't out of range
        x *= sign;
        if (res.noReduciton)
            return _Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::cos,SIN_COS_FOLDING_RATIO, accuracy>>(x);
        const int quad = res.quad >= 0 ? res.quad : -res.quad;
        // split function period into 8 equal parts shifted by Pi/4
        switch (quad & PiPi3by4_2Pi)
        {
        case Zero_Pi4:
            return _Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::cos,SIN_COS_FOLDING_RATIO, accuracy>>(x);
        case Pi4_Pi2:
            return _Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::sin,SIN_COS_FOLDING_RATIO, accuracy>>(QUARTER_PI - x);
        case Pi2_Pi3by4:
            return -_Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::sin,SIN_COS_FOLDING_RATIO, accuracy>>(x);
        case Pi3by4_Pi:
            return -_Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::cos,SIN_COS_FOLDING_RATIO, accuracy>>(QUARTER_PI - x);
        case PiZero_Pi4:
            return -_Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::cos,SIN_COS_FOLDING_RATIO, accuracy>>(x);
        case PiPi4_Pi2:
            return -_Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::sin,SIN_COS_FOLDING_RATIO, accuracy>>(QUARTER_PI - x);
        case PiPi2_Pi3by4:
            return _Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::sin,SIN_COS_FOLDING_RATIO, accuracy>>(x);
        case PiPi3by4_2Pi:
            return _Internal::generic_inner_table<T,_Internal::LUTInfo<double,Trigonometrix::cos,SIN_COS_FOLDING_RATIO, accuracy>>(QUARTER_PI - x);
        }
    }
    assert(false && "invalid range");
    return x;
}

// wrapper function to handle interger arguments
template<typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr auto cos(T x) noexcept requires (std::is_integral_v<T>)
{
    return cos<double,accuracy,polyApprox>(double(x));
}

template <typename T, std::size_t accuracy, bool polyApprox>
constexpr T sin(T x) noexcept requires(std::is_floating_point_v<T>)
{
    static_assert (accuracy < SIN_COS_ACC_MAP_COUNT, "invalid accuracy");
    if (x == std::numeric_limits<T>::infinity())
        return std::numeric_limits<T>::signaling_NaN();
    if constexpr (polyApprox)
    {
        const _Internal::ReductionRes res = x > RANGE_REDUCTION_SWITCH ? _Internal::payneHayekRangeReduce(x) : _Internal::addRangeReduce(x, HALF_PI, INV_HALF_PI);
        if (res.noReduciton)
            return _Internal::sin_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        const int sign = res.quad >= 0 ? 1 : -1;
        x *= sign;
        switch ((res.quad*sign) & Pi3by2_2Pi)
        {
        case Zero_Pi2:
            return sign*_Internal::sin_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        case Pi2_Pi:
            return sign*_Internal::cos_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        case Pi_Pi3by2:
            return -sign*_Internal::sin_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        case Pi3by2_2Pi:
            return -sign*_Internal::cos_inner_polinomial<T,SIN_COS_ACC_MAP[accuracy]>(x);
        }
    }
    else
    {
        return cos<T,accuracy,polyApprox>(HALF_PI - x);
    }
    assert(false && "invalid range");
    return x;
}

// wrapper function to handle interger arguments
template <typename T, std::size_t accuracy, bool polyApprox>
constexpr auto sin(T x) noexcept requires (std::is_integral_v<T>)
{
    return sin<double,accuracy,polyApprox>(double(x));
}

template <typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr T cosDeg(T degrees) noexcept requires(std::is_floating_point_v<T>)
{
    if (degrees == std::numeric_limits<T>::infinity())
        return std::numeric_limits<T>::signaling_NaN();
    return cos<T,accuracy,polyApprox>(degToRad(degrees));
}

template <typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr T sinDeg(T degrees) noexcept requires(std::is_floating_point_v<T>)
{
    if (degrees == std::numeric_limits<T>::infinity())
        return std::numeric_limits<T>::signaling_NaN();
    return sin<T,accuracy,polyApprox>(degToRad(degrees));
}

// wrapper function to handle interger arguments
template <typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr auto cosDeg(T degrees) noexcept requires(std::is_integral_v<T>)
{
    return cos<double,accuracy,polyApprox>(degToRad(double(degrees)));
}

// wrapper function to handle interger arguments
template <typename T, std::size_t accuracy = sinCosAcc<T>, bool polyApprox = true>
constexpr auto sinDeg(T degrees) noexcept requires(std::is_integral_v<T>)
{
    return sin<double,accuracy,polyApprox>(degToRad(double(degrees)));
}

/*
 * Polynomial approximation of Tangent(x), where:
 * there are two modes, fast has max relative error of 0.0033 and slow - 1e-7,
 * although runtime is about the same, with slight advantage of the former.
 */
template <typename T, bool fast = true>
constexpr T tan(T x) noexcept requires(std::is_floating_point_v<T>)
{
    const int sign = x >= 0 ? 1 : -1;
    const _Internal::ReductionRes res = x > RANGE_REDUCTION_SWITCH ? _Internal::payneHayekRangeReduce(x) : _Internal::addRangeReduce(x, QUARTER_PI, INV_QUARTER_PI);
    x *= sign;
    if (x == 0 && (res.quad == 2*sign || res.quad == 6*sign))// result approaches infinity for args Pi/2 and 3Pi/2
        return std::numeric_limits<T>::infinity();
    switch ((res.quad*sign) & PiPi3by4_2Pi)
    {
    case Zero_Pi4:
         return sign*_Internal::tan_inner_polynomial<T,fast>(x);
    case Pi4_Pi2:
        return sign/_Internal::tan_inner_polynomial<T,fast>(QUARTER_PI - x);
    case Pi2_Pi3by4:
        return -sign/_Internal::tan_inner_polynomial<T,fast>(x);
    case Pi3by4_Pi:
        return -sign*_Internal::tan_inner_polynomial<T,fast>(QUARTER_PI - x);
    case PiZero_Pi4:
        return sign*_Internal::tan_inner_polynomial<T,fast>(x);
    case PiPi4_Pi2:
        return sign/_Internal::tan_inner_polynomial<T,fast>(QUARTER_PI - x);
    case PiPi2_Pi3by4:
        return -sign/_Internal::tan_inner_polynomial<T,fast>(x);
    case PiPi3by4_2Pi:
        return -sign*_Internal::tan_inner_polynomial<T,fast>(QUARTER_PI - x);
    }
    assert(false && "invalid range");
    return x;
}

// wrapper function to handle interger arguments
template <typename T, bool fast = true>
constexpr auto tan(T x) noexcept requires(std::is_integral_v<T>)
{
    return tan<double,fast>(double(x));
}

template <typename T, bool fast = true>
constexpr auto tanDeg(T degrees) noexcept requires(std::is_floating_point_v<T>)
{
    return tan<T,fast>(degToRad(degrees));
}

// wrapper function to handle interger arguments
template <typename T, bool fast = true>
constexpr auto tanDeg(T degrees) noexcept requires(std::is_integral_v<T>)
{
    return tan<double,fast>(degToRad(double(degrees)));
}


/*
 * Polynomial approximation for ArcTangent(x), returns atan in range [-pi/2,pi/2].
 * It has two modes - fast but lower max accuracy, and vice versa, although this affects
 * only calculations with args lower than switch value, where approximation switches to linear.
 * Fast version has max. absolute error - 0.014, slow version - 0.009, basically
 * former has 1 digit precision in fractional part and latter has 2 digits, but
 * also almost two times faster than slow version in the region before approx. switch.
 */
template <typename T, bool fast = true>
constexpr T atan(T x) noexcept requires(std::is_floating_point_v<T>)
{
    if (x == 0)
        return x;
    const int sign = x >= 0 ? 1 : -1;
    x *= sign;
    if constexpr(fast)// use lower degree polynomial
    {
        if (x > ATAN_APPROX_SWITCH_DEGREE_3)
            return sign * std::min(HALF_PI, ATAN_LINEAR_DEGREE_3_A * x + ATAN_LINEAR_DEGREE_3_B);
        else
        {
            constexpr auto polySize = ATAN_DEGREE_3.size();
            T res = ATAN_DEGREE_3[polySize-1];
            for (int i = polySize-2; i >= 0; --i)
                res = res * x + ATAN_DEGREE_3[i];

            return sign * res;
        }
    }
    else
    {
        if (x > ATAN_APPROX_SWITCH_DEGREE_8)
            return sign * std::min(HALF_PI, ATAN_LINEAR_DEGREE_8_A * x + ATAN_LINEAR_DEGREE_8_B);
        else
        {
            constexpr auto polySize = ATAN_DEGREE_8.size();
            T res = ATAN_DEGREE_8[polySize-1];
            for (int i = polySize-2; i >= 0; --i)
                res = res * x + ATAN_DEGREE_8[i];

            return sign * res;
        }
    }
}

// wrapper function to handle interger arguments
template <typename T, bool fast = true>
constexpr auto atan(T x) noexcept requires(std::is_integral_v<T>)
{
    return atan<double,fast>(x);
}


}
