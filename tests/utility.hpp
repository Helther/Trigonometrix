#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <functional>


//============================== Accuracy suite ============================//

// enchanced Kahan summation algorithm for negating loss of significance
template<typename T> T preciseSum(const std::vector<T>& input)
{
    static_assert (std::is_floating_point<T>(), "Useful only with floating-point numbers");
    T sum = 0;
    T correction = 0;
    if (input.size() == 0)
        return sum;

    for(const T& i : input)
    {
        T newSum = sum + i;
        if (std::abs(sum) >= std::abs(i))
            correction += (sum - newSum) + i; // compensate lower digits of input
        else
            correction += (i - newSum) + sum; // compensate lower digits of sum
        sum = newSum;
    }

    return sum + correction;
}

template<typename T> T absoluteError(const T& measure, const T& control)
{
    return measure - control;
}

template<typename T> T relativeError(const T& measure, const T& control, bool* exception = nullptr)
{
    if (std::abs(control) == 0)
    {
        if (exception)
            *exception = true;
        return 0;
    }

    return std::abs(1 - std::abs(measure / control));
}

// Root-Mean-Square error
template<typename T> T rmsError(const std::vector<T>& measure, const std::vector<T>& control)
{
    static_assert (std::is_integral<T>() || std::is_floating_point<T>(), "NaN");
    size_t count = measure.size();
    assert(count != 0 && count == control.size() && "invalid arrays sizes!");
    T divided = 0;
    if (std::is_floating_point<T>())
    {// use float optimized sum then
        std::vector<T> sum;
        for(size_t i = 0; i < count; ++i)
            sum.push_back(std::pow(measure[i] - control[i], 2));
        divided = preciseSum(sum);

    }
    else
    {// fallback
        for(size_t i = 0; i < count; ++i)
            divided += (std::pow(measure[i] - control[i], 2));
    }

    return std::sqrt(divided / count);
}

template<typename T> T absoluteMaxError(const std::vector<T>& measure, const std::vector<T>& control)
{
    assert(measure.size() == control.size() && "arrays sizes differ!");
    T maximum = 0;
    for(size_t i = 0; i < measure.size(); ++i)
        if (std::abs(measure[i] - control[i]) > maximum)
            maximum = measure[i] - control[i];

    return maximum;
}

template<typename T> T absoluteAverageError (const std::vector<T>& measure, const std::vector<T>& control)
{
    static_assert (std::is_integral<T>() || std::is_floating_point<T>(), "NaN");
    size_t count = measure.size();
    assert(count != 0 && count == control.size() && "invalid arrays sizes!");
    T divided = 0;
    if (std::is_floating_point<T>())
    {// use float optimized sum then
        std::vector<T> sum;
        for(size_t i = 0; i < count; ++i)
            sum.push_back(measure[i] - control[i]);
        divided = preciseSum(sum);
    }
    else
    {// fallback
        for(size_t i = 0; i < count; ++i)
            divided += (measure[i] - control[i]);
    }

    return divided / count;
}

template<typename T> T relativeMaxError(const std::vector<T>& measure, const std::vector<T>& control)
{
    assert(measure.size() == control.size() && "arrays sizes differ!");
    T maximum = 0;
    for(size_t i = 0; i < measure.size(); ++i)
        if (control[i] != 0)
            maximum = std::max(relativeError(measure[i], control[i]), maximum);

    return maximum;
}

template<typename T> T relativeAverageError (const std::vector<T>& measure, const std::vector<T>& control)
{
    static_assert (std::is_integral<T>() || std::is_floating_point<T>(), "NaN");
    size_t count = measure.size();
    assert(count != 0 && count == control.size() && "invalid arrays sizes!");
    T divided = 0;
    size_t currentCount = count;
    if (std::is_floating_point<T>())
    {// use float optimized sum then
        std::vector<T> sum;
        for(size_t i = 0; i < count; ++i)
        {
            bool except = false;
            const T res = relativeError(measure[i], control[i], &except);
            if (except)
                --currentCount;
            else
                sum.push_back(res);
        }
        divided = preciseSum(sum);
    }
    else
    {// fallback
        for(size_t i = 0; i < count; ++i)
        {
            bool except = false;
            const T res = relativeError(measure[i], control[i], &except);
            if (except)
                --currentCount;
            else
                divided += res;
        }
    }

    return divided / currentCount;
}

//============================ END Accuracy suite ===========================//


//============================ Benchmarking suite ===========================//

template<typename T>
void accuracyBench(T start, T range, T step, T(*MeasureFunc)(T), T(*ControlFunc)(T), const char* name)
{
    std::cout << std::endl <<"=========== Accuracy Benchmark " << name << " ===========" << std::endl;
    std::vector<T> measure;
    std::vector<T> control;
    for(T i = start; i < range; i+=step)
    {
        measure.push_back(MeasureFunc(i));
        control.push_back(ControlFunc(i));
    }
    std::cout << name << std::endl;
    std::cout << "abs err: " << absoluteAverageError(measure, control) << std::endl;
    std::cout << "rel err: " << relativeAverageError(measure, control) << std::endl;
    std::cout << "rms err: " << rmsError(measure, control) << std::endl;
    std::cout << "max abs err: " << absoluteMaxError(measure, control) << std::endl;
    std::cout << "max rel err: " << relativeMaxError(measure, control) << std::endl;
}
// TODO add args for choo distribution
template<typename T>
void accuracyBenchRand(T min, T max, int size, T(*MeasureFunc)(T), T(*ControlFunc)(T), std::random_device& rnd, const char* name)
{
    assert(min < max && "invalid args");
    std::cout << std::endl <<"=========== Accuracy Benchmark Random uniform " << name << " ============" << std::endl;
    std::vector<T> measure;
    std::vector<T> control;
    measure.resize(size);
    control.resize(size);
    for(int i = 0; i < size; ++i)
    {
        std::uniform_real_distribution<T> dist(min, max);
        T arg = dist(rnd);
        measure[i] = (MeasureFunc(arg));
        control[i] = (ControlFunc(arg));
    }
    std::cout << name << std::endl;
    std::cout << "abs err: " << absoluteAverageError(measure, control) << std::endl;
    std::cout << "rel err: " << relativeAverageError(measure, control) << std::endl;
    std::cout << "rms err: " << rmsError(measure, control) << std::endl;
    std::cout << "max abs err: " << absoluteMaxError(measure, control) << std::endl;
    std::cout << "max rel err: " << relativeMaxError(measure, control) << std::endl;
}
template<typename T>
void speedBench(T start, T range, T step, T(*MeasureFunc)(T), T(*ControlFunc)(T), const char* name)
{
    std::cout << std::endl <<"=========== Speed Benchmark " << name << " ============" << std::endl;
    auto startTime = std::chrono::high_resolution_clock::now();
    for(T i = start; i < range; i+=step)
        MeasureFunc(i);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Trigonometrix duration in mS: " << duration.count() << std::endl;
    startTime = std::chrono::high_resolution_clock::now();
    for(T i = start; i < range; i+=step)
        ControlFunc(i);
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Std duration in mS: " << duration.count() << std::endl;
}

template<typename T>
void speedBenchRand(T min, T max, int size, T(*MeasureFunc)(T), T(*ControlFunc)(T), std::random_device& rnd, const char* name)
{
    assert(min < max && "invalid args");
    std::cout << std::endl <<"=========== Speed Benchmark Random uniform " << name << " ============" << std::endl;
    std::vector<T> data;
    data.resize(size);
    for(int i = 0; i < size; ++i)
    {
        std::uniform_real_distribution<T> dist(min, max);
        data[i] = dist(rnd);
    }
    auto startTime = std::chrono::high_resolution_clock::now();
    for(T& i : data)
    MeasureFunc(i);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Trigonometrix duration in mS: " << duration.count() << std::endl;
    startTime = std::chrono::high_resolution_clock::now();
    for(T& i : data)
    ControlFunc(i);
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Std duration in mS: " << duration.count() << std::endl;
}
