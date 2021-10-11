#include "../src/trigonometry.hpp"
#include "utility.hpp"

//============================= constexpr tests ==============================//
inline constexpr auto test1 = Trigonometrix::sin<float,false>(float(M_PI));
inline constexpr auto test2 = Trigonometrix::sin<double,false>(M_PI);
inline constexpr auto test3 = Trigonometrix::sin<float,true>(float(M_PI));
inline constexpr auto test4 = Trigonometrix::sin<double,false>(M_PI);
inline constexpr auto test5 = Trigonometrix::cos<float,false>(float(M_PI));
inline constexpr auto test6 = Trigonometrix::cos<double,false>(M_PI);
inline constexpr auto test7 = Trigonometrix::cos<float,true>(float(M_PI));
inline constexpr auto test8 = Trigonometrix::cos<double,false>(M_PI);


//====================== accuracy and speed tests ============================//

inline constexpr auto runCount = 1000000;
inline constexpr auto rangeVal = 20000*M_PI;
inline constexpr auto stepVal = 0.01;
inline constexpr auto periodRange = 2*M_PI;

template <typename T>
void accuracyValuesSin()
{
    std::cout << std::endl <<"============== 0, Pi/4, Pi/2,...2*Pi Accuracy range test" << " ==============" << std::endl;
    T start = 0;
    T step = M_PI_4;
    accuracyBench(start,T(periodRange),step,&Trigonometrix::sinRad<T, false>, &std::sin,"table implementation");
    accuracyBench(start,T(periodRange),step,&Trigonometrix::sinRad<T, true>, &std::sin,"polynomial implementation");
}

template <typename T>
void accuracyValuesCos()
{
    std::cout << std::endl <<"============== 0, Pi/4, Pi/2,...2*Pi Accuracy range test" << " ==============" << std::endl;
    T start = 0;
    T step = M_PI_4;
    accuracyBench(start,T(periodRange),step,&Trigonometrix::cosRad<T, false>, &std::cos,"table implementation");
    accuracyBench(start,T(periodRange),step,&Trigonometrix::cosRad<T, true>, &std::cos,"polynomial implementation");
}

void accuracyRangeTestsSin(std::random_device& r)
{
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sinRad<float, false>, &std::sin,"float table implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, false>, &std::sin,"double table implementation");
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sinRad<float, true>, &std::sin,"float test polynomial implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, true>, &std::sin,"double test polynomial implementation");

    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sinRad<float, false>, &std::sin,r,"float RAND table implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sinRad<double, false>, &std::sin,r,"double RAND table implementation");
    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount,&Trigonometrix::sinRad<float, true>, &std::sin,r,"float test RAND polynomial implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sinRad<double, true>, &std::sin,r,"double test RAND polynomial implementation");
}

void speedTestsSin(std::random_device& r)
{
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sinRad<float, false>, &std::sin,"float test table implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, false>, &std::sin,"double test table implementation");
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sinRad<float, true>, &std::sin,"float test polynomial implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, true>, &std::sin,"double test polynomial implementation");

    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sinRad<float, false>, &std::sin,r,"float table implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sinRad<double, false>, &std::sin,r,"double table implementation");
    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sinRad<float, true>, &std::sin,r,"float polynomial implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sinRad<double, true>, &std::sin,r,"double polynomial implementation");
}

void accuracyRangeTestsCos(std::random_device& r)
{
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sinRad<float, false>, &std::sin,"float table implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, false>, &std::sin,"double table implementation");
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sinRad<float, true>, &std::sin,"float test polynomial implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, true>, &std::sin,"double test polynomial implementation");

    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sinRad<float, false>, &std::sin,r,"float RAND table implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sinRad<double, false>, &std::sin,r,"double RAND table implementation");
    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount,&Trigonometrix::sinRad<float, true>, &std::sin,r,"float test RAND polynomial implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sinRad<double, true>, &std::sin,r,"double test RAND polynomial implementation");
}

void speedTestsCos(std::random_device& r)
{
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sinRad<float, false>, &std::sin,"float test table implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, false>, &std::sin,"double test table implementation");
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sinRad<float, true>, &std::sin,"float test polynomial implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sinRad<double, true>, &std::sin,"double test polynomial implementation");

    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sinRad<float, false>, &std::sin,r,"float table implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sinRad<double, false>, &std::sin,r,"double table implementation");
    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sinRad<float, true>, &std::sin,r,"float polynomial implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sinRad<double, true>, &std::sin,r,"double polynomial implementation");
}

void polyAccuracyTests(bool sin)
{
    std::cout << std::endl <<"============== Polinomials Accuracy test by number of terms" << " ==============" << std::endl;
    double start = 0;
    double step = stepVal;
    // have to do it the old way
    if (sin)
    {
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,0>, &std::sin, std::string("sin, number of terms: " + std::to_string(2)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,1>, &std::sin, std::string("sin, number of terms: " + std::to_string(3)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,2>, &std::sin, std::string("sin, number of terms: " + std::to_string(4)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,3>, &std::sin, std::string("sin, number of terms: " + std::to_string(5)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,4>, &std::sin, std::string("sin, number of terms: " + std::to_string(6)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,5>, &std::sin, std::string("sin, number of terms: " + std::to_string(7)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,6>, &std::sin, std::string("sin, number of terms: " + std::to_string(8)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sinRad<double,true,7>, &std::sin, std::string("sin, number of terms: " + std::to_string(9)).c_str());
    }
    else
    {
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,0>, &std::cos,std::string("cos, number of terms: " + std::to_string(2)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,1>, &std::cos,std::string("cos, number of terms: " + std::to_string(3)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,2>, &std::cos,std::string("cos, number of terms: " + std::to_string(4)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,3>, &std::cos,std::string("cos, number of terms: " + std::to_string(5)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,4>, &std::cos,std::string("cos, number of terms: " + std::to_string(6)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,5>, &std::cos,std::string("cos, number of terms: " + std::to_string(7)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,6>, &std::cos,std::string("cos, number of terms: " + std::to_string(8)).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cosRad<double,true,7>, &std::cos,std::string("cos, number of terms: " + std::to_string(9)).c_str());
    }
}


int main()
{
    std::random_device r;
    std::cout << std::endl <<"============== Sine Benchmark " << " ==============" << std::endl;
    accuracyRangeTestsSin(r);
    accuracyValuesSin<float>();
    accuracyValuesSin<double>();
    speedTestsSin(r);
    polyAccuracyTests(true);
    std::cout << std::endl <<"============== END Sine Benchmark " << " ==============" << std::endl;
    std::cout << std::endl <<"============== Cosine Benchmark " << " ==============" << std::endl;
    accuracyRangeTestsCos(r);
    accuracyValuesCos<float>();
    accuracyValuesCos<double>();
    speedTestsCos(r);
    polyAccuracyTests(false);
    std::cout << std::endl <<"============== END Cosine Benchmark " << " ==============" << std::endl;
    return 0;
}
