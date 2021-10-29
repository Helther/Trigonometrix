#include "../src/trigonometry.hpp"
#include "utility.hpp"

//============================= constexpr tests ==============================//

template <typename T>
constexpr void constexprTest()
{// TODO add tan,atan,asin,acos
    constexpr auto test1 = Trigonometrix::sin<T,false>(T(0));
    constexpr auto test2 = Trigonometrix::sin<T,true>(T(0));
    constexpr auto test3 = Trigonometrix::cos<T,false>(T(0));
    constexpr auto test4 = Trigonometrix::cos<T,true>(T(0));
    constexpr auto test5 = Trigonometrix::sinDeg<T,false>(T(0));
    constexpr auto test6 = Trigonometrix::sinDeg<T,true>(T(0));
    constexpr auto test7 = Trigonometrix::cosDeg<T,false>(T(0));
    constexpr auto test8 = Trigonometrix::cosDeg<T,true>(T(0));
    constexpr auto test9 = Trigonometrix::tan<T,true>(T(0));
    constexpr auto test10 = Trigonometrix::tanDeg<T,true>(T(0));
    constexpr auto test11 = Trigonometrix::atan<T,true>(T(0));
    constexpr auto test12 = Trigonometrix::atan<T,false>(T(0));
    constexpr auto test13 = Trigonometrix::acos<T>(T(0));
    constexpr auto test14 = Trigonometrix::asin<T>(T(0));
    constexpr auto test15 = Trigonometrix::tan<T,false>(T(0));
    constexpr auto test116 = Trigonometrix::tanDeg<T,false>(T(0));
}


//====================== accuracy and speed tests ============================//

inline constexpr auto runCount = 1000000;
inline constexpr auto rangeVal = 5000*M_PI;
inline constexpr auto stepVal = 0.01;
inline constexpr auto periodRange = 2*M_PI;

template <typename T>
void accuracyValuesSin()
{
    const char* type;
    if constexpr(std::is_same_v<T,float>)
        type = "float";
    else
        type = "double";
    std::cout << std::endl <<"============== 0, Pi/4, Pi/2,...2*Pi Accuracy range test " << type << " ==============" << std::endl;
    T start = 0;
    T step = M_PI_4;
    accuracyBench(start,T(periodRange),step,&Trigonometrix::sin<T, false>, &std::sin,"table implementation");
    accuracyBench(start,T(periodRange),step,&Trigonometrix::sin<T, true>, &std::sin,"polynomial implementation");
}

template <typename T>
void accuracyValuesCos()
{
    const char* type;
    if constexpr(std::is_same_v<T,float>)
        type = "float";
    else
        type = "double";
    std::cout << std::endl <<"============== 0, Pi/4, Pi/2,...2*Pi Accuracy range test " << type << " ==============" << std::endl;
    T start = 0;
    T step = M_PI_4;
    accuracyBench(start,T(periodRange),step,&Trigonometrix::cos<T, false>, &std::cos,"table implementation");
    accuracyBench(start,T(periodRange),step,&Trigonometrix::cos<T, true>, &std::cos,"polynomial implementation");
}

void accuracyRangeTestsSin(std::random_device& r)
{
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sin<float, false>, &std::sin,"float table implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, false>, &std::sin,"double table implementation");
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sin<float, true>, &std::sin,"float test polynomial implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, true>, &std::sin,"double test polynomial implementation");

    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sin<float, false>, &std::sin,r,"float RAND table implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sin<double, false>, &std::sin,r,"double RAND table implementation");
    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount,&Trigonometrix::sin<float, true>, &std::sin,r,"float test RAND polynomial implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sin<double, true>, &std::sin,r,"double test RAND polynomial implementation");
}

void speedTestsSin(std::random_device& r)
{
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sin<float, false>, &std::sin,"float test table implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, false>, &std::sin,"double test table implementation");
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sin<float, true>, &std::sin,"float test polynomial implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, true>, &std::sin,"double test polynomial implementation");

    /* unreliable
    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sin<float, false>, &std::sin,r,"float table implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sin<double, false>, &std::sin,r,"double table implementation");
    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sin<float, true>, &std::sin,r,"float polynomial implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sin<double, true>, &std::sin,r,"double polynomial implementation");
    */
}

void accuracyRangeTestsCos(std::random_device& r)
{
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sin<float, false>, &std::sin,"float table implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, false>, &std::sin,"double table implementation");
    accuracyBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sin<float, true>, &std::sin,"float test polynomial implementation");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, true>, &std::sin,"double test polynomial implementation");

    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sin<float, false>, &std::sin,r,"float RAND table implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sin<double, false>, &std::sin,r,"double RAND table implementation");
    accuracyBenchRand(float(-rangeVal),float(rangeVal),runCount,&Trigonometrix::sin<float, true>, &std::sin,r,"float test RAND polynomial implementation");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::sin<double, true>, &std::sin,r,"double test RAND polynomial implementation");
}

void speedTestsCos(std::random_device& r)
{
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal), &Trigonometrix::sin<float, false>, &std::sin,"float test table implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, false>, &std::sin,"double test table implementation");
    speedBench(float(-rangeVal),float(rangeVal),float(stepVal),&Trigonometrix::sin<float, true>, &std::sin,"float test polynomial implementation");
    speedBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::sin<double, true>, &std::sin,"double test polynomial implementation");
    /* unreliable
    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sin<float, false>, &std::sin,r,"float table implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sin<double, false>, &std::sin,r,"double table implementation");
    speedBenchRand(float(-rangeVal),float(rangeVal),runCount, &Trigonometrix::sin<float, true>, &std::sin,r,"float polynomial implementation");
    speedBenchRand(-rangeVal,rangeVal,runCount, &Trigonometrix::sin<double, true>, &std::sin,r,"double polynomial implementation");
    */
}

void SinCosPolyAccuracyTests(bool sin)
{
    std::cout << std::endl <<"============== Polynomials Accuracy test by number of terms" << " ==============" << std::endl;
    double start = 0;
    double step = stepVal;
    // have to do it the old way
    if (sin)
    {
        accuracyBench(start,periodRange,step,&Trigonometrix::sin<double,0>, &std::sin, std::string("sin, number of terms: " +
                                                                                                        std::to_string(std::get<0>(SIN_POLIES[1]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sin<double,1>, &std::sin, std::string("sin, number of terms: " +
                                                                                                        std::to_string(std::get<0>(SIN_POLIES[2]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sin<double,2>, &std::sin, std::string("sin, number of terms: " +
                                                                                                        std::to_string(std::get<0>(SIN_POLIES[3]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sin<double,4>, &std::sin, std::string("sin, number of terms: " +
                                                                                                        std::to_string(std::get<0>(SIN_POLIES[4]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sin<double,6>, &std::sin, std::string("sin, number of terms: " +
                                                                                                        std::to_string(std::get<0>(SIN_POLIES[5]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sin<double,7>, &std::sin, std::string("sin, number of terms: " +
                                                                                                        std::to_string(std::get<0>(SIN_POLIES[6]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::sin<double,9>, &std::sin, std::string("sin, number of terms: " +
                                                                                                        std::to_string(std::get<0>(SIN_POLIES[7]))).c_str());
    }
    else
    {
        accuracyBench(start,periodRange,step,&Trigonometrix::cos<double,0>, &std::cos,std::string("cos, number of terms: " +
                                                                                                       std::to_string(std::get<0>(COS_POLIES[1]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cos<double,1>, &std::cos,std::string("cos, number of terms: " +
                                                                                                       std::to_string(std::get<0>(COS_POLIES[2]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cos<double,2>, &std::cos,std::string("cos, number of terms: " +
                                                                                                       std::to_string(std::get<0>(COS_POLIES[3]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cos<double,4>, &std::cos,std::string("cos, number of terms: " +
                                                                                                       std::to_string(std::get<0>(COS_POLIES[4]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cos<double,6>, &std::cos,std::string("cos, number of terms: " +
                                                                                                      std::to_string(std::get<0>(COS_POLIES[5]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cos<double,7>, &std::cos,std::string("cos, number of terms: " +
                                                                                                       std::to_string(std::get<0>(COS_POLIES[6]))).c_str());
        accuracyBench(start,periodRange,step,&Trigonometrix::cos<double,9>, &std::cos,std::string("cos, number of terms: " +
                                                                                                       std::to_string(std::get<0>(COS_POLIES[7]))).c_str());
    }
}

void tanTests(std::random_device& r)
{
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::tan<double,true>, &std::tan,"tan, fast version");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::tan<double,false>, &std::tan,"tan, slow version");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::tan<double, true>, &std::tan,r,"RAND tan, fast version");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::tan<double, false>, &std::tan,r,"RAND tan, slow version");
    speedBench<double,std::chrono::microseconds>(-rangeVal,rangeVal,stepVal,&Trigonometrix::tan<double,true>, &std::tan,"tan, fast version");
    speedBench<double,std::chrono::microseconds>(-rangeVal,rangeVal,stepVal,&Trigonometrix::tan<double,false>, &std::tan,"tan, slow version");
}

void atanTests(std::random_device& r)
{
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::atan<double,true>, &std::atan,"atan, fast version");
    accuracyBench(-rangeVal,rangeVal,stepVal,&Trigonometrix::atan<double,false>, &std::atan,"atan, slow version");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::atan<double, true>, &std::atan,r,"RAND atan, fast version");
    accuracyBenchRand(-rangeVal,rangeVal,runCount,&Trigonometrix::atan<double, false>, &std::atan,r,"RAND atan, slow version");
    std::cout << std::endl <<"============== Speed test inside higher degree polynomial arg range" << " ==============" << std::endl;
    speedBench<double,std::chrono::microseconds>(-double(ATAN_APPROX_SWITCH_DEGREE_3),ATAN_APPROX_SWITCH_DEGREE_3,0.0001,&Trigonometrix::atan<double,true>, &Trigonometrix::atan,"atan, fast version");
    speedBench<double,std::chrono::microseconds>(-double(ATAN_APPROX_SWITCH_DEGREE_8),ATAN_APPROX_SWITCH_DEGREE_8,0.0001,&Trigonometrix::atan<double,false>, &Trigonometrix::atan,"atan, slow version");
    std::cout << std::endl <<"============== Speed test for general arg range" << " ==============" << std::endl;
    speedBench<double,std::chrono::microseconds>(-rangeVal, rangeVal, stepVal,&Trigonometrix::atan<double,true>, &Trigonometrix::atan,"atan, fast version");
    speedBench<double,std::chrono::microseconds>(-rangeVal, rangeVal, stepVal,&Trigonometrix::atan<double,false>, &Trigonometrix::atan,"atan, slow version");

}

void asinTests(std::random_device& r)
{
    accuracyBench(-1.,1.,0.0001,&Trigonometrix::asin, &std::asin,"asin");
    accuracyBenchRand(-1.,1.,runCount,&Trigonometrix::asin, &std::asin,r,"RAND asin");
    speedBench<double,std::chrono::nanoseconds>(-1.,1.,0.0001,&Trigonometrix::asin, &Trigonometrix::asin,"asin");
}

void acosTests(std::random_device& r)
{
    accuracyBench(-1.,1.,0.0001,&Trigonometrix::acos, &std::acos,"acos");
    accuracyBenchRand(-1.,1.,runCount,&Trigonometrix::acos, &std::acos,r,"RAND acos");
    speedBench<double,std::chrono::nanoseconds>(-1.,1.,0.0001,&Trigonometrix::acos, &Trigonometrix::acos,"acos");
}


int main()
{
    const char* sep = "========================================================================";
    const char* sepBrackets = "========================";
    std::random_device r;
    constexprTest<float>();
    constexprTest<double>();
    constexprTest<int>();
    std::cout << std::endl << sep << std::endl << sepBrackets << " Sine Benchmark " << sepBrackets << std::endl;
    accuracyRangeTestsSin(r);
    accuracyValuesSin<float>();
    accuracyValuesSin<double>();
    speedTestsSin(r);
    SinCosPolyAccuracyTests(true);
    std::cout << std::endl << sep << std::endl << sepBrackets << " END Sine Benchmark " << sepBrackets << std::endl;

    std::cout << std::endl << sep << std::endl << sepBrackets << " Cosine Benchmark " << sepBrackets << std::endl;
    accuracyRangeTestsCos(r);
    accuracyValuesCos<float>();
    accuracyValuesCos<double>();
    speedTestsCos(r);
    SinCosPolyAccuracyTests(false);
    std::cout << std::endl << sep << std::endl << sepBrackets << " END Cosine Benchmark " << sepBrackets << std::endl;

    std::cout << std::endl << sep << std::endl << sepBrackets << " Tangent Benchmark " << sepBrackets << std::endl;
    tanTests(r);
    std::cout << std::endl << sep << std::endl << sepBrackets << " END Tangent Benchmark " << sepBrackets << std::endl;

    std::cout << std::endl << sep << std::endl << sepBrackets << " Arc Tangent Benchmark " << sepBrackets << std::endl;
    atanTests(r);
    std::cout << std::endl << sep << std::endl << sepBrackets << " END Arc Tangent Benchmark " << sepBrackets << std::endl;

    std::cout << std::endl << sep << std::endl << sepBrackets << " Arc Sine Benchmark " << sepBrackets << std::endl;
    asinTests(r);
    std::cout << std::endl << sep << std::endl << sepBrackets << " END Arc SineBenchmark " << sepBrackets << std::endl;

    std::cout << std::endl << sep << std::endl << sepBrackets << " Arc Cosine Benchmark " << sepBrackets << std::endl;
    acosTests(r);
    std::cout << std::endl << sep << std::endl << sepBrackets << " END Arc Cosine Benchmark " << sepBrackets << std::endl;
    return 0;
}
