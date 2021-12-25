# Trigonometrix
Fast constexpr C++ implementation of trigonometric functions approximations
## Features
 * constexpr time evaluated  functions
 * compile-time options for implementation(for speed vs accuracy tradeoffs)
 * constexpr lookup table generator
 * testing and benchmarking utilities
 * test cases for accuracy and speed benchmarks
 * SSE2/FMA implementations of some functions(where it's appropriate)
## Description
Include header files and build in your project or run the provided tests.  
For sine and cosine there are two implementations to compare against each other: lookup-table-based and polynom-based.  
LUT is generated at compile-time, given desired accuracy (which derives it's size). 
Sizes where chosen based on desired accuracy (function was folded over 1/8 of the period, so as tables). Additionally gradients between data points 
were calculated for linear approximation to boost accuracy.  
The other implementation uses MiniMax Polynomial Approximations generated for different degrees of accuracy. 
The more accurate the longer it takes to compute. There is an option to set the desired accuracy at compile-time.  
There are also fast implementations using FMA instruction set to speed up polynom computation and reduce rounding error.  

Tangent, and arc-functions aren't suited as well as periodic sine/cosine for straight minimax polynomial approximation,
so they provide an option for fast/slow versions.
## Dependencies
1. A C++ compiler that supports C++20 standart.
The following compilers should work:

  * [gcc 10+](https://gcc.gnu.org/)

  * [clang 11+](https://clang.llvm.org/)

2. [CMake 3.5+](https://cmake.org/)

### Build
Target platforms are Windows and Linux. You can build main target via CMake
From build directory:
```
 $ cmake -DCMAKE_BUILD_TYPE=%TypeName% -DCMAKE_PREFIX_PATH=%QTDIR% -G "%Generator name%" .
 $ cmake --build . --target constexprTrig
``` 

## Testing
CMake's ctest was used for testing.
Following tests were conducted against the standart STL functions:
* accuracy tests over the large range of arguments for diffrent implementations over 1e+6+ runs
* speed tests over the large range of arguments for different implementations over 1e+6+ runs
* both random and fixed argument values versions of same tests
* accuracy tests for different degrees of polinomial approximation and table sizes

Results of test run perfomed with clang-11 and gcc-10 with -O3 option(w/o fast math optimization) 
are saved in [test_results.txt](test_results.txt)


In conclusion tests showed lower than expected accuracy for both implementations (likely due to rounding 
while reducing argument range, but it was expected for large args). Also perfomance is on par(in case of polynom) 
or slighly worse than STL in case of table implementation, so further optimizations are needed.
Other functions showed overall better perfomance compared to stl at a higher cost of lower accuracy.  
SSE benchmark results varied heavily between compilers, with clang somehow outperforming SSE version  
with two separate calls to sin and cos.  
Test results snippet(averages, because actual resuls may vary greatly depending on the compiler):
```
sine
=========== Speed Benchmark float test table implementation ============
1.03 times slower than std
=========== Speed Benchmark double test table implementation ============
1.72 times faster than std
=========== Speed Benchmark float test polynomial implementation ============
1.22 times slower than std
=========== Speed Benchmark double test polynomial implementation ============
1.68 times faster than std

cosine
=========== Speed Benchmark float test table implementation ============
1.03 times slower than std
=========== Speed Benchmark double test table implementation ============
1.77 times faster than std
=========== Speed Benchmark float test polynomial implementation ============
1.25 times slower than std
=========== Speed Benchmark double test polynomial implementation ============
1.65 times faster than std

tangent
=========== Speed Benchmark tan, fast version ============
5.06 times faster than std
=========== Speed Benchmark tan, slow version ============
4.41 times faster than std

arc tangent
=========== Speed Benchmark atan, fast version ============
1.58 times faster than std
=========== Speed Benchmark atan, slow version ============
1.73 times faster than std

arc sine
=========== Speed Benchmark asin ============
1.50 times faster than std

arc cosine
=========== Speed Benchmark acos ============
2.20 times faster than std


SSE enchanced functions
=========== Speed Benchmark sin ============
1.58 times faster than std
=========== Speed Benchmark cos ============
1.64 times faster than std
=========== Speed Benchmark sin-cos ============
3.03 times faster than std
```

### Building the tests
ctests are build via CMake by specifying test target.
Example:
```
 $ cmake --build . --target BenchmarkTests
```

### Running the tests
To run built tests
```
 $ ctest -C Debug
```

## TODO: 
* Add comparison tests with other libraries
* Add Unit tests for accuracy

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
