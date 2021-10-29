# Trigonometrix
Constexpr C++ implementation of trigonometric functions approximations
## Features
 * constexpr time evaluated  functions
 * compile-time options for implementation(for speed vs accuracy tradeoffs)
 * lookup table generator
 * testing and benchmarking utilities
 * test cases for accuracy and speed benchmarks
## Description
Include header files and build in your project or run the provided tests.  
For sine and cosine there are two implementations to compare against each other: lookup-table-based and polynom-based.  
LUT is generated into source files given desired accuracy(which derives it's size). 
Sizes where chosen based on used type precision, with accuracy in mind (they are quite large, 
even though the function was folded over 1/8 of the period). Additionally gradients between data points 
were calculated for linear approximation to boost accuracy.  
The other implementation uses MiniMax Polynomial Approximations generated for different degrees of accuracy. 
The more accurate the longer it takes to compute. There is an option to set the desired accuracy at compile-time.  

Tangent, and arc-functions aren't suited as well as periodic sine/cosine for straight minimax polynomial approximation,
so they provide an option for fast/slow versions.
## Dependencies
1. A C++ compiler that supports C++20 standart.
The following compilers should work:

  * [gcc 9+](https://gcc.gnu.org/)

  * [clang 10+](https://clang.llvm.org/)

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
* accuracy tests for different degrees of polinomial approximation  

Results of test run perfomed with clang-11 and gcc-10 with -O3 option(w/o fast math optimization) 
are saved in [test_results.txt](test_results.txt)


In conclusion tests showed lower than expected accuracy for both implementations (likely due to rounding 
while reducing argument range, but it was expected for large args). Also perfomance is on par(in case of polynom) 
or slighly worse than STL, but it is also expected for such large tables, so further optimizations will come.
Other functions showed overall better perfomance compared to stl at a higher cost of accuracy.  
Test results snippet:
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
=========== Speed Benchmark acos ============
2.20 times faster than std
```

### Building the tests
ctests are build via CMake by specifying test target.
```
 $ cmake --build . --target SinCosTests
```

### Running the tests
To run built tests
```
 $ ctest -C Debug
```

## TODO: 
* Add comparison tests with other libraries
* add constexpr generation and usage of lookup tables of variable sizes
* Add tests for FMA contracted versions of approximations
* Add implementations using vector SSE intrinsics and corresponding tests 


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
