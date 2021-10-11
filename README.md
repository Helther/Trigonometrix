# constexprTrig
Constexpr C++ implementation of trigonometric functions approximations
## Features
 * constexpr time evaluated sine and cosine functions
 * two compile-time options for implementation
 * lookup table files generator
 * testing and benchmarking utilities
 * test cases for sine and cosine
## Description
Inlude header files and build in your project or run the provided tests.  
There are two implementations to compare agains each other: lookup-table-based and polynom-based.  
LUT is generated into source files given desired accuracy(which derives it's size). 
Sizes where chosen based on used type precision, with accuracy in mind (they are quite large, 
even though the function was folded over 1/8 of the period). Additionally gradients between data points 
were calculated for linear approximation to boost accuracy.  
The other implementation uses MiniMax Polynomial Approximations generated for different degrees of accuracy. 
The more accurate the longer it takes to compute. There is an option to set the desired accuracy at compile-time.
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
For testing was used CMake's ctest.
Following tests were conducted against the standart STL functions:
* accuracy tests over the large range of arguments for both implementations over 1e+6+ runs
* speed tests over the large range of arguments for both implementations over 1e+6+ runs
* both random and fixed argument values versions of same tests
* accuracy tests for different degrees of polinomial approximation  

Results of test run perfomed with clang++-11 with -O3 option(w/o fast math optimization) 
are saved in [test_results.txt](test_results.txt)


In conclusion tests showed lower than expected accuracy for both implementations (likely due to rounding 
while reducing argument range, but it was expected for large args). Also perfomance is on par(in case of polynom) 
or slighly worse than STL, but it is also expected for such large tables, so further optimizations will come.

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
* Update documentation
* Add comparison tests with other libraries
* add constexpr generation and usage of lookup tables of variable sizes
* Add tests for FMA contracted versions of approximations
* Add implementations using vector SSE intrinsics and corresponding tests 


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
