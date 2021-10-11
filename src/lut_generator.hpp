#pragma once

#include <iostream>
#include <ios>
#include <fstream>
#include <cmath>
#include <cassert>
#include <iomanip>
#include "trigonometry.hpp"
#include <vector>

inline constexpr auto FLOAT_FNAME = "float_table.hpp";
inline constexpr auto DOUBLE_FNAME = "double_table.hpp";
inline constexpr auto FT_SIN_NAME = "SIN_TABLE_F";
inline constexpr auto DT_SIN_NAME = "SIN_TABLE_D";
inline constexpr auto FT_SIN_GRAD_NAME = "SIN_GRAD_F";
inline constexpr auto DT_SIN_GRAD_NAME = "SIN_GRAD_D";
inline constexpr auto FT_COS_NAME = "COS_TABLE_F";
inline constexpr auto DT_COS_NAME = "COS_TABLE_D";
inline constexpr auto FT_COS_GRAD_NAME = "COS_GRAD_F";
inline constexpr auto DT_COS_GRAD_NAME = "COS_GRAD_D";
inline constexpr int FLOAT_PREC_DIGITS = 12;
inline constexpr int DOUBLE_PREC_DIGITS = 19;

using namespace std;


class FileRedirectStream
{
public:
    FileRedirectStream(const char* filename, int precision) :
      fileStream(filename),
      cout_buff(cout.rdbuf())
    {
        assert(fileStream && "can't open file");
        cout.rdbuf(fileStream.rdbuf());
        flags = cout.flags();
        cout << hexfloat << setprecision(precision);
        fileStream << "#pragma once\n\n";
    };
    ~FileRedirectStream()
    {
        if (fileStream)
           fileStream.close();
        cout.rdbuf(cout_buff);
        cout.flags(flags);
    }

    // for regular output of variables and so on
    template<typename T> FileRedirectStream& operator<<(const T& something)
    {
        cout << something;
        return *this;
    }
    // for manipulators like std::endl
    typedef std::ostream& (*stream_function)(std::ostream&);
    FileRedirectStream& operator<<(stream_function func)
    {
        func(cout);
        return *this;
    }

private:
    ofstream fileStream;
    basic_streambuf<char>* cout_buff;
    ios_base::fmtflags flags;
};

static int tableSizeFromAcc(double relError)
{
    return int(M_PI / acos(1 - relError) / 8) + 1;
}

template<typename T, bool isSin>
static void writeTable(FileRedirectStream& file, int size, const char* countStr)
{
    const char* decl = "inline constexpr ";
    file << decl;
    if constexpr (is_same<T, float>())
    {
        file << "float ";
        if constexpr (isSin)
            file << FT_SIN_NAME;
        else
            file << FT_COS_NAME;
    }
    else
    {
        file << "double ";
        if constexpr (isSin)
            file << DT_SIN_NAME;
        else
            file << DT_COS_NAME;
    }
    file << "[" << countStr << ']' << " = {\n";
    T startArg = 0.;
    T currentArg = startArg;
    T step = QUARTER_PI / size;
    vector<T> table;
    table.resize(size);
    for(T i = 0; i < size; ++i)
    {
        T val;
        if constexpr (isSin)
            val = sin(currentArg);
        else
            val = cos(currentArg);
        table[i] = val;
        currentArg += step;
        file << val << ",\n";
    }
    file << "};\n\n";

    // gradient values
    file << decl;
    if constexpr (is_same<T, float>())
    {
        file << "float ";
        if constexpr (isSin)
            file << FT_SIN_GRAD_NAME;
        else
            file << FT_COS_GRAD_NAME;
    }
    else
    {
        file << "double ";
        if constexpr (isSin)
            file << DT_SIN_GRAD_NAME;
        else
            file << DT_COS_GRAD_NAME;
    }
    file << '[' << countStr << ']' << " = {\n";
    T prevVal = table[0];
    for(int i = 1; i < size-1; ++i)
    {
        file << (table[i] - prevVal) << ",\n";
        prevVal = table[i];
    }
    file << "};\n";
}

template<typename T>
static void generateSinCosTable()
{
    static_assert(is_floating_point<T>() && (is_same<T, float>() || is_same<T, double>()));
    int size;
    if constexpr (is_same<T, float>())
        size = tableSizeFromAcc(1E-9);
    else
    size = tableSizeFromAcc(1E-11);// should be 1E-17 actually, but it would be too enormous
    assert(size > 0 && "invalid table size");
    const char* fileName;
    int precision;
    if constexpr (is_same<T, float>())
    {
        fileName = FLOAT_FNAME;
        precision = FLOAT_PREC_DIGITS;
    }
    else
    {
        fileName = DOUBLE_FNAME;
        precision = DOUBLE_PREC_DIGITS;
    }
    FileRedirectStream file(fileName, precision);
    const int strBufSize = 8;
    char sizeStr[strBufSize];
    sprintf(sizeStr, "%d", size);
    const char* countStr;
    if constexpr (is_same<T, float>())
        countStr = "TABLE_COUNT_F";
    else
        countStr = "TABLE_COUNT_D";
    const char* decl = "inline constexpr ";
    file << decl << "int " << countStr << " = " << sizeStr << ";\n\n";
    writeTable<T, true>(file, size, countStr);
    file << '\n';
    writeTable<T, false>(file, size, countStr);
}
