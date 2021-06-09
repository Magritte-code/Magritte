#pragma once


#include <vector>
using std::vector;
#include <list>
using std::list;
#include <string>
using std::string;
using std::to_string;
#include <iostream>
using std::cout;
using std::endl;

#include "paracabs.hpp"
namespace pc = paracabs;

// Default Real and Size types
typedef double   Real;
typedef uint32_t Size;

using Vector3D = pc::datatypes::Vector3D <double>;

using AcceleratorThreads = pc::accelerator::AcceleratorThreads;
using        HostThreads = pc::multi_threading::HostThreads;

template <typename type>
using Vector = pc::datatypes::Vector<type>;
template <typename type>
using Matrix = pc::datatypes::Matrix<type>;
template <typename type>
using Tensor = pc::datatypes::Tensor<type>;
template <typename type, typename XThreads>
using TP = pc::datatypes::TP<type, XThreads>;
template <typename type, typename XThreads>
using VectorTP = pc::datatypes::VectorTP<type, XThreads>;
template <typename type, typename XThreads>
using MatrixTP = pc::datatypes::MatrixTP<type, XThreads>;


const Real ONE  = 1.0;
const Real TWO  = 2.0;
const Real HALF = 0.5;


// Vectors of Size
typedef vector<Size>  Size1;
typedef vector<Size1> Size2;
typedef vector<Size2> Size3;
typedef vector<Size3> Size4;

// Vectors of Real
typedef vector<Real>  Real1;
typedef vector<Real1> Real2;
typedef vector<Real2> Real3;
typedef vector<Real3> Real4;

// Vectors of bool
typedef vector<bool>  Bool1;
typedef vector<Bool1> Bool2;
typedef vector<Bool2> Bool3;

// Vectors of char
typedef vector<char>  Char1;
typedef vector<Char1> Char2;
typedef vector<Char2> Char3;

// Vectors of int
typedef vector<int>  Int1;
typedef vector<Int1> Int2;
typedef vector<Int2> Int3;

// Vectors of size_t
typedef vector<size_t>  Size_t1;
typedef vector<Size_t1> Size_t2;
typedef vector<Size_t2> Size_t3;
typedef vector<Size_t3> Size_t4;

// Vectors of long
typedef vector<long>  Long1;
typedef vector<Long1> Long2;
typedef vector<Long2> Long3;
typedef vector<Long3> Long4;

// Vectors of double
typedef vector<double>  Double1;
typedef vector<Double1> Double2;
typedef vector<Double2> Double3;
typedef vector<Double3> Double4;
typedef vector<Double4> Double5;

// Vectors of string
typedef vector<string>   String1;
typedef vector<String1>  String2;
typedef vector<String2>  String3;
