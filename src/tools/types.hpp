#pragma once

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;
#include <list>
using std::list;
#include <string>
using std::string;
using std::to_string;
#include <Eigen/Core>
using Eigen::VectorXd;
using Eigen::MatrixXd;
#include <cmath>
using std::sqrt;

#include "paracabs.hpp"
namespace pc = paracabs;

// Default Real and Size types
typedef long double Real;
typedef uint32_t   Size;

using Vector3D = pc::datatypes::Vector3D <double>;

template <typename type>
using Vector = pc::datatypes::Vector <type>;
template <typename type>
using Matrix = pc::datatypes::Matrix <type>;
template <typename type>
using Tensor = pc::datatypes::Tensor <type>;


const Real one       = 1.0;
const Real two       = 2.0;
const Real three     = 3.0;
const Real half      = 0.5;
const Real quart     = 0.25;
const Real sqrt2     = sqrt(2.0);
const Real sqrt3     = sqrt(3.0);
const Real inv_sqrt2 = 1.0/sqrt2;


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

//// Vectors of Eigen::Vector3d
//typedef vector<Vector3d>   Vector3d1;
//typedef vector<Vector3d1>  Vector3d2;
//typedef vector<Vector3d2>  Vector3d3;
//
// Vectors of Eigen::VectorXd
typedef vector<VectorXd>   VectorXd1;
typedef vector<VectorXd1>  VectorXd2;
typedef vector<VectorXd2>  VectorXd3;

typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorXr;
typedef vector<VectorXr>  VectorXr1;
typedef vector<VectorXr1> VectorXr2;
typedef vector<VectorXr2> VectorXr3;

// Vectors of Eigen::MatrixXd
typedef vector<MatrixXd>   MatrixXd1;
typedef vector<MatrixXd1>  MatrixXd2;
typedef vector<MatrixXd2>  MatrixXd3;

typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;
typedef vector<VectorXr>  MatrixXr1;
typedef vector<VectorXr1> MatrixXr2;
typedef vector<VectorXr2> MatrixXr3;
