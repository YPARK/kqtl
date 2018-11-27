// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "convergence.hh"
#include "options.hh"
#include "parameters.hh"
#include "rcpp_util.hh"
#include "regression.hh"
#include "regression_factored.hh"
#include "residual.hh"
#include "sgvb_regression_inference.hh"
#include "tuple_util.hh"
#include "kqtl_model.hh"

#ifndef RCPP_KQTL_HH_
#define RCPP_KQTL_HH_

using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
using Vec = Eigen::Matrix<float, Eigen::Dynamic, 1>;
using SpMat = Eigen::SparseMatrix<float, Eigen::ColMajor>;
using Scalar = Mat::Scalar;
using Index = Mat::Index;

#include "rcpp_kqtl_util.hh"
#include "rcpp_kernel_util.hh"
#include "rcpp_zscore_util.hh"
#include "rcpp_kqtl_regression.hh"

#endif
