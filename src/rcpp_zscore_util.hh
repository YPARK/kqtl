#include "rcpp_util.hh"

#ifndef RCPP_ZSCORE_UTIL_HH_
#define RCPP_ZSCORE_UTIL_HH_

///////////////////////////////////////////////////
// standardize z-scores w.r.t. covariance matrix //
// R = V * D^2 * V'                              //
///////////////////////////////////////////////////

template <typename Derived, typename Derived2, typename Derived3>
Mat standardize_zscore(const Eigen::MatrixBase<Derived>& _zscore,
                       const Eigen::MatrixBase<Derived2>& Vt,
                       const Eigen::MatrixBase<Derived3>& D) {
  Mat Z = _zscore;
  Mat Y = D.cwiseInverse().asDiagonal() * Vt * Z;
  Mat xx = D.asDiagonal() * Vt * Mat::Ones(Vt.cols(), 1);
  Mat rr(Y.rows(), 1);
  Scalar xx_sum = xx.cwiseProduct(xx).sum();
  // Scalar n = Z.rows();
  Scalar denom = Y.rows();

  for (Index k = 0; k < Z.cols(); ++k) {
    Scalar xy = Y.col(k).cwiseProduct(xx).sum();
    Scalar mu = xy / xx_sum;
    rr = Y.col(k) - xx * mu;
    Scalar tau = rr.cwiseProduct(rr).sum() / denom + 1e-8;

    // TLOG("Standardize mu : " << mu << " xy : " << xy << " tau : " << tau
    //                          << " denom : " << denom);

    Z.col(k) = Vt.transpose() * D.asDiagonal() * rr / std::sqrt(tau);
  }

  return Z;
}

//////////////////////////////////////////////
// center z-scores w.r.t. covariance matrix //
// R = V * D^2 * V'                         //
//////////////////////////////////////////////

template <typename Derived, typename Derived2, typename Derived3>
Mat center_zscore(const Eigen::MatrixBase<Derived>& _zscore,
                  const Eigen::MatrixBase<Derived2>& Vt,
                  const Eigen::MatrixBase<Derived3>& D) {
  Mat Z = _zscore;
  Mat Y = D.cwiseInverse().asDiagonal() * Vt * Z;
  Mat xx = D.asDiagonal() * Vt * Mat::Ones(Vt.cols(), 1);
  Scalar xx_sum = xx.cwiseProduct(xx).sum();
  // Scalar denom = Y.rows();

  for (Index k = 0; k < Z.cols(); ++k) {
    Scalar xy = Y.col(k).cwiseProduct(xx).sum();
    Scalar mu = xy / xx_sum;
    // TLOG("Center mu : " << mu << " xy : " << xy);
    Z.col(k) = Vt.transpose() * D.asDiagonal() * (Y.col(k) - xx * mu);
  }

  return Z;
}

#endif
