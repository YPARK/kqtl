#include "tuple_util.hh"

#ifndef RCPP_KERNEL_UTIL_HH_
#define RCPP_KERNEL_UTIL_HH_

/////////////////////////////////////////////
// basic pipeline of kernel construction   //
//                                         //
// 1. input "design" or "adjacency" matrix //
// 2. output eigen decomposition           //
/////////////////////////////////////////////

////////////////////////////////////////////////////////
// Normalized Laplacian (Ng, Jordan)                  //
//                                                    //
//     L = I - D^(-1/2) * W * D^(-1/2)                //
//                                                    //
// Re-scaled Laplacian (Chaudhuri, Chung, Tsiatas)    //
//                                                    //
//     L(tau) = I - D(tau)^(-1/2) * W * D(tau)^(-1/2) //
//     D(tau) = D + tau I                             //
//                                                    //
// When tau -> 0, it's the same                       //
////////////////////////////////////////////////////////

std::tuple<Mat, Mat, Mat> normalized_laplacian(const Mat& W,
                                               const options_t& opt) {
  Mat ret = W;
  Mat Dleft;
  Mat Dright;

  if (W.cols() != W.rows()) {
    ELOG("The given matrix is not symmetric.");
    return std::tuple<Mat, Mat, Mat>(ret, Dleft, Dright);
  }

  // remove diagonal elements and construct diagonal matrix
  const Index n = W.rows();
  for (Index j = 0; j < n; ++j) {
    ret(j, j) = 0.0;
  }

  // calculate degree
  Dleft = ret.transpose() * Mat::Ones(n, 1);
  Dright = ret * Mat::Ones(n, 1);

  const Scalar tau = opt.laplacian_tau() * (Dleft.mean() + Dright.mean()) *
                     static_cast<Scalar>(0.5);

  const Scalar one_val = 1.0;

  auto safe_inv_sq = [&](const auto& x) {
    return std::sqrt(one_val / (x + tau));
  };

  ret = (Dleft.unaryExpr(safe_inv_sq).asDiagonal()) * ret;
  ret = ret * (Dright.unaryExpr(safe_inv_sq).asDiagonal());

  Mat eye(n, n);
  eye.setIdentity(n, n);
  ret = eye - ret;

  return std::tuple<Mat, Mat, Mat>(ret, Dleft, Dright);
}

/////////////////////////////////////////////////
// perform eigen decomposition of the X matrix //
/////////////////////////////////////////////////

std::tuple<Mat, Mat> do_eigen_decomp(const Mat& X, const options_t& opt) {
  Eigen::SelfAdjointEigenSolver<Mat> es(X);
  Mat S = es.eigenvalues();
  Mat Vt = es.eigenvectors().transpose();

  const Scalar TOL = opt.eigen_tol();

  // set negative eigenvalues and vectors to zero
  for (Index j = 0; j < Vt.rows(); ++j) {
    if (S(j) <= TOL) {
      if (opt.verbose()) WLOG("Ignoring small eigen values ... " << S(j));
      S(j) = 0.0;
      Vt.row(j) *= 0.0;
    }
  }

  return std::tuple<Mat, Mat>(S, Vt);
}

//////////////////////////////////////////////////////////////////
// Compute kernels in graph data                                //
//                                                              //
// 1. compute Laplacian                                         //
// 2. apply reglurization function to eigen values              //
//    (NOTE: lambda >= 0 due to PSD)                            //
//                                                              //
// Examples of regularization functions (Smola & Kondor)        //
//                                                              //
//   r(lambda) = 1/(1 + sig2 * lamda) --> regularized Laplacian //
//                                                              //
//   r(lambda) = exp(- scale * lambda) --> diffusion (scale<1/d)//
//                                                              //
//   r(lambda) = (a I - lambda)^p  --> p-step random walk (a>2) //
//                                                              //
//   r(lambda) = cos lambda * pi/4  --> inverse cosine          //
//////////////////////////////////////////////////////////////////

std::tuple<Mat, Mat, Mat, Mat> do_gdk_eigen(const Mat& W,
                                            const options_t& opt) {
  Mat L, Dleft, Dright;
  Mat S, Vt;

  std::tie(L, Dleft, Dright) = normalized_laplacian(W, opt);
  std::tie(S, Vt) = do_eigen_decomp(L, opt);

  const Scalar maxS = S.maxCoeff();
  const Scalar maxBeta = 1.0 / maxS;

  const Scalar beta = std::min(opt.gdk_beta(), maxBeta);
  const Scalar TOL = opt.eigen_tol();

  if (opt.gdk_beta() > maxBeta) {
    WLOG("Beta is too big for this graph: " << opt.gdk_beta()
                                            << " --> rescaled: " << maxBeta);
  }

  if (maxS < TOL) {
    WLOG("Underlying Laplacian may not be PSD");
  }

  auto diffusion_op = [&](const auto& lambda) {
    return std::exp(-beta * lambda);
  };

  Mat eS = S.unaryExpr(diffusion_op);

  return std::tuple<Mat, Mat, Mat, Mat>(eS, Vt, S, L);
}

/////////////////////////////////////////////////////////
// singular value decomposition to construct LD matrix //
/////////////////////////////////////////////////////////

template <typename Derived>
std::tuple<Mat, Mat, Mat> do_svd(const Eigen::MatrixBase<Derived>& X,
                                 const options_t& opt) {
  using RowVec = typename Eigen::internal::plain_row_type<Mat>::type;
  // Center X matrix and divided by sqrt(n-1)
  // covariance = X'X
  Mat Xsafe = X;
  if (opt.std_ld()) {
    standardize(Xsafe);
    if (opt.verbose()) TLOG("Standardized matrix");
  } else {
    center(Xsafe);
    if (opt.verbose()) TLOG("Centered matrix");
  }

  is_obs_op<Mat> obs_op;
  const RowVec num_obs = X.unaryExpr(obs_op).colwise().sum();

  for (Index j = 0; j < Xsafe.cols(); ++j) {
    const Scalar nj = std::max(static_cast<Scalar>(2.0), num_obs(j));
    Xsafe.col(j) = Xsafe.col(j) / std::sqrt(nj - 1.0);
  }

  const Scalar n = static_cast<Scalar>(Xsafe.rows());
  Xsafe *= n;  // prevent underflow
  const Scalar TOL = opt.eigen_tol();

  if (opt.verbose()) TLOG("Start SVD ... ");
  Eigen::JacobiSVD<Mat> svd;
  svd.setThreshold(TOL);
  svd.compute(Xsafe, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if (opt.verbose()) TLOG("Done SVD");

  // Walk through eigen spectrum and pick components
  Index num_comp = 0;
  Mat d2vec = svd.singularValues() / n;
  d2vec = d2vec.cwiseProduct(d2vec);
  for (num_comp = 0; num_comp < svd.matrixV().cols(); ++num_comp) {
    if (d2vec(num_comp) < TOL) break;
  }

  if (opt.verbose()) TLOG("Included number of components : " << num_comp);

  if (num_comp < 1) {
    ELOG("0 Number of SVD component!");
  }

  // Don't forget to rescale
  Mat dvec_out(num_comp, 1);
  dvec_out = svd.singularValues().head(num_comp) / n;
  Mat Vt(num_comp, Xsafe.cols());
  Vt = svd.matrixV().leftCols(num_comp).transpose();
  Mat U(Xsafe.rows(), num_comp);
  U = svd.matrixU().leftCols(num_comp);

  return std::tuple<Mat, Mat, Mat>(U, dvec_out, Vt);
}

#endif
