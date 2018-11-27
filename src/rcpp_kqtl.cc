#include "rcpp_kqtl.hh"

RcppExport SEXP rcpp_laplacian(SEXP W_sexp, SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const Mat>::type W(W_sexp);
  Rcpp::List options_list(options_sexp);

  options_t opt;
  set_options_from_list(options_list, opt);

  Mat L, Dleft, Dright;
  std::tie(L, Dleft, Dright) = normalized_laplacian(W, opt);

  return Rcpp::List::create(Rcpp::_["D.left"] = Dleft,
                            Rcpp::_["D.right"] = Dright, Rcpp::_["L"] = L);

  END_RCPP
}

RcppExport SEXP rcpp_gdk_eigen(SEXP W_sexp, SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const Mat>::type W(W_sexp);
  Rcpp::List options_list(options_sexp);
  options_t opt;
  set_options_from_list(options_list, opt);

  Mat eS, V, S, L;
  std::tie(eS, V, S, L) = do_gdk_eigen(W, opt);

  return Rcpp::List::create(Rcpp::_["values"] = eS, Rcpp::_["vectors"] = V,
                            Rcpp::_["L.values"] = S, Rcpp::_["L"] = L);

  END_RCPP
}

RcppExport SEXP rcpp_kqtl(SEXP effect_sexp, SEXP effect_se_sexp,  //
                          SEXP vt_sexp, SEXP d2_sexp,             //
                          SEXP c_sexp, SEXP c_delta_sexp,         //
                          SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const Mat>::type effect(effect_sexp);
  Rcpp::traits::input_parameter<const Mat>::type effect_se(effect_se_sexp);
  Rcpp::traits::input_parameter<const Mat>::type D2(d2_sexp);
  Rcpp::traits::input_parameter<const Mat>::type Vt(vt_sexp);
  Rcpp::traits::input_parameter<const Mat>::type C(c_sexp);
  Rcpp::traits::input_parameter<const Mat>::type Cdelta(c_delta_sexp);
  Rcpp::List options_list(options_sexp);

  options_t opt;
  set_options_from_list(options_list, opt);
  return Rcpp::wrap(impl_fit_kqtl(effect, effect_se, Vt, D2, C, Cdelta, opt));
  END_RCPP
}

RcppExport SEXP rcpp_fac_kqtl(SEXP effect_sexp, SEXP effect_se_sexp,  //
                              SEXP vt_sexp, SEXP d2_sexp,             //
                              SEXP c_sexp, SEXP c_delta_sexp,
                              SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const Mat>::type effect(effect_sexp);
  Rcpp::traits::input_parameter<const Mat>::type effect_se(effect_se_sexp);
  Rcpp::traits::input_parameter<const Mat>::type D2(d2_sexp);
  Rcpp::traits::input_parameter<const Mat>::type Vt(vt_sexp);
  Rcpp::traits::input_parameter<const Mat>::type C(c_sexp);
  Rcpp::traits::input_parameter<const Mat>::type Cdelta(c_delta_sexp);
  Rcpp::List options_list(options_sexp);

  options_t opt;
  set_options_from_list(options_list, opt);
  return Rcpp::wrap(
      impl_fit_fac_kqtl(effect, effect_se, Vt, D2, C, Cdelta, opt));
  END_RCPP
}
