#ifndef RCPP_KQTL_UTIL_HH_
#define RCPP_KQTL_UTIL_HH_


/////////////////
// set options //
/////////////////

void set_options_from_list(Rcpp::List& _list, options_t& opt) {
  if (_list.containsElementNamed("tau.lb"))
    opt.TAU_LODDS_LB = Rcpp::as<Scalar>(_list["tau.lb"]);

  if (_list.containsElementNamed("tau.ub"))
    opt.TAU_LODDS_UB = Rcpp::as<Scalar>(_list["tau.ub"]);

  if (_list.containsElementNamed("pi.lb"))
    opt.PI_LODDS_LB = Rcpp::as<Scalar>(_list["pi.lb"]);

  if (_list.containsElementNamed("pi.ub"))
    opt.PI_LODDS_UB = Rcpp::as<Scalar>(_list["pi.ub"]);

  if (_list.containsElementNamed("tau")) {
    opt.TAU_LODDS_LB = Rcpp::as<Scalar>(_list["tau"]);
    opt.TAU_LODDS_UB = Rcpp::as<Scalar>(_list["tau"]);
  }
  if (_list.containsElementNamed("pi")) {
    opt.PI_LODDS_LB = Rcpp::as<Scalar>(_list["pi"]);
    opt.PI_LODDS_UB = Rcpp::as<Scalar>(_list["pi"]);
  }
  if (_list.containsElementNamed("tol"))
    opt.VBTOL = Rcpp::as<Scalar>(_list["tol"]);
  if (_list.containsElementNamed("gammax"))
    opt.GAMMAX = Rcpp::as<Scalar>(_list["gammax"]);
  if (_list.containsElementNamed("decay"))
    opt.DECAY = Rcpp::as<Scalar>(_list["decay"]);
  if (_list.containsElementNamed("rate"))
    opt.RATE0 = Rcpp::as<Scalar>(_list["rate"]);
  if (_list.containsElementNamed("adam.m"))
    opt.RATE_M = Rcpp::as<Scalar>(_list["adam.m"]);
  if (_list.containsElementNamed("adam.v"))
    opt.RATE_V = Rcpp::as<Scalar>(_list["adam.v"]);
  if (_list.containsElementNamed("nsample"))
    opt.NSAMPLE = Rcpp::as<Index>(_list["nsample"]);

  if (_list.containsElementNamed("nboot"))
    opt.NBOOT = Rcpp::as<Index>(_list["nboot"]);

  if (_list.containsElementNamed("nboot.var"))
    opt.NBOOT_VAR = Rcpp::as<Index>(_list["nboot.var"]);

  if (_list.containsElementNamed("scale.var"))
    opt.SCALE_VAR_CALC = Rcpp::as<bool>(_list["scale.var"]);

  if (_list.containsElementNamed("num.duplicate.sample"))
    opt.N_DUPLICATE_SAMPLE = Rcpp::as<Index>(_list["num.duplicate.sample"]);

  if (_list.containsElementNamed("num.strat.size"))
    opt.N_STRAT_SIZE = Rcpp::as<Index>(_list["num.strat.size"]);

  if (_list.containsElementNamed("nsubmodel"))
    opt.N_SUBMODEL_MED = Rcpp::as<Index>(_list["nsubmodel"]);

  if (_list.containsElementNamed("num.submodel"))
    opt.N_SUBMODEL_MED = Rcpp::as<Index>(_list["num.submodel"]);

  if (_list.containsElementNamed("submodel.size"))
    opt.N_SUBMODEL_SIZE = Rcpp::as<Index>(_list["submodel.size"]);

  if (_list.containsElementNamed("print.interv"))
    opt.INTERV = Rcpp::as<Index>(_list["print.interv"]);

  if (_list.containsElementNamed("print.interval"))
    opt.INTERV = Rcpp::as<Index>(_list["print.interval"]);

  if (_list.containsElementNamed("nthread"))
    opt.NTHREAD = Rcpp::as<Index>(_list["nthread"]);
  if (_list.containsElementNamed("num.thread"))
    opt.NTHREAD = Rcpp::as<Index>(_list["num.thread"]);
  if (_list.containsElementNamed("k")) opt.K = Rcpp::as<Index>(_list["k"]);
  if (_list.containsElementNamed("K")) opt.K = Rcpp::as<Index>(_list["K"]);
  if (_list.containsElementNamed("re.k"))
    opt.RE_K = Rcpp::as<Index>(_list["re.k"]);
  if (_list.containsElementNamed("RE.K"))
    opt.RE_K = Rcpp::as<Index>(_list["RE.K"]);
  if (_list.containsElementNamed("vbiter"))
    opt.VBITER = Rcpp::as<Index>(_list["vbiter"]);
  if (_list.containsElementNamed("verbose"))
    opt.VERBOSE = Rcpp::as<bool>(_list["verbose"]);
  if (_list.containsElementNamed("random.effect"))
    opt.WITH_RANDOM_EFFECT = Rcpp::as<bool>(_list["random.effect"]);

  if (_list.containsElementNamed("ld.matrix"))
    opt.WITH_LD_MATRIX = Rcpp::as<bool>(_list["ld.matrix"]);

  if (_list.containsElementNamed("laplacian.tau"))
    opt.LAPLACIAN_TAU = Rcpp::as<bool>(_list["laplacian.tau"]);

  if (_list.containsElementNamed("gdk.beta"))
    opt.GDK_BETA = Rcpp::as<float>(_list["gdk.beta"]);

  if (_list.containsElementNamed("do.rescale"))
    opt.DO_RESCALE = Rcpp::as<bool>(_list["do.rescale"]);
  if (_list.containsElementNamed("rescale"))
    opt.DO_RESCALE = Rcpp::as<bool>(_list["rescale"]);
  if (_list.containsElementNamed("do.stdize"))
    opt.STD_LD = Rcpp::as<bool>(_list["do.stdize"]);
  if (_list.containsElementNamed("svd.init"))
    opt.MF_SVD_INIT = Rcpp::as<bool>(_list["svd.init"]);
  if (_list.containsElementNamed("mu.min"))
    opt.MU_MIN = Rcpp::as<Scalar>(_list["mu.min"]);
  if (_list.containsElementNamed("right.nn"))
    opt.MF_RIGHT_NN = Rcpp::as<bool>(_list["right.nn"]);
  if (_list.containsElementNamed("right.nonneg"))
    opt.MF_RIGHT_NN = Rcpp::as<bool>(_list["right.nonneg"]);
  if (_list.containsElementNamed("jitter"))
    opt.JITTER = Rcpp::as<Scalar>(_list["jitter"]);
  if (_list.containsElementNamed("rseed"))
    opt.RSEED = Rcpp::as<Scalar>(_list["rseed"]);
  if (_list.containsElementNamed("eigen.tol"))
    opt.EIGEN_TOL = Rcpp::as<Scalar>(_list["eigen.tol"]);
  if (_list.containsElementNamed("sample.size"))
    opt.SAMPLE_SIZE = Rcpp::as<Scalar>(_list["sample.size"]);
  if (_list.containsElementNamed("med.sample.size"))
    opt.M_SAMPLE_SIZE = Rcpp::as<Scalar>(_list["med.sample.size"]);
  if (_list.containsElementNamed("med.lodds.cutoff"))
    opt.MED_LODDS_CUTOFF = Rcpp::as<Scalar>(_list["med.lodds.cutoff"]);

  if (_list.containsElementNamed("do.hyper"))
    opt.DO_HYPER = Rcpp::as<bool>(_list["do.hyper"]);

  if (_list.containsElementNamed("do.finemap.unmed"))
    opt.DO_FINEMAP_UNMEDIATED = Rcpp::as<bool>(_list["do.finemap.unmed"]);

  if (_list.containsElementNamed("do.finemap.direct"))
    opt.DO_FINEMAP_UNMEDIATED = Rcpp::as<bool>(_list["do.finemap.direct"]);

  if (_list.containsElementNamed("do.var.calc"))
    opt.DO_VAR_CALC = Rcpp::as<bool>(_list["do.var.calc"]);

  if (_list.containsElementNamed("do.direct.estimation"))
    opt.DO_DIRECT_EFFECT = Rcpp::as<bool>(_list["do.direct.estimation"]);

  if (_list.containsElementNamed("do.control.backfire"))
    opt.DO_CONTROL_BACKFIRE = Rcpp::as<bool>(_list["do.control.backfire"]);

  if (_list.containsElementNamed("do.med.two.step"))
    opt.DO_MED_TWO_STEP = Rcpp::as<bool>(_list["do.med.two.step"]);

  if (_list.containsElementNamed("de.propensity")) {
    opt.DO_DIRECT_EFFECT_PROPENSITY = Rcpp::as<bool>(_list["de.propensity"]);
    opt.DO_DIRECT_EFFECT_FACTORIZATION = false;
  }

  if (_list.containsElementNamed("de.factorization")) {
    opt.DO_DIRECT_EFFECT_FACTORIZATION =
        Rcpp::as<bool>(_list["de.factorization"]);
    opt.DO_DIRECT_EFFECT_PROPENSITY = false;
  }

  if (_list.containsElementNamed("factorization.model")) {
    opt.DE_FACTORIZATION_MODEL = Rcpp::as<int>(_list["factorization.model"]);
  }

  if (_list.containsElementNamed("out.resid"))
    opt.OUT_RESID = Rcpp::as<bool>(_list["out.resid"]);
  if (_list.containsElementNamed("out.residual"))
    opt.OUT_RESID = Rcpp::as<bool>(_list["out.residual"]);
  if (_list.containsElementNamed("multivar.mediator"))
    opt.MULTI_MED_EFFECT = Rcpp::as<bool>(_list["multivar.mediator"]);
}

////////////////////////////
// preprocessing z-scores //
////////////////////////////

template <typename Derived>
Mat calc_effect_sqrt(const Eigen::MatrixBase<Derived>& effect,
                     const Eigen::MatrixBase<Derived>& effect_se,
                     const Scalar sample_size) {
  const Scalar one_val = 1.0;
  Mat effect_sqrt = effect_se;
  if (sample_size < one_val) {
    WLOG("Ingoring summary-statistics sample size");
  } else {
    effect_sqrt = (effect.cwiseProduct(effect) / sample_size +
                   effect_se.cwiseProduct(effect_se))
                      .cwiseSqrt();
  }

  return effect_sqrt;
}

template <typename Derived>
std::tuple<Mat, Mat, Mat> preprocess_effect(
    const Eigen::MatrixBase<Derived>& _effect,
    const Eigen::MatrixBase<Derived>& _effect_se, const Scalar sample_size) {
  // 1. characterize NaN or infinite elements
  const Scalar zero_val = 0.0;
  const Scalar one_val = 1.0;
  is_obs_op<Mat> is_obs;
  Mat obs_mat =
      _effect.unaryExpr(is_obs).cwiseProduct(_effect_se.unaryExpr(is_obs));

  if (obs_mat.sum() < (obs_mat.rows() * obs_mat.cols())) {
    WLOG("Make sure all the effect size and standard errors are observed!");
    WLOG("Results may become biased due to the different level of missingness.")
  }

  Mat effect, effect_se;
  remove_missing(_effect, effect);
  remove_missing(_effect_se, effect_se);

  // effect_sqrt = sqrt(effect^2 /n + effect_se^2)
  // effect_z = effect / effect_sqrt
  // weight = 1/effect_sqrt

  Mat effect_sqrt =
      calc_effect_sqrt(effect, effect_se, sample_size).cwiseProduct(obs_mat);

  // 2. safely inverse
  auto safe_inverse = [&](const Scalar& x) {
    if (x <= zero_val) return zero_val;
    return one_val / x;
  };

  auto safe_division = [&](const Scalar& a, const Scalar& b) {
    if (b <= zero_val) return zero_val;
    return a / b;
  };

  Mat weight = effect_sqrt.unaryExpr(safe_inverse).cwiseProduct(obs_mat);
  Mat effect_z =
      effect.binaryExpr(effect_sqrt, safe_division).cwiseProduct(obs_mat);

  return std::make_tuple(effect_z, effect_sqrt, weight);
}

//////////////////////////////////
// output effect size to R list //
//////////////////////////////////

template <typename T>
Rcpp::List param_rcpp_list(const T& param) {
  return impl_param_rcpp_list(param, sgd_tag<T>());
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param, const tag_param_spike_slab) {
  return Rcpp::List::create(Rcpp::_["theta"] = mean_param(param),
                            Rcpp::_["theta.var"] = var_param(param),
                            Rcpp::_["lodds"] = log_odds_param(param));
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param, const tag_param_spike_gamma) {
  return Rcpp::List::create(Rcpp::_["theta"] = mean_param(param),
                            Rcpp::_["theta.var"] = var_param(param),
                            Rcpp::_["lodds"] = log_odds_param(param));
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param, const tag_param_mixture) {
  return Rcpp::List::create(Rcpp::_["theta"] = mean_param(param),
                            Rcpp::_["theta.var"] = var_param(param),
                            Rcpp::_["lodds"] = log_odds_param(param));
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param,
                                const tag_param_col_spike_slab) {
  return Rcpp::List::create(Rcpp::_["theta"] = mean_param(param),
                            Rcpp::_["theta.var"] = var_param(param),
                            Rcpp::_["lodds"] = log_odds_param(param));
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param,
                                const tag_param_col_spike_gamma) {
  return Rcpp::List::create(Rcpp::_["theta"] = mean_param(param),
                            Rcpp::_["theta.var"] = var_param(param),
                            Rcpp::_["lodds"] = log_odds_param(param));
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param, const tag_param_col_slab) {
  return Rcpp::List::create(Rcpp::_["theta"] = mean_param(param),
                            Rcpp::_["theta.var"] = var_param(param));
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param, const tag_param_col_slab_zero) {
  return Rcpp::List::create(Rcpp::_["theta.var"] = var_param(param));
}

template <typename T>
Rcpp::List impl_param_rcpp_list(const T& param, const tag_param_slab) {
  return Rcpp::List::create(Rcpp::_["theta"] = mean_param(param),
                            Rcpp::_["theta.var"] = var_param(param));
}

#endif
