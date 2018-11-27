################################################################
#' Compute Normalized Graph Laplacian matrix
#'
#' @export
#' @name normalized.laplacian
#'
#' @usage normalized.laplacian(W, tau = 1)
#'
#' @param W weighted adjacency matrix (n x n)
#' @param tau degree-correction parameter (1 x 1)
#'
#' @return Normalized Laplacian Matrix
#'
#' @examples
#'
#' W = abs(matrix(rnorm(3*3), 3, 3))
#' W = W %*% t(W)
#'
#' L = kqtl::normalized.laplacian(W, tau = 1)
#' print(L)
#'
#' diag(W) <- 0
#' d = apply(W, 2, sum)
#' d.tau.inv.sqrt = diag(1/sqrt(d + mean(d)))
#'
#' .L = - d.tau.inv.sqrt  %*% W %*% d.tau.inv.sqrt
#' .L = .L + diag(3)
#' print(.L)
#' 
normalized.laplacian <- function(W, tau = 1) {

    options <- list(laplcian.tau = tau)

    return(.Call('rcpp_laplacian', W, options, PACKAGE = 'kqtl'))
}


################################################################
#' Compute Eigen decomposition of the Graph Diffusion Kernel
#'
#' @export
#' @name normalized.gdk.eigen
#'
#' @usage normalized.gdk.eigen(W, tau = 1, beta = .5)
#'
#' @param W weighted adjacency matrix (n x n)
#' @param tau degree-correction parameter (0, infinity)
#' @param beta rate parameter for the diffusion [0, 1]
#'
#' @return eigen values (D2) and vectors (Vt) where Kernel K = t(Vt)*diag(D2)*Vt
#'
#' @examples
#'
#' W = abs(matrix(rnorm(3*3), 3, 3))
#' W = W %*% t(W)
#'
#' eig.out = kqtl::normalized.gdk.eigen(W, tau = 1, beta = .5)
#' print(eig.out)
#' 
#' Vt = sweep(eig.out$vectors, 1, sqrt(eig.out$values), `*`)
#' K = t(Vt) %*% Vt # graph diffusion kernel
#' print(K)
#' 
normalized.gdk.eigen <- function(W, tau = 1, beta = .5) {

    stopifnot(beta < 1)

    options <- list(laplcian.tau = tau, gdk.beta = beta)

    return(.Call('rcpp_gdk_eigen', W, options, PACKAGE = 'kqtl'))
}




################################################################
#' Variational inference of kernelized QTL model
#'
#' @export
#' @name fit.kqtl
#'
#' @usage fit.kqtl(effect, effect.se, kernel.eigen)
#'
#' @param effect Marginal effect size matrix (SNP x trait)
#' @param effect.se Marginal effect size standard error matrix (SNP x trait)
#' @param kernel.eigen Eigen decomposition of Kernel matrix
#' @param multi.C multivariate SNP confounding factors (SNP x confounder; default: NULL)
#' @param univar.C univariate SNP confounding factors (SNP x confounder; default: NULL)
#' @param factored Fit factored QTL model (default: FALSE)
#' @param options A list of inference/optimization options.
#' @param do.hyper Hyper parameter tuning (default: FALSE)
#' @param do.rescale Rescale z-scores by standard deviation (default: FALSE)
#' @param tau Fixed value of tau
#' @param pi Fixed value of pi
#' @param tau.lb Lower-bound of tau (default: -10)
#' @param tau.ub Upper-bound of tau (default: -4)
#' @param pi.lb Lower-bound of pi (default: -4)
#' @param pi.ub Upper-bound of pi (default: -1)
#' @param tol Convergence criterion (default: 1e-4)
#' @param gammax Maximum precision (default: 1000)
#' @param rate Update rate (default: 1e-2)
#' @param decay Update rate decay (default: 0)
#' @param jitter SD of random jitter for mediation & factorization (default: 0.1)
#' @param nsample Number of stochastic samples (default: 10)
#' @param vbiter Number of variational Bayes iterations (default: 2000)
#' @param verbose Verbosity (default: TRUE)
#' @param k Rank of the factored model (default: 1)
#' @param svd.init initialize by SVD (default: TRUE)
#' @param right.nn non-negativity in factored effect (default: FALSE)
#' @param mu.min mininum non-negativity weight (default: 0.01)
#' @param print.interv Printing interval (default: 10)
#' @param nthread Number of threads during calculation (default: 1)
#' @param out.residual estimate residual z-scores (default: FALSE)
#' @param do.var.calc variance calculation (default: FALSE)
#' @param nboot Number of bootstraps followed by finemapping (default: 0)
#' @param nboot.var Number of bootstraps for variance estimation (default: 10)
#' @param scale.var Scaled variance calculation (default: TRUE)
#' @param min.se Minimum level of SE (default: 1e-4)
#' @param rseed Random seed
#'
#'
#' @return a list of variational inference results
#'
#' @author Yongjin Park, \email{ypp@@csail.mit.edu}, \email{yongjin.peter.park@@gmail.com}
#'
#' @details
#'
#' @examples
#'
fit.kqtl <- function(effect,              # marginal effect : y ~ x
                     effect.se,           # marginal se : y ~ x
                     kernel.eigen,        # eigen decomposition
                     multi.C = NULL,      # covariate matrix (before LD)
                     univar.C = NULL,     # covariate (multified by LD)
                     factored = FALSE,    # Factored multiple traits
                     options = list(),
                     do.hyper = FALSE,
                     do.rescale = FALSE,
                     tau = NULL,
                     pi = NULL,
                     tau.lb = -10,
                     tau.ub = -4,
                     pi.lb = -4,
                     pi.ub = -1,
                     tol = 1e-4,
                     gammax = 1e3,
                     rate = 1e-2,
                     decay = 0,
                     jitter = 1e-1,
                     nsample = 10,
                     vbiter = 2000,
                     verbose = TRUE,
                     k = 1,
                     svd.init = TRUE,
                     right.nn = FALSE,
                     mu.min = 1e-2,
                     print.interv = 10,
                     nthread = 1,
                     out.residual = FALSE,
                     do.var.calc = FALSE,
                     nboot = 0,
                     nboot.var = 10,
                     scale.var = TRUE,
                     min.se = 1e-4,
                     rseed = NULL) {

    stopifnot(is.matrix(effect))
    stopifnot(is.matrix(effect.se))
    stopifnot(all(dim(effect) == dim(effect.se)))

    V.t <- kernel.eigen$vectors
    D2 <- matrix(as.numeric(kernel.eigen$values), ncol = 1)

    stopifnot(is.matrix(V.t))
    stopifnot(nrow(V.t) == nrow(D2))
    stopifnot(ncol(V.t) == nrow(effect))

    ## SNP confounding factors
    if(is.null(multi.C)) {
        p <- dim(effect)[1]
        multi.C <- matrix(1/p, p, 1)
    }

    if(is.null(univar.C)) {
        p <- dim(effect)[1]
        univar.C <- matrix(1/p, p, 1)
    }

    stopifnot(is.matrix(multi.C))
    stopifnot(dim(effect)[1] == dim(multi.C)[1])

    stopifnot(is.matrix(univar.C))
    stopifnot(dim(effect)[1] == dim(univar.C)[1])

    ## Override options
    opt.vars <- c('do.hyper', 'do.rescale', 'tau', 'pi', 'tau.lb',
                  'tau.ub', 'pi.lb', 'pi.ub', 'tol', 'gammax',
                  'rate', 'decay',
                  'jitter', 'nsample', 'vbiter', 'verbose',
                  'k', 'svd.init', 'right.nn', 'mu.min',
                  'print.interv', 'nthread',
                  'out.residual', 'min.se',
                  'rseed', 'do.var.calc', 'scale.var',
                  'nboot', 'nboot.var')

    .eval <- function(txt) eval(parse(text = txt))
    for(v in opt.vars) {
        val <- .eval(v)
        if(!(v %in% names(options)) && !is.null(val)) {
            options[[v]] <- val
        }
    }
    options[['sample.size']] <- n

    ## call R/C++ functions ##
    if(factored) {
        return(.Call('rcpp_fac_kqtl', effect, effect.se, V.t, D2, multi.C, univar.C, options, PACKAGE = 'kqtl'))
    } else {
        return(.Call('rcpp_kqtl', effect, effect.se, V.t, D2, multi.C, univar.C, options, PACKAGE = 'kqtl'))
    }
}

