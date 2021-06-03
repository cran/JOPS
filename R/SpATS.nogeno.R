#' Two-dimensional P-spline smoothing
#'
#' @description Two-dimensional smoothing of scattered data points with tensor product P-splines.
#'
#' @param response a character string with the name of the variable that contains the response variable of interest.
#' @param spatial a right hand \code{\link{formula}} object specifying the spatial P-Spline model. See \code{\link{SAP}} and \code{\link{PSANOVA}} for more details about how to specify the spatial trend.
#' @param fixed an optional right hand \code{\link{formula}} object specifying the fixed effects.
#' @param random an optional right hand \code{\link{formula}} object specifying the random effects. Currently, only sets of independent and identically distributed random effects can be incorporated.
#' @param data a data frame containing the variables.
#' @param family object of class \code{\link{family}} specifying the distribution and link function.
#' @param offset an optional numerical vector containing an a priori known component to be included in the linear predictor during fitting.
#' @param weights an optional numerical vector of weights to be used in the fitting process. By default, the weights are considered to be one.
#' @param control a list of control values.
#'
#' @return A list with the following components:
#' \item{call}{the matched call.}
#' \item{data}{the original supplied data argument with a new column with the weights used during the fitting process.}
#' \item{model}{a list with the model components: response, spatial, fixed and/or random.}
#' \item{fitted}{a numeric vector with the fitted values.}
#' \item{residuals}{a numeric vector with deviance residuals.}
#' \item{psi}{a two-length vector with the values of the dispersion parameters at convergence. For Gaussian responses both elements coincide, being the (REML) estimate of dispersion parameter. For non-Gaussian responses, the result depends on the argument \code{update.psi} of the \code{\link{controlSpATS}} function. If this argument was specified to \code{FALSE} (the default), the first component of the vector corresponds to the default value used for the dispersion parameter (usually 1). The second element, correspond to the (REML) estimate of the dispersion parameter at convergence. If the argument \code{update.psi} was specified to \code{TRUE}, both components coincide (as in the Gaussian case).}
#' \item{var.comp}{a numeric vector with the (REML) variance component estimates. This vector contains the variance components associated with the spatial trend, as well as those related with the random model terms.}
#' \item{eff.dim}{a numeric vector with the estimated effective dimension (or effective degrees of freedom) for each model component (spatial, fixed and/or random).}
#' \item{dim}{a numeric vector with the (model) dimension of each model component (spatial, fixed and/or random). This value corresponds to the number of parameters to be estimated.}
#' \item{dim.nom}{a numeric vector with the (nominal) dimension of each component (spatial, fixed and/or random). For the random terms of the model, this value corresponds to upper bound for the effective dimension (i.e., the maximum effective dimension a random term can achive). This nominal dimension is \eqn{rank[X, Z_k] - rank[X]}, where \eqn{Z_k} is the design matrix of the \eqn{k}th random factor and \eqn{X} is the design matrix of the fixed part of the model. In most cases (but not always), the nominal dimension corresponds to the model dimension minus one, ``lost'' due to the implicit constraint that ensures the mean of the random effects to be zero.}
#' \item{nobs}{number of observations used to fit the model.}
#' \item{niterations}{number of iterations EM-algorithm.}
#' \item{deviance}{the (REML) deviance at convergence (i.e., \eqn{-2} times the restricted log-likelihood).}
#' \item{coeff}{a numeric vector with the estimated fixed and random effect coefficients.}
#' \item{terms}{a list with the model terms: response, spatial, fixed and/or random. The information provided here is useful for printing and prediction purposes.}
#' \item{vcov}{inverse of the coefficient matrix of the mixed models equations. The inverse is needed for the computation of standard errors. For computational issues, the inverse is returned as a list: C22_inv corresponds to the coefficient matrix associated with the spatial, the fixed and the random components.}
#'
#' @details This function is a modified version of the function  \code{\link{SpATS}} in the package \code{SpATS}. The difference is that genotypes have been removed.
#'
#' @author  Maria-Xose Rodriguez-Alvarez and Paul Eilers
#'
#' @references Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics, 23, 52 - 71. https://doi.org/10.1016/j.spasta.2017.10.003.
#'
#' @export
#'
#' @examples
#' # Get the data
#' library(SemiPar)
#' data(ethanol)
#'
#' # Fit the PS-ANOVA model
#' ps2d <- SpATS.nogeno(response = "NOx",
#'                      spatial = ~PSANOVA(E, C,  nseg = c(20, 20), nest.div = c(2, 2)),
#'                      data = ethanol,
#'                      control = list(maxit = 100, tolerance = 1e-05,
#'                                     monitoring = 0, update.psi = FALSE))
#'
#' # Report effective dimensions, if desired
#' # print(summary(ps2d))
#'
#' # Compute component surface and their sum on a fine grid
#' Tr = obtain.spatialtrend(ps2d, grid = c(100, 100))
#'
#' # Plot surface and contours
#' image(Tr$row.p, Tr$col.p, Tr$fit, col = terrain.colors(100), xlab = 'C', ylab = 'E')
#' contour(Tr$row.p, Tr$col.p, Tr$fit, add = TRUE, col = 'blue')
#' points(ethanol$C, ethanol$E, pch = '+')

SpATS.nogeno <- function(response, spatial, fixed = NULL,
                         random = NULL, data, family = gaussian(),
                         offset = 0, weights = NULL,
                         control = list(maxit = 100)) {

	if (control$monitoring) start = proc.time()[3]

	weights <- as.vector(weights)
	if(is.null(weights)) weights = rep(1, nrow(data))
	if(length(offset) == 1) offset <- rep(offset, nrow(data))

	if(inherits(fixed, "character")) fixed <- as.formula(fixed)
	if(inherits(random, "character"))	random <- as.formula(random)
	if(inherits(spatial, "character")) spatial <- as.formula(spatial)

	sf <- interpret.SpATS.formula(spatial)

	# NAs in the covariates
	model.terms <- c(sf$x.coord, sf$y.coord,
    if(!is.null(fixed)) attr(terms.formula(fixed), "term.labels"),
    if(!is.null(random)) attr(terms.formula(random),"term.labels"))
	na.ind <- apply(is.na(data[,model.terms]), 1, any)
	na.pos <- (1:nrow(data))[!na.ind]
	weights <- weights*(!na.ind)*(!is.na(data[, response]))

	data.na <- data[!na.ind,]
	weights.na <- weights[!na.ind]
	offset.na <-  offset[!na.ind]

	# NAs in the response
	na.res <- is.na(data[,response])

	y <- data.na[,response]
	nobs <- length(y[weights.na != 0])

	# Construct design matrices
	MMns <- NULL
    init.var <- gg <- dim <- list()
	int <- rep(1, nrow(data))
    dim.int <- c(Intercept = 1)
    attr(dim.int, "random") <- FALSE
    attr(dim.int, "spatial") <- FALSE
    MMns <- cbind(MMns, Intercept = int)
    dim <- c(dim, list(dim.int))
    if (!is.null(fixed)) {
        fixed.part <- construct.fixed.part(formula = fixed, data = data)
        MMns <- cbind(MMns, fixed.part$X)
        dim <- c(dim, list(fixed.part$dim))
    }
    if (!is.null(random)) {
        random.part <- construct.random.part(formula = random,
            data = data)
        MMns <- cbind(MMns, random.part$Z)
        dim <- c(dim, list(random.part$dim))
        gg <- c(gg, list(random.part$g))
        init.var <- c(init.var, list(rep(0.001, length(random.part$init.var))))
    }
    spat.part <- construct.2d.pspline(formula = spatial, data = data, na.res = na.res)
    MMns <- cbind(MMns, spat.part$X, spat.part$Z)
    dim <- c(dim, list(spat.part$dim$fixed), list(spat.part$dim$random))
    gg <- c(gg, list(spat.part$g))
    init.var <- c(init.var, list(spat.part$init.var))

    g <- construct.capital.lambda(gg)
    random. <- unlist(lapply(dim, get.attribute, "random"))
    spatial. <- unlist(lapply(dim, get.attribute, "spatial"))
    spat.part$terms$fixed$pos <-
       create.position.indicator(unlist(dim), !random. & spatial.)
    spat.part$terms$random$pos <-
       create.position.indicator(unlist(dim), random. & spatial.)
    res <- list()
    res$fixed.part <- if (!is.null(fixed)){
        fixed.part
    } else{
    	NULL
    }
    res$random.part <- if (!is.null(random)){
        random.part
    } else {
    	NULL
    }
    res$spat.part <- spat.part
    res$terms$spatial <- spat.part$terms
    res$terms$fixed <- if (!is.null(fixed)){
        fixed.part$terms
    } else {
    	NULL
    }
    res$terms$random <- if (!is.null(random)){
        random.part$terms
    } else {
    	NULL
    }
    res$g <- g
    res$dim <- dim
    res$init.var <- init.var
    res$MM <- list(MMns = MMns)

	MM <- res
	ldim <- MM$dim
	random. <- unlist(lapply(ldim, get.attribute, "random"))
	spatial. <- unlist(lapply(ldim, get.attribute, "spatial"))
	dim <- unlist(ldim)
	g <- MM$g

	# Nominal dimension
	dim.nom <- obtain.nominal.dimension(MM$MM$MMns,
      dim, random., spatial., weights.na)

	random.pos <- create.position.indicator(dim, random.)
	fixed.pos <- create.position.indicator(dim, !random.)
	df.fixed <- sum(dim[!random.])

	X <- MMns[,fixed.pos,drop = FALSE]
	Z <- MMns[,random.pos,drop = FALSE]

	np <- c(ncol(X), ncol(Z))
	D <- diag(c(rep(0,np[1]), rep(1,sum(np[2]))))
	mustart <- etastart <- NULL
	eval(family$initialize)
	mu <- mustart
	eta <- family$linkfun(mustart)

	# Initialize variance components
	la <- c(1, unlist(MM$init.var))
	# Initialize deviance and psi
	devold <- 1e10
	psi <- la[1]

	if(control$monitoring > 1) {
		cat("Effective dimensions\n")
		cat("-------------------------\n")
		cat(sprintf("%1$3s %2$12s","It.","Deviance"), sep = "")
		cat(sprintf("%12s", names(g)), sep = "")
		cat("\n")
	}
	for (iter in 1:control$maxit) {
		deriv <- family$mu.eta(eta)
		z <- (eta - offset.na) + (y - mu)/deriv
		w <- as.vector(deriv ^ 2 / family$variance(mu))
		w <- w * weights.na
		z[!weights.na] <- 0

 		XtW. = t(X * w)
		XtX. = XtW. %*% X
		XtZ. = XtW. %*% Z
		ZtX. = t(XtZ.)
		ZtW. = t(Z*w)
		ZtZ. = ZtW. %*% Z
		Xty. = XtW. %*% z
		Zty. = ZtW. %*% z
		yty. <- sum(w * z ^ 2)
		ZtXtZ = rbind(XtZ., ZtZ.)
		u <- c(Xty., Zty.)

		if(control$monitoring) start1 <- proc.time()[3]
		for (it in 1:control$maxit) {
			# Build penalty matrix: block diagonal matrix
			Ginv <- vector(length = sum(dim[random.], na.rm = TRUE))
			for (i in 1:length(g)) {
				Ginv <- Ginv + (1/la[i+1])*g[[i]]
			}
			G = 1 / Ginv
      V <- rbind(cbind(XtX., t(ZtX. * G)),
                 cbind(ZtX., t(ZtZ. * G)))
			H <- (1 / la[1]) * V + D
			Hinv <- solve(H)
			b.aux <- (1 / la[1]) * Hinv %*% u

			b.fixed <- b.aux[1:np[1]]
			b.random <- G*b.aux[-(1:np[1])]

			b <- rep(0, sum(np))
			b[fixed.pos] <- b.fixed
			b[random.pos] <- b.random

			dZtPZ <- 1 / la[1] * apply((t(Hinv[-(1:np[1]),]) * ZtXtZ), 2, sum)
			ed <- tau <- vector(mode="list")
			for (i in 1:length(g)) {
				g.inv.d <- (1 / la[i+1]) * g[[i]]
				ed[[i]] <- sum(dZtPZ * (g.inv.d * G ^ 2))
				ed[[i]] <- ifelse(ed[[i]] == 0, 1e-50, ed[[i]])
				tau[[i]] <- sum(b.random ^ 2 * g[[i]]) / ed[[i]]
				tau[[i]] <- ifelse(tau[[i]] == 0, 1e-50, tau[[i]])
			}
			ssr = yty. - t(c(b.fixed, b.random)) %*% (2 * u - V %*% b.aux)
			# Compute deviance
			dev <- determinant(H)$modulus +
                sum(log(la[1] * 1 / w[w != 0])) +
                ssr / la[1] + sum(b.random ^ 2 * Ginv)
			psinew <- as.numeric((ssr/(nobs - sum(unlist(ed)) - df.fixed)))
			if(family$family == "gaussian" | control$update.psi) {
				psi2 <- psinew
			} else {
				psi2 <- 1
			}
			# New variance components and convergence check
			lanew <- c(psi2, unlist(tau))
			dla = abs(devold - dev)
			if(control$monitoring > 1) {
				cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
				cat(sprintf("%12.3f", unlist(ed)), sep = "")
				cat('\n')
			}
			if (dla < control$tolerance) break
			la = lanew
			psi = psinew
			devold = dev
		}
		if (control$monitoring) {
			end1 <- proc.time()[3]
			cat("Timings:\nSpATS", (end1-start1), "seconds\n")
		}
		eta.old <- eta
		eta <- MMns%*%b + offset.na
		mu <- family$linkinv(eta)
		# Convergence criterion: linear predictor
		tol <- sum((eta - eta.old)^2)/sum(eta^2)
		if (tol < control$tolerance | (family$family == "gaussian" & family$link== "identity")) break
	}
	var.comp <- la[-1]
	eff.dim <- unlist(ed)
	names(var.comp) <- names(eff.dim) <- names(g)

	# Effective dimension (fixed + random)
	eff.dim <- c(dim[!random.], eff.dim)

	attr(dim, "random") <- random.
	attr(dim, "spatial") <- spatial.

	fitted <- rep(NA, nrow(data))
	fitted[!na.ind] <- mu

	# Deviance residuals
	dev.residuals <- family$dev.resids(data[,response], fitted, weights)
    s <- attr(dev.residuals, "sign")
    if (is.null(s))
        s <- sign(data[,response] - fitted)
    dev.residuals <- sqrt(pmax(dev.residuals, 0))*s

	if (control$monitoring) {
		end <- proc.time()[3]
		cat("All process", (end - start), "seconds\n")
	}
	res <- list()
	res$call <- match.call()
	res$data <- cbind(data, weights = weights)
	res$model <- list(response = response, spatial = spatial, fixed = fixed, random = random)
	res$fitted <- fitted
	res$residuals <- dev.residuals
	res$psi <- c(la[1], psi)
	res$var.comp <- var.comp
	res$eff.dim <- eff.dim
	res$dim <- dim
	res$dim.nom <- dim.nom
	res$nobs <- nobs
	res$deviance <- dev
	res$coeff <- b
	random.coeff <- rep(FALSE, length(b))
	random.coeff[create.position.indicator(dim, random.)] <- TRUE
	attr(res$coeff, "random") <- random.coeff
	# Terms
	res$terms <- MM$terms
	class(res) <- "SpATS"
	res
}
