"factor.model.stat" <-
function (x, weights=seq(1, 3, length.out=nobs), output="full", center=TRUE, 
	frac.var=.5, iter.max=1, nfac.miss=1, full.min=20, reg.min=40, 
	sd.min=20, quan.sd=.90, tol=1e-3, zero.load=FALSE, 
	range.factors=c(0, Inf), constant.returns.okay=FALSE,
	specific.floor=0.1, floor.type="quantile", verbose=2)
{
        fun.copyright <- "Placed in the public domain 2006-2012 by Burns Statistics Ltd."
	fun.version <- "factor.model.stat 013"

	subfun.ssd <- function(z, weights, sd.min) {
		nas <- is.na(z)
		if(any(nas)) {
			if(sum(!nas) < sd.min) return(NA)
			sum(weights[!nas] * z[!nas]^2) / sum(weights[!nas])
		} else {
			sum(weights * z^2)
		}
	}

	#
	# start of main function
	#

	if(is.data.frame(x)) {
		x <- as.matrix(x)
	} else {
		if(!is.matrix(x)) stop("'x' needs to be a matrix")
	}
	if(!is.numeric(x)) stop("'x' needs to be numeric")
	x[!is.finite(x)] <- NA

	# for use in finance, try to check it is returns and not prices
	if(verbose >= 1 && min(x, na.rm=TRUE) >= 0) {
		warning(paste("minimum of values in 'x' is",
			min(x, na.rm=TRUE), "are you giving a price",
			"matrix rather than a return matrix?",
			"(warning suppressed if verbose < 1)"))
	}

	xna <- is.na(x)
	allmis <- rowSums(xna) == ncol(x)
	if(any(allmis)) {
		x <- x[!allmis, , drop=FALSE]
		xna <- is.na(x)
	}
	num.mis <- colSums(xna)

	if(any(num.mis > 0)) {
		if(sum(num.mis == 0) < full.min) 
			stop("not enough columns without missing values")
		if(!length(dimnames(x)[[2]])) 
			stop("'x' needs column names when missing values exist")

		max.miss <- max(num.mis)
		lnfm <- length(nfac.miss)
		if(lnfm == 0) stop("'nfac.miss' must have positive length")
		nfac.miss <- round(nfac.miss)
		if(any(nfac.miss < 0)) 
			stop("negative values in 'nfac.miss'")
		if(lnfm < max.miss) {
			nfac.miss <- c(nfac.miss, rep(nfac.miss[lnfm],
				max.miss - lnfm))
		}
	}

	if(!is.character(output) || length(output) != 1) {
                stop(paste("'output' should be a single character string",
                        "-- given has mode", mode(output), "and length",
                        length(output)))
        }
        output.menu <- c("full", "factor", "systematic", "specific")
        output.num <- pmatch(output, output.menu, nomatch=0)
        if(output.num == 0) {
                stop(paste("unknown or ambiguous input for 'output'",
                        "-- the allowed choices are:",
                        paste(output.menu, collapse=", ")))
        }
        output <- output.menu[output.num]

	nassets <- ncol(x)
	nobs <- nrow(x)
	if(!is.numeric(weights)) {
		stop(paste("'weights' must be numeric -- given has mode",
			mode(weights), "and length", length(weights)))
	}
	if(length(weights) != nobs) {
		if(length(weights) == nobs + sum(allmis)) {
			weights <- weights[!allmis]
		} else if(length(weights) == 1 && weights > 0) {
			weights <- rep(1, nobs)
		} else {
			stop(paste("bad value for 'weights'",
				"-- must be a single positive number",
				"(meaning equal weighting) or have length",
				"equal to the number of observations"))
		}
	}
	if(any(weights < 0)) {
		stop(paste(sum(weights < 0), "negative value(s) in 'weights'"))
	}
	weights <- weights / sum(weights)

	if(is.logical(center)) {
		if(center) {
			center <- colSums(x * weights, na.rm=TRUE)
		} else {
			center <- rep(0, nassets)
		}
	} else if(length(center) != nassets) stop("wrong length for 'center'")
	x <- sweep(x, 2, center, "-")

	sdev <- sqrt(apply(x, 2, subfun.ssd, weights=weights, sd.min=sd.min))
	if(any(sdev <= 0, na.rm=TRUE)) {
		sdzero <- !is.na(sdev) & sdev <= 0
		if(constant.returns.okay) {
		    sdev[sdzero] <- NA
		    if(verbose >= 1) {
			warning(paste(sum(sdzero),
				"asset(s) with constant returns:",
				paste(dimnames(x)[[2]][sdzero], collapse=", "),
				"(warning suppressed with verbose < 1)"))
		    }
		} else {
			stop(paste(sum(sdzero),
				"asset(s) with constant returns:",
				paste(dimnames(x)[[2]][sdzero], collapse=", ")))
		}
	}
	if(any(is.na(sdev))) {
		sdev[is.na(sdev)] <- quantile(sdev, quan.sd, na.rm=TRUE)
	}
	x <- scale(x, scale=sdev, center=FALSE)

	xw <- sqrt(weights) * x
	fullxw <- xw[, num.mis == 0, drop=FALSE]
	fx.svd <- svd(fullxw, nu=0)
	cumvar <- cumsum(fx.svd$d^2) / sum(fx.svd$d^2)
	nfac <- sum(cumvar < frac.var) + 1
	if(nfac > max(range.factors)) {
		nfac <- max(range.factors)
	} else if(nfac < min(range.factors)) {
		nfac <- min(range.factors)
	}
	if(nfac > length(cumvar)) nfac <- length(cumvar)
	fseq <- 1:nfac
	loadings <- scale(fx.svd$v[, fseq, drop=FALSE], 
		scale=1/fx.svd$d[fseq], center=FALSE)

	if(iter.max > 0) {
		cormat <- t(fullxw) %*% fullxw
                uniqueness <- 1 - rowSums(loadings^2)
                uniqueness[uniqueness < 0] <- 0
                uniqueness[uniqueness > 1] <- 1
		start <- uniqueness
		converged <- FALSE
                for(i in 1:iter.max) {
                        cor.red <- cormat
                        diag(cor.red) <- diag(cor.red) - uniqueness
                        t.eig <- eigen(cor.red)
                        t.val <- t.eig$value[fseq]
                        t.val[t.val < 0] <- 0
                        loadings <- scale(t.eig$vector[, fseq, drop=FALSE], 
				center=FALSE, scale=1/sqrt(t.val))
                        uniqueness <- 1 - rowSums(loadings^2)
                        uniqueness[uniqueness < 0] <- 0
                        uniqueness[uniqueness > 1] <- 1
                        if(all(abs(uniqueness - start) < tol)) {
                                converged <- TRUE
                                break
                        }
                        start <- uniqueness
                }
	}
	dimnames(loadings) <- list(dimnames(fullxw)[[2]], NULL)

	if(any(num.mis > 0)) {
		# calculate loadings for columns with NAs
		floadings <- loadings
		if(zero.load) {
			loadings <- array(0, c(nassets, nfac))
		} else {
			meanload <- colMeans(floadings)
			loadings <- t(array(meanload, c(nfac, nassets)))
		}
		dimnames(loadings) <- list(dimnames(x)[[2]], NULL)
		loadings[dimnames(floadings)[[1]], ] <- floadings
		scores <- fullxw %*% floadings
		dsquare <- fx.svd$d[1:nfac]^2
		nfac.miss[nfac.miss > nfac] <- nfac
		
		for(i in (1:nassets)[num.mis > 0 & nobs - num.mis > reg.min]) {
			t.nfac <- nfac.miss[ num.mis[i] ]
			if(t.nfac == 0) next
			t.okay <- !is.na(xw[, i])
			t.seq <- 1:t.nfac
			t.load <- lsfit(xw[t.okay, i], scores[t.okay, t.seq], 
				intercept=FALSE)$coef / dsquare[t.seq]
			loadings[i, t.seq] <- t.load
			NULL
		}
	}

	comm <- rowSums(loadings^2)
	if(any(comm > 1)) {
		# adjust loadings where communalities too large
		toobig <- comm > 1
		if(verbose >= 2) {
			toobigwn <- comm > 1 + 1e-10
			anam <- dimnames(loadings)[[1]]
			if(!length(anam)) {
				anam <- paste("V", 1:nrow(loadings))
			}
			if(sum(toobigwn)) {
			   warning(paste(sum(toobig), "asset(s) being adjusted",
				"from negative specific variance",
				"-- the assets are:", paste(anam[toobig],
				collapse=", "), 
				"(warning suppressed with verbose < 2)"))
			}
		}
		loadings[toobig,] <- loadings[toobig,] / sqrt(comm[toobig])
		comm[toobig] <- 1
	}

	uniqueness <- 1 - comm
	if(!is.numeric(specific.floor) || length(specific.floor) != 1) {
		stop(paste("'specific.floor' should be a single number",
			"-- given has mode", mode(specific.floor),
			"and length", length(specific.floor)))
	}
	if(specific.floor > 0) {
		if(!is.character(floor.type) || length(floor.type) != 1) {
			stop(paste("'floor.type' must be a single character",
				"string -- given has mode", mode(floor.type),
				"and length", length(floor.type)))
		}
		floor.menu <- c("quantile", "fraction")
		floor.num <- pmatch(floor.type, floor.menu, nomatch=0)
		if(floor.num == 0) {
			stop(paste("unknown or ambiguous input for",
				"'floor.type' -- valid choices are:",
				paste(floor.menu, collapse=", ")))
		}
		floor.type <- floor.menu[floor.num]
		switch(floor.type,
			"quantile" = {
				uf <- quantile(uniqueness, specific.floor)
				uniqueness[uniqueness < uf] <- uf
			},
			"fraction" = {
				uniqueness[uniqueness < specific.floor] <- 
					specific.floor
			}
		)
	}
	
	switch(output, 
		full= {
			ans <- loadings %*% t(loadings)
			ans <- t(sdev * ans) * sdev
			diag(ans) <- diag(ans) + uniqueness * sdev^2
			attr(ans, "number.of.factors") <- ncol(loadings)
			attr(ans, "timestamp") <- date()
		},
		systematic=,
		specific=,
		factor={
			ans <- list(loadings=loadings, 
				uniquenesses= uniqueness, sdev=sdev, 
				timestamp=date(), call=match.call())
			class(ans) <- "statfacmodBurSt"
		})
	if(output == "systematic" || output == "specific") {
		fitted(ans, output=output)
	} else {
		ans
	}
}

