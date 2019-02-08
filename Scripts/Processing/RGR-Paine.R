### Appendix 1: R script                   ###
### Paine et al "Best practice for         ###
### growth analysis: model implementations ###
### and growth rate calculations to foster ###
### ecological inference"                  ###
### 3 June 2011                            ###

#### Load libraries and set universal params ####
library(nlme)
library(ggplot2)
library(RColorBrewer)

# Set some global variables for all analyses that follow. 
COL     <- brewer.pal(8, "Dark2")
n.preds <- 100                 # At how many timepoints should predicted biomass and growth rates be calculated?

#### Prep data sets for analysis ####
dat_nonasymp <-
structure(list(X = c(7L, 7L, 7L, 14L, 14L, 14L, 21L, 21L, 21L, 
28L, 28L, 28L, 42L, 42L, 42L, 56L, 56L, 56L, 70L, 70L, 83L, 70L, 
83L, 83L), Y = c(0.000442, 0.000498, 0.000512, 0.002286, 0.002692, 
0.00306, 0.0117, 0.0203, 0.0216, 0.0514, 0.0646, 0.0936, 0.891, 
0.91, 1.173, 4.775, 4.95, 6.082, 11.156, 11.35, 12.539, 13.712, 
15.07, 15.58), logY = c(-7.72420067588658, -7.60491048093962, 
-7.57718593292477, -6.08095171360952, -5.91747086719966, -5.78934036301785, 
-4.44816643717843, -3.8971343929344, -3.83506196429202, -2.96811710652102, 
-2.73954086819358, -2.36872489549859, -0.115410851511328, -0.0943106794712413, 
0.159564569671338, 1.56339397393269, 1.5993875765806, 1.80533358925517, 
2.41197746976524, 2.42921774392741, 2.5288437872084, 2.61827136185542, 
2.7127060126384, 2.7459880404426)), .Names = c("X", "Y", "logY"
), row.names = c(NA, -24L), class = "data.frame")

dat_asymp <-
structure(list(X = c(14L, 14L, 14L, 35L, 35L, 35L, 70L, 133L, 
70L, 161L, 98L, 161L, 133L, 98L, 189L, 189L, 189L, 98L, 133L, 
161L), Y = c(0.031, 0.087, 0.09, 0.261, 0.291, 0.437, 2.104, 
2.736, 2.814, 5.832, 2.979, 9, 6.309, 3.483, 6.103, 5.655, 5.889, 
4.462, 8.392, 7.043)), .Names = c("X", "Y"), class = "data.frame", row.names = c(NA, 
-20L))

dat_asymp_spp <-
structure(list(X = c(14L, 14L, 14L, 70L, 70L, 35L, 35L, 35L, 
189L, 133L, 161L, 70L, 189L, 189L, 133L, 161L, 133L, 98L, 98L, 
98L, 14L, 14L, 14L, 35L, 35L, 35L, 70L, 133L, 70L, 161L, 98L, 
161L, 133L, 98L, 189L, 189L, 189L, 98L, 133L, 161L), Y = c(1.125, 
1.517, 2.667, 2.875, 4.668, 5.862, 5.469, 6.702, 4.331, 9.404, 
11.712, 10.392, 7.602, 8.897, 10.055, 17.204, 14.298, 10.078, 
8.19, 11.075, 0.031, 0.087, 0.09, 0.261, 0.291, 0.437, 2.104, 
2.736, 2.814, 5.832, 2.979, 9, 6.309, 3.483, 6.103, 5.655, 5.889, 
4.462, 8.392, 7.043), species = structure(c(2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L), .Label = c("CER", "GER"), class = "factor")), .Names = c("X", 
"Y", "species"), class = "data.frame", row.names = c(NA, -40L
))


Xes_nonasymp    <- data.frame(X = seq(min(dat_nonasymp$X), max(dat_nonasymp$X), length = n.preds))
Xes_asymp        <- data.frame(X = seq(min(dat_asymp$X),     max(dat_asymp$X),     length = n.preds))
Xes_annuals_spp  <- data.frame(X = seq(min(dat_asymp_spp$X), max(dat_asymp_spp$X), length = n.preds))
dat_asymp_spp   <- groupedData(Y ~ X | species, dat_asymp_spp) # this is needed for Figure 3


#### FUNCTIONS FOR ANALYSIS ####
# Main program code starts at line 1044 below

# Self-start function for linear fit. THis is not necessary, since lm() will do a perfectly adequate job. It's included only for consistency with the other models.
Init.lin <- function(mCall,LHS,data) {
 	xy    <- sortedXyData(mCall[["X"]],LHS,data)
	a     <-coef(lm(y ~ x, xy))[2]  # Use the slope from a linear fit to the data as an initial guess for r
	M0    <- min(xy$y)		        # Use the minimum y value as an initial guess for M0
 	value <- c(M0, a)
 	names(value) <- mCall[c("a", "M0")]
 	return(value)
}
fmla.lin <- as.formula("~M0 + a*X")
SSlin    <- selfStart(fmla.lin, initial = Init.lin, parameters = c("a", "M0"))

# Self-start function for linear fit forced through the origin. THis is not necessary, since lm() will do a perfectly adequate job. It's included only for consistency with the other models.
Init.linno <- function(mCall,LHS,data) {
 	xy    <- sortedXyData(mCall[["X"]],LHS,data)
	a     <-coef(lm(y ~ x, xy))[2]  # Use the slope from a linear fit to the data as an initial guess for r
 	value <- a
 	names(value) <- mCall[c("a")]
 	return(value)
	}
fmla.linno  <- as.formula("~a*X")
SSlinno   <- selfStart(fmla.linno, initial = Init.linno, parameters = c("a"))

# Self-start function for exponential fit. THis is not necessary, since lm(log(y)~ ...) will do a perfectly adequate job. It's included only for consistency with the other models.
Init.exp <- function(mCall,LHS,data) {
 	xy    <- sortedXyData(mCall[["X"]],LHS,data)
	r     <- (coef(lm(log(y) ~ x, xy))[2]) # Use the slope from a linear fit to the data as an initial guess for r
	M0    <- min(xy$y)		             # Use the minimum y value as an initial guess for A
 	value <- c(M0, r)
 	names(value) <- mCall[c("M0", "r")]
 	return(value)
	}
fmla.exp <- as.formula("~M0*exp(r*X)")

SS.exp   <- selfStart(fmla.exp, initial = Init.exp, parameters = c("M0", "r"))

# Self-start function for power-law fit. This self-start for this model is not provided by Pinhero and Bates (2000), so I wrote one. Unfortunately, it often does not lead to convergence. Any suggestions would be welcome
fmla.pow <- as.formula("~(M0^(1-beta) + r*x*(1-beta))^(1/(1-beta))")
Init.pow <- function(mCall, LHS, data){
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    if(nrow(xy) < 4) {stop("Too few distinct x values to fit a power-law")}
	r    <-coef(lm(log(y) ~ x, xy))[2]    # Use the slope from a log fit to the data as an initial guess for r
	M0   <- min(xy$y)		              # Use the minimum y value as an initial guess for M0
	beta <- 0.9	                          # give initial guess of beta as 0.9. don't use beta = 1, as it's undefined (1/0)
    value <- c(M0, r, beta) 
    names(value) <- mCall[c("M0", "r", "beta")]
 	return(value)
	}
SS.pow  <- selfStart(fmla.pow, initial = Init.pow, parameters = c("M0", "r", "beta"))


# Self-start function for 2-parameter power-law fit (initial mass fixed by used). This self-start for this model is not provided by Pinhero and Bates (2000), so I wrote one. 
Init.pow2 <- function(mCall, LHS, data){
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    if(nrow(xy) < 4) {stop("Too few distinct x values to fit a power-law")}
	r    <-coef(lm(log(y) ~ x, xy))[2]    # Use the slope from a log fit to the data as an initial guess for r
	beta <- 0.9	                           # give initial guess of beta as 0.9. don't use beta = 1, as it's undefined (1/0)
    value <- c(r, beta) 
    names(value) <- mCall[c("r", "beta")]
 	return(value)
	}


# Calculate the likelihood of any set of parameter values given the data AND A POWER-LAW model
ll.pow <- function(M0, r, beta, x, y){
	N <- length(x)
	y.pred <- (M0^(1-beta) + r*x*(1-beta))^(1/(1-beta))
	resid <- y - y.pred
	logLik <- -N * (log(2 * pi) + 1 - log(N) + log(sum(resid^2)))/2
	return(logLik)
	}

# Brute-force search for parameter estimates fr the POWER-LAW model (Version 1)
# THis will run many times faster on a Mac or unix machine, because they can take advantage of the multicore package. On windows, only one core ca be used at a time, and it will run much slower. 
grid.search.pow <- function(lim_lo, lim_hi, len, X, Y){
	seq_1 <- seq(lim_lo[1], lim_hi[1], length = len[1])
	seq_2 <- seq(lim_lo[2], lim_hi[2], length = len[2])
	seq_3 <- seq(lim_lo[3], lim_hi[3], length = len[3])
	n_1 <- length(seq_1)
	n_2 <- length(seq_2)
	n_3 <- length(seq_3)
	if(.Platform$OS.type != "windows"){# as of June 2011, these packages do not exist for Windows
		pow.out <- foreach(ii=1:n_1, .combine = rbind) %dopar% {
			pow.out.ii <- data.frame(matrix(NA, ncol = 4, nrow = n_2 * n_3))
			names(pow.out.ii) <- c("M0", "r", "beta", "LL")
			cat(ii, "\t")
			counter <- 0
			for(jj in 1:n_2){
				for(kk in 1:n_3){
					counter <- counter + 1
					pred <- ll.pow(seq_1[ii], seq_2[jj], seq_3[kk], X, Y)
					pow.out.ii[counter, ]  <- c(seq_1[ii], seq_2[jj], seq_3[kk], pred)
				}
			}
		pow.out.ii
		}
	} else {
		pow.out <- data.frame(matrix(NA, ncol = 4, nrow = n_1*n_2*n_3))
		names(pow.out) <- c("M0", "r", "beta", "LL")
		counter <- 0
		for(ii in 1:n_1){
			cat(ii, "\t")
			for(jj in 1:n_2){
				for(kk in 1:n_3){
					counter <- counter + 1
					pred <- ll.pow(seq_1[ii], seq_2[jj], seq_3[kk], X, Y)
					pow.out[counter, ]  <- c(seq_1[ii], seq_2[jj], seq_3[kk], pred)
					}
				}
			}		
		
		}	
	return(pow.out)
	}



# The following three functions transform the parameters of the logistic, gompertz and monomolecular models (respectively), to put them onto a scale most easily comprable with the other models.
transform_param.logis <- function(coef){
	K = coef[1]
	r = 1/(coef[3])
	M0 =  K/(1 + exp(coef[2]/coef[3])) #untransform best-fit parameters to K, r and M0
	if(is.data.frame(K)){
		out <- cbind(K, r, M0)
		} else {
		out <- c(K, r, M0)
		}
	names(out) <- c("K", "r", "M0")
	return(out)
	}

transform_param.gomp <- function(coef){
	K  <- coef[1]
	M0 <- K/exp(coef[2])
	r  <- -log(coef[3])
	out <- c(K, r, M0)
	names(out) <- c("K", "r", "M0")
	return(out)
	}


transform_param.mono <- function(coef){
	K  <- coef[1]
	M0 <- coef[2] 
	r  <- exp(coef[3])
	out <- c(K, r, M0)
	names(out) <- c("K", "r", "M0")
	return(out)
	}

# this function returns confidence envelopes around growth trajectories, and growth rates. 
summarizer <- function(dat, alpha){
	n <- length(dat)
	quantiles <- c(alpha/2, 1-(alpha/2))
	CIs <- data.frame(matrix(NA, ncol(dat[[1]]), n*2))
	names(CIs) <- paste(rep(names(dat), each = 2), c("lo", "hi"), sep = ".")
	for(i in 1:n){
		CIs[,(2*i-1):(2*i)] <- t(apply(dat[[i]],    2, quantile, quantiles, na.rm = T))
		}
	return(CIs)
	}


# The following seven functions return predictions of biomass and growth rates.
output.lin.nls <- function(fit, times, CI = F, alpha = 0.05){
	params <- coef(fit)
	names(params) <- NULL
	a <- params[1]
	if(length(params) > 1){
		M0 <- params[2]
		} else {
		M0 <- 0
		}
	fitted <- fit$m$fitted()
	resid  <- fit$m$resid()
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(M0, 2)) + .(round(a, 3)) *t))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    =  M0 + a*times,
		AGR  = a)
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm  <-  a/rates$M
	if(CI == T){
		cov  <- summary(fit)$cov
		x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
		if(ncol(x) == 1){x$M0 = 0}
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <- x.i$M0 + x.i$a*times
			AGR[i,]   <- rep(x.i$a, length(times))
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$a/M[i,]
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.lin.gnls <- function(fit, times, CI = F, alpha = 0.05){
	params <- coef(fit)
	names(params) <- NULL
	a <- params[1]
	if(length(params) > 1){
		M0 <- params[2]
		} else {
		M0 <- 0
		}
	fitted <- fit$fitted
	resid  <- fit$residual
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(M0, 2)) + .(round(a, 3)) *t))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    =  M0 + a*times,
		AGR  = a)
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm  <-  a/rates$M
	if(CI == T){
		cov  <- fit$varBeta
		x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
		if(ncol(x) == 1){x$M0 = 0}
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <- x.i$M0 + x.i$a*times
			AGR[i,]   <- rep(x.i$a, length(times))
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$a/M[i,]
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.exp.nls <- function(fit, times, CI = F, alpha = 0.05){
	params <- coef(fit)
	names(params) <- NULL
	M0 <- params[1]; r <- params[2]
	fitted <- fit$m$fitted()
	resid  <- fit$m$resid()
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(.(round(M0, 3)) * e^{.(round(r, 3))*t})
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    =   M0*exp(r*times),
		AGR  = r*M0*exp(r*times),
		RGRt = r,
		RGRm = r
		)
	if(CI == T){
		cov  <- summary(fit)$cov
		x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <-       x.i$M0*exp(x.i$r*times)
			AGR[i,]   <- x.i$r*x.i$M0*exp(x.i$r*times)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$r
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		names(params) <- c("M0", "r")
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.exp.gnls <- function(fit, times, CI = F, alpha = 0.05){
	params <- coef(fit)
	names(params) <- NULL
	M0 <- params[1]; r <- params[2]
	fitted <- fit$fitted
	resid  <- fit$resid
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(.(round(M0, 3)) * e^{.(round(r, 3))*t})
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    =   M0*exp(r*times),
		AGR  = r*M0*exp(r*times),
		RGRt = r,
		RGRm = r
		)
	if(CI == T){
		cov  <- fit$varBeta
		x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <-       x.i$M0*exp(x.i$r*times)
			AGR[i,]   <- x.i$r*x.i$M0*exp(x.i$r*times)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$r
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		names(params) <- c("M0", "r")
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.pow <- function(X, Y, params, times){
	M0 <- params$M0; r <- params$r; beta <- params$beta
	fitted <- (M0^(1-beta)   + r  * X *(1-beta))^(1/(1-beta))
	resid  <- Y - fitted
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
	mss  <- sum((fitted - mean(fitted))^2)
	rss  <- sum(resid^2)
	R2   <- mss/(mss + rss)
	rmse <- sqrt(rss)
	N <- length(X)
	logLik <- -N * (log(2 * pi) + 1 - log(N) + log(sum(resid^2)))/2
	AIC  <- -2 * logLik  + 2 * 3 # three parameters
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	temp <- M0^(1-beta) + r*times*(1-beta)
	rates = data.frame(
		times = times,
		M    =     temp^(1/(1-beta)),
		AGR  = r * temp^(beta/(1-beta)))
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm  <-  r * rates$M^(beta - 1)
	out <- list(params = params[-4], summary = summary, equation = eq, data = data, rates = rates)
	return(out)
	}

output.pow2.nls <- function(fit, times, CI = F, alpha = 0.05, M0){
	params <- c(coef(fit), M0 = M0)
	r = params[1]; beta = params[2]
	fitted <- fit$m$fitted()
	resid  <- fit$m$resid()
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
	mss  <- sum((fitted - mean(fitted))^2)
	rss  <- sum(resid^2)
	R2   <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	temp <- M0^(1-beta) + r*times*(1-beta)
	rates = data.frame(
		times = times,
		M    =     temp^(1/(1-beta)),
		AGR  = r * temp^(beta/(1-beta)))
	rates$RGRt  <- rates$AGR/rates$M
	rates$RGRm  <- r * rates$M^(beta-1)
	if(CI == T){
		cov <- summary(fit)$cov
		cov <- cbind(rbind(cov, 0), 0) # variance-covariance for M0 is 0 - it's fixed. 
		x <- data.frame(rmvnorm(n=1000, mean=params, sigma=cov))
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			temp.i    <- x.i$M0^(1-x.i$beta) + x.i$r * times *(1-x.i$beta)
			M[i,]     <-         temp.i^(1/(1-x.i$beta))
			AGR[i,]   <- x.i$r * temp.i^(x.i$beta/(1-x.i$beta))
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$r * M[i,]^(x.i$beta - 1)
			}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.pow2.gnls <- function(fit, times, CI = F, alpha = 0.05, M0){
	params <- c(coef(fit), M0 = M0)
	r = params[1]; beta = params[2]
	fitted <- fit$fitted
	resid  <- fit$residual
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
	mss  <- sum((fitted - mean(fitted))^2)
	rss  <- sum(resid^2)
	R2   <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	temp <- M0^(1-beta) + r*times*(1-beta)
	rates = data.frame(
		times = times,
		M    =     temp^(1/(1-beta)),
		AGR  = r * temp^(beta/(1-beta)))
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm  <-  r * rates$M^(beta-1)
	if(CI == T){
		cov <- fit$varBeta
		cov <- cbind(rbind(cov, 0), 0) # variance-covariance for M0 is 0 - it's fixed. 
		x <- data.frame(rmvnorm(n=1000, mean=params, sigma=cov))
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			temp.i   <- x.i$M0^(1-x.i$beta) + x.i$r * times *(1-x.i$beta)
			M[i,]    <-         temp.i^(1/(1-x.i$beta))
			AGR[i,]  <- x.i$r * temp.i^(x.i$beta/(1-x.i$beta))
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,] <- x.i$r * M[i,]^(x.i$beta - 1)
			}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.mono.nls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	coef <- coef(fit)
	params <- transform_param.mono(coef)
	K = params[1]; r = params[2]; M0 = params[3]
	fitted <- fit$m$fitted()
	resid  <- fit$m$resid()
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(r, 4)) - .(round(K - M0, 2)) * e^{.(round(r, 4))*t}))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M     = K-(exp(-r*times)*(K-M0)),
		AGR   = r*exp(-r*times)*(K-M0)) #dM/dt
		rates$RGRt <- rates$AGR/rates$M
		rates$RGRm <- (r*(K-rates$M))/rates$M
	if(LOG ==T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- (r*(K-rates$M))
		rates$AGR  <- rates$AGR*exp(rates$M)
		}
	if(CI == T){
		cov <- summary(fit)$cov
		y <- x <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
		x$K  <- y$Asym
		x$r  <- exp(y$lrc)
		x$M0 <- y$R0 #untransform best-fit parameters to K, r and M0
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <-  x.i$K-exp(-x.i$r*times)*(x.i$K-x.i$M0)
			AGR[i,]   <-  x.i$r*exp(-x.i$r*times)*(x.i$K-x.i$M0)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$r*(x.i$K - M[i,])/M[i,]
			if(LOG ==T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- (x.i$r*(x.i$K-M[i,]))
				AGR[i,]  <- AGR[i,]*exp(M[i,])
				}
			}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}
	
output.mono.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	coef <- coef(fit)
	params <- transform_param.mono(coef)
	K = params[1]; r = params[2]; M0 = params[3]
	fitted <- fit$fitted
	resid  <- fit$residual
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(r, 4)) - .(round(K - M0, 2)) * e^{.(round(r, 4))*t}))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M     = K-(exp(-r*times)*(K-M0)),
		AGR   = r*exp(-r*times)*(K-M0)) #dM/dt
		rates$RGRt <- rates$AGR/rates$M
		rates$RGRm <- (r*(K-rates$M))/rates$M
	if(LOG ==T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- r*(K-rates$M)
		rates$AGR  <- rates$AGR*exp(rates$M)
		}
	if(CI == T){
		cov <- fit$varBeta
		y <- x <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
		x$K  <- y$Asym
		x$r  <- exp(y$lrc)
		x$M0 <- y$R0 #untransform best-fit parameters to K, r and M0
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]    <-  x.i$K-exp(-x.i$r*times)*(x.i$K-x.i$M0)
			AGR[i,]  <-  x.i$r*exp(-x.i$r*times)*(x.i$K-x.i$M0)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,] <- x.i$r*(x.i$K - M[i,])/M[i,]
			if(LOG ==T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- (x.i$r*(x.i$K-M[i,]))
				AGR[i,]  <- AGR[i,]*exp(M[i,])
				}
			}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.logis.nls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	coef <- coef(fit)
	params <- transform_param.logis(coef)
	K = params[1]; r = params[2]; M0 = params[3]
	fitted <- fit$m$fitted()
	resid  <- fit$m$resid()
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(M0*r, 4)) /(.(round(M0, 3)) + .(round(K-M0, 2)) * e^{.(round(-r, 3))*t})))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    = (M0*K)/(M0+(K-M0)*exp(-r*times)),
		AGR  = (r*M0*K*(K-M0)*exp(-r*times))/(M0+(K-M0)*exp(-r*times))^2)
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm <- r*(1 - rates$M/K)
	if(LOG ==T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- r*rates$M*(1-rates$M/K)
		rates$AGR  <- rates$AGR*exp(rates$M)	
		}
	if(CI == T){
		cov   <- summary(fit)$cov
		x <- y <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
		x$K  <- y$Asym
		x$r  <- 1/y$xmid
		x$M0 <- y$Asym/(1 + exp(y$xmid/y$scal)) #untransform best-fit parameters to K, r and M0
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <- (x.i$M0*x.i$K)/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))
			AGR[i,]   <- (x.i$r*x.i$M0*x.i$K*(x.i$K-x.i$M0)*exp(-x.i$r*times))/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))^2
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <-  x.i$r*(1 - M[i,]/x.i$K)
			if(LOG ==T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- x.i$r*M[i,]*(1 - M[i,]/x.i$K)
				AGR[i,]  <- AGR[i,]*exp(M[i,])
				}
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}
	
output.logis.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	coef <- coef(fit)
	params <- transform_param.logis(coef)
	K = params[1]; r = params[2]; M0 = params[3]
	fitted <- fit$fitted
	resid  <- fit$residuals
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(M0*r, 4)) /(.(round(M0, 3)) + .(round(K-M0, 2)) * e^{.(round(-r, 3))*t})))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    = (M0*K)/(M0+(K-M0)*exp(-r*times)),
		AGR  = (r*M0*K*(K-M0)*exp(-r*times))/(M0+(K-M0)*exp(-r*times))^2)
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm <- r*(1 - rates$M/K)
	if(LOG ==T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- r*rates$M*(1-rates$M/K)
		rates$AGR  <- rates$AGR*exp(rates$M)
		}
	if(CI == T){
		cov   <- fit$varBeta
		x <- y <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
		x$K  <- y$Asym
		x$r  <- 1/y$xmid
		x$M0 <- y$Asym/(1 + exp(y$xmid/y$scal)) #untransform best-fit parameters to K, r and M0
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]     <- (x.i$M0*x.i$K)/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))
			AGR[i,]   <- (x.i$r*x.i$M0*x.i$K*(x.i$K-x.i$M0)*exp(-x.i$r*times))/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))^2
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <-  x.i$r*(1 - M[i,]/x.i$K)
			if(LOG ==T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- x.i$r*M[i,]*(1 - M[i,]/x.i$K)
				AGR[i,]  <- AGR[i,]*exp(M[i,])
				}
			}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}
	
# We fit the logistic function unsing nlsList to compare fits between two species. This comparison is done in nlsList, which has a somewhat different structure than nls. So, we need a new function to extract predictions
output.logis.nlsList <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	coef <- coef(fit)
	params <- transform_param.logis(coef)
	rates <- list()
	groups <- rownames(params)
	n.groups <- nrow(coef)

	# compute rates for each group seperately
	rates <- list()
	for(i in 1:(n.groups)){
		K <- params[i,1]; r <- params[i,2]; M0 <- params[i,3]
		rates[[i]] = data.frame(
			times = times,
			M    = (M0*K)/(M0+(K-M0)*exp(-r*times)),
			AGR  = (r*M0*K*(K-M0)*exp(-r*times))/(M0+(K-M0)*exp(-r*times))^2
			)
		rates[[i]]$RGRt <- rates[[i]]$AGR/rates[[i]]$M
		rates[[i]]$RGRm <- r*(1 - rates[[i]]$M/K)
		if(LOG == T){
			rates[[i]]$RGRt <- rates[[i]]$AGR
			rates[[i]]$RGRm <- r*rates[[i]]$M*(1-rates[[i]]$M/K)
			rates[[i]]$AGR  <- rates[[i]]$AGR*exp(rates[[i]]$M)
			}
		# commute CIs for each group's estaimates, if desired
		if(CI == T){
			cov   <- summary(fit)$cov[i,,]
			x <- y <- data.frame(rmvnorm(n=1000, mean=c(coef[i, 1], coef[i, 2], coef[i, 3]), sigma=cov))
			x$K  <- y[,1]
			x$r  <- 1/y[,3]
			x$M0 <- y[,1]/(1 + exp(y[,2]/y[,3])) 
			M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
			for(j in 1:nrow(x)){
				K <- x[j,4]; r <- x[j,5]; M0 <- x[j,6]
				M[j,]     <- (M0*K)/(M0+(K-M0)*exp(-r*times))
				AGR[j,]   <- (r*M0*K*(K-M0)*exp(-r*times))/(M0+(K-M0)*exp(-r*times))^2
				RGRt[j,]  <- AGR[j,]/M[j,]
				RGRm[j,]  <- r*(1 - M[j,]/K)
				if(LOG ==T){
					RGRt[j,] <- AGR[j,]
					RGRm[j,] <- r*M[j,]*(1 - M[j,]/K)
					AGR[j,]  <- AGR[j,]*exp(M[j,])
					}
				}
			}
			CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
			rates[[i]]  <-  cbind(rates[[i]], CIs)
		}
	names(rates) <- rownames(params)

	# now compute differences among groups
	diffs <- list()
	for(i in 1:(n.groups-1)){
		Ki <- params[i,1]; ri <- params[i,2]; M0i <- params[i,3]
		for(j in (i+1):n.groups){
			Kj <- params[j,1]; rj <- params[j,2]; M0j <- params[j,3]
			diffs.ij = data.frame(
				times = times,
				diffM     = rates[[i]]$M    - rates[[j]]$M,
				diffAGR   = rates[[i]]$AGR  - rates[[j]]$AGR,
				diffRGRt  = rates[[i]]$RGRt - rates[[j]]$RGRt
				)
				# comparing RGRm has to be done on a biomass basis. So it needs special treatment. First, we aev to know what range of biomasses are shared between groups.
				Mmin <- max(min(rates[[i]]$M), min(rates[[j]]$M)) # yieds the range of overlapping masses between groups i & j
				Mmax <- min(max(rates[[i]]$M), max(rates[[j]]$M))
				diffs.ij$Mseq <- seq(Mmin, Mmax, length = 100)
				if(LOG == F){
					diffs.ij$diffRGRm  <- ri*(1 - diffs.ij$Mseq/Ki) - rj*(1 - diffs.ij$Mseq/Kj)
					} else{
					diffs.ij$diffRGRm <- ri*diffs.ij$Mseq*(1 - diffs.ij$Mseq/Ki) - rj*Mseq*(1 - diffs.ij$Mseq/Kj)
				}		
				
			}
		if(CI == T){
			# get params for group i
			covi   <- summary(fit)$cov[i,,]
			xi <- yi <- data.frame(rmvnorm(n=1000, mean=c(coef[i, 1], coef[i, 2], coef[i, 3]), sigma=covi))
			xi$K  <- yi[,1]
			xi$r  <- 1/yi[,3]
			xi$M0 <- yi[,1]/(1 + exp(yi[,2]/yi[,3])) 

			# get params for group j
			covj   <- summary(fit)$cov[j,,]
			xj <- yj <- data.frame(rmvnorm(n=1000, mean=c(coef[j, 1], coef[j, 2], coef[j, 3]), sigma=covj))
			xj$K  <- yj[,1]
			xj$r  <- 1/yj[,3]
			xj$M0 <- yj[,1]/(1 + exp(yj[,2]/yj[,3])) 
			
			# now compute diffs for each random set of drawn parameters
			Mi <- Mj <- AGRi <- AGRj <- RGRti <- RGRtj <- RGRmi <- RGRmj <- diffM <- diffAGR <- diffRGRt <- diffRGRm <- matrix(NA, ncol = length(times), nrow = nrow(xi))
			for(k in 1:nrow(xi)){
				Ki <- xi[k,4]; ri <- xi[k,5]; M0i <- xi[k,6]
				Kj <- xj[k,4]; rj <- xj[k,5]; M0j <- xj[k,6]
				Mi[k,]    <- (M0i*Ki)/(M0i+(Ki-M0i)*exp(-ri*times))
				Mj[k,]	  <- (M0j*Kj)/(M0j+(Kj-M0j)*exp(-rj*times))
				AGRi[k,]  <- (ri*M0i*Ki*(Ki-M0i)*exp(-ri*times))/(M0i+(Ki-M0i)*exp(-ri*times))^2
				AGRj[k,]  <- (rj*M0j*Kj*(Kj-M0j)*exp(-rj*times))/(M0j+(Kj-M0j)*exp(-rj*times))^2
				RGRti[k,] <- AGRi[k,]/Mi[k,]
				RGRtj[k,] <- AGRj[k,]/Mj[k,]
				RGRmi[k,] <- ri*(1 - diffs.ij$Mseq/Ki)
				RGRmj[k,] <- rj*(1 - diffs.ij$Mseq/Kj)
				if(LOG == T){
					RGRti[k,] <- AGRi[k,]
					RGRtj[k,] <- AGRj[k,]
					RGRmi[k,] <- ri*diffs.ij$Mseq*(1 - diffs.ij$Mseq/Ki)
					RGRmj[k,] <- rj*diffs.ij$Mseq*(1 - diffs.ij$Mseq/Kj)
					AGRi[k,]  <- AGRi[k,]*exp(Mi[k,])
					AGRj[k,]  <- AGRj[k,]*exp(Mj[k,])
					}
				diffM[k,]    <- Mi[k,]    - Mj[k,]
				diffAGR[k,]  <- AGRi[k,]  - AGRj[k,]
				diffRGRt[k,] <- RGRti[k,] - RGRtj[k,]
				diffRGRm[k,] <- RGRmi[k,] - RGRmj[k,]
				}
			CIs <- summarizer(list(diffM = diffM, diffAGR = diffAGR, diffRGRt = diffRGRt, diffRGRm = diffRGRm), alpha)
			diffs[[paste(groups[i], groups[j], sep = "_")]]  <-  cbind(diffs.ij, CIs)
			} else{
			diffs[[paste(groups[i], groups[j], sep = "_")]]  <-  diffs.ij
			}
		} # end loop over pairwise combinations of groups
	
	out <- list(params = params, rates = rates, diffs = diffs)
	return(out)
	}

output.fpl.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	params <- coef <- coef(fit)
	M0 = params[1]; K = params[2]; xmid = params[3]; r = params[4]
	fitted <- fit$fitted
	resid  <- fit$residuals
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(M0, 4))+ (.(round(K-M0, 4)) /1+e^{(.(round(xmid, 3))*t)/.(round(r, 3))})))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	exp <- exp((xmid-times)/r)
	rates = data.frame(
		times = times,
		M    = M0+(K-M0)/(1+exp),
		AGR  = ((K-M0)*exp)/(r*(1+exp)^2)
		)
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm <- (M0-rates$M)*(K-rates$M)/(r*rates$M*(M0-K))
	if(LOG == T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- (M0-rates$M)*(K-rates$M)/(r*(M0-K))
		rates$AGR  <- rates$AGR*exp(rates$M)
		}
	if(CI == T){
		cov   <- fit$varBeta
		x <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
		names(x) <- c("M0", "K", "xmid", "r")
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			exp <- exp((x.i$xmid-times)/x.i$r)
			M[i,]     <- x.i$M0 + (x.i$K-x.i$M0)/(1+exp)
			AGR[i,]   <- ((x.i$K-x.i$M0)*exp)/(x.i$r*(1+exp)^2)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,]  <- x.i$r*(1 - (M[i,]/x.i$K)^(1/x.i$xmid))
#			RGRm[i,]  <- ((x.i$M0-M[i,])*(x.i$K-M[i,]))/(M[i,]*x.i$r*(x.i$M0-x.i$K))
			if(LOG == T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- (x.i$M0-M[i,])*(x.i$K-M[i,])/(x.i$r*(x.i$M0-x.i$K))
				AGR[i,]  <- AGR[i,]*exp(M[i,]) 
				}
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}
		
output.gomp.nls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	coef <- coef(fit)
	params <- transform_param.gomp(coef)
	K = params[1]; r = params[2]; M0 = params[3]
	fitted <- fit$m$fitted()
	resid  <- fit$m$resid()
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(K, 2)) %*% .(round(M0 / K, 5)) ^ e^{.(round(r, 3))*t}))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    = K*((M0/K)^exp(-r*times)),
		AGR  = r*K*exp(-r*times)*log(K/M0)*(M0/K)^exp(-r*times))
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm <- r*log(K/M)
	if(LOG == T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- r*rates$M*log(K/rates$M)
		rates$AGR  <- rates$AGR*exp(rates$M)
		}
	if(CI == T){
		cov  <- summary(fit)$cov
		x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
		x$K  <- x$Asym
		x$M0 <- x$K/exp(x$b2)
		x$r  <- -log(x$b3)
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]    <- x.i$K*((x.i$M0/x.i$K)^exp(-x.i$r*times))
			AGR[i,]  <- x.i$r*x.i$K*exp(-x.i$r*times)*log(x.i$K/x.i$M0)*(x.i$M0/x.i$K)^exp(-x.i$r*times)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,] <- x.i$r*log(x.i$K/M[i,])
			if(LOG ==T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- x.i$r*M[i,]*log(x.i$K/M[i,])
				AGR[i,]  <- AGR[i,]*exp(M[i,])
				}
			}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}

output.gomp.gnls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
	coef <- coef(fit)
	params <- transform_param.gomp(coef)
	K = params[1]; r = params[2]; M0 = params[3]
	fitted <- fit$fitted
	resid  <- fit$residual
	data <- data.frame(fitted = fitted, resid = resid)
	eq   <- bquote(paste(.(round(K, 2)) %*% .(round(M0 / K, 5)) ^ e^{.(round(r, 3))*t}))
	mss <- sum((fitted - mean(fitted))^2)
	rss <- sum(resid^2)
	R2  <- mss/(mss + rss)
	rmse <- sqrt(rss)
	AIC <- AIC(fit)
	summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
	rates = data.frame(
		times = times,
		M    = K*((M0/K)^exp(-r*times)),
		AGR  = r*K*exp(-r*times)*log(K/M0)*(M0/K)^exp(-r*times))
	rates$RGRt <- rates$AGR/rates$M
	rates$RGRm <- r*log(K/rates$M)
	if(LOG == T){
		rates$RGRt <- rates$AGR
		rates$RGRm <- r*rates$M*log(K/rates$M)
		rates$AGR  <- rates$AGR*exp(rates$M)
		}
	if(CI == T){
		cov  <- fit$varBeta
		x    <- data.frame(rmvnorm(n=1000, mean=coef(fit), sigma=cov))
		x$K  <- x$Asym
		x$M0 <- x$K/exp(x$b2)
		x$r  <- -log(x$b3)
		M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
		for(i in 1:nrow(x)){
			x.i  <- x[i,]
			M[i,]    <- x.i$K*((x.i$M0/x.i$K)^exp(-x.i$r*times))
			AGR[i,]  <- x.i$r*x.i$K*exp(-x.i$r*times)*log(x.i$K/x.i$M0)*(x.i$M0/x.i$K)^exp(-x.i$r*times)
			RGRt[i,] <- AGR[i,]/M[i,]
			RGRm[i,] <- x.i$r*log(x.i$K/M[i,])
			if(LOG ==T){
				RGRt[i,] <- AGR[i,]
				RGRm[i,] <- x.i$r*M[i,]*log(x.i$K/M[i,])
				AGR[i,]  <- AGR[i,]*exp(M[i,])
				}
		}
		CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
		} else {
		out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
		}
	return(out)
	}


#### MAIN PROGRAM CODE ####

# Which fitting routine to use? nls(), gnls(), nlsList(), or nlme()?
# If you have no treatment effects, ie, a model of the form y ~ x, use gnls. It can do everything that nls can do, but additionally allows variance modelling to deal with heteroscedasticity. Also, gnls allows you to choose which parameters should vary among treatment groups, and which should be global. See Pinhero & Bates, pg 401, for an example. However, occasionally gnls will fail to converge. In those cases, nls can still be useful.
# Use nlsList to fit a series of nls models (this is rarely useful in practice. Rather, you probably want to use nlme)
# nlme can do all of the above, as well as allowing for fixed and (nested) random effects. However, it requires a somewhat different syntax for specifying the model.
# See Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", Springer for much more detail on these topics.
###

#### ASYMPTOTIC FITS ####
## Monomolecular fit ####
# Fit to Cerastium data
tmp.mono <- getInitial(Y ~ SSasymp(X, Asym, R0, lrc), data = dat_asymp)
fit.mono <- gnls(Y ~ SSasymp(X, Asym, R0, lrc), data = dat_asymp)
out.mono <- output.mono.gnls(fit.mono, Xes_asymp$X, CI = T)

# gnls fit with variance function simply does not converge
#fit.mono <- gnls(Y ~ SSasymp(X, Asym, R0, lrc), data = dat_asymp, weights= varExp(form = ~ fitted(.))), start = c(Asym = 200, R0 = -1, lrc = -8), control = list(nlsTol = 1, returnObject = T))
#out.mono <- output.mono.gnls(fit.mono, Xes_asymp$X, CI = T)

# Fit to log-Holcus data
tmp.mono     <- getInitial(logY ~ SSasymp(X, Asym, R0, lrc), data = dat_nonasymp)
fit.mono.log <- gnls(logY ~ SSasymp(X, Asym, R0, lrc), data = dat_nonasymp)
out.mono.log <- output.mono.gnls(fit.mono.log, Xes_nonasymp$X, CI = T, LOG = T)

## Logistic fit ####
# Fit to Cerastium data
#tmp.logis  <- getInitial(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp)
#fit.logis  <- nls( Y ~ SSlogis(X, Asym, xmid, scal), trace = F, control = list(maxiter=500), data = dat_asymp)
#summary(fit.logis)
#out.logis <- output.logis.nls(fit.logis, Xes_asymp$X, T)

tmp.logis  <- getInitial(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp)
fit.logis <- gnls(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp, weights= varExp(form = ~ fitted(.)))
out.logis <- output.logis.gnls(fit.logis, Xes_asymp$X, CI = T)

# Fit to log-Holcus data
tmp.logis.log <- getInitial(logY ~ SSlogis(X, Asym, xmid, scal), data = dat_nonasymp)
fit.logis.log <- gnls(logY ~ SSlogis(X, Asym, xmid, scal), data = dat_nonasymp, weights= varExp(form = ~ fitted(.)))
out.logis.log <- output.logis.gnls(fit.logis.log, Xes_nonasymp$X, CI = T, LOG = T)

## 4-parameter logistic fit ####
# Fit to Cerastium data
#tmp.fpl  <- getInitial(Y ~ SSfpl(X, Asym, xmid, scal), data = dat_asymp)
#fit.fpl  <- nls( Y ~ SSfpl(X, Asym, xmid, scal), trace = F, control = list(maxiter=500), data = dat_asymp)
#summary(fit.fpl)
#out.fpl <- output.fpl.nls(fit.fpl, Xes_asymp$X, T)

tmp.fpl  <- getInitial(Y ~ SSfpl(X, A, B, xmid, scal), data = dat_asymp)
fit.fpl <- gnls(Y ~ SSfpl(X, A, B, xmid, scal), data = dat_asymp, weights= varExp(form = ~ fitted(.)))
out.fpl <- output.fpl.gnls(fit.fpl, Xes_asymp$X, CI = T)

# Fit to log-Holcus data
tmp.fpl.log <- getInitial(logY ~ SSfpl(X, A, B, xmid, scal), data = dat_nonasymp)
fit.fpl.log <- gnls(logY ~ SSfpl(X, A, B, xmid, scal), data = dat_nonasymp, weights= varExp(form = ~ fitted(.)))
out.fpl.log <- output.fpl.gnls(fit.fpl.log, Xes_nonasymp$X, CI = T, LOG = T)

## Gompertz ####
# Fit to Cerastium data
#tmp.gomp <- getInitial(Y ~ SSgompertz(X, Asym, b2, b3), data = dat_asymp)
#fit.gomp <- nls(Y ~ SSgompertz(X, Asym, b2, b3), trace = F, control = nls.control(maxiter=500), algorithm = "port", lower = c(Asym = 0, b2 = 0, b3 = 0), data = dat_asymp)
#out.gomp <- output.gomp.nls(fit.gomp, Xes_asymp$X, T)
#summary(fit.gomp)

tmp.gomp <- getInitial(Y ~ SSgompertz(X, Asym, b2, b3), data = dat_asymp)
fit.gomp <- gnls(Y ~ SSgompertz(X, Asym, b2, b3), data = dat_asymp, weights= varExp(form = ~ fitted(.)))
out.gomp <- output.gomp.gnls(fit.gomp, Xes_asymp$X, CI = T)

# Fit to log-Holcus data
# Biomass data must be converted to milligrams from grams prior to log-transformation in order to fit - this is a trick suggested by Mark Rees. It seems to make the gompertz model more likely to converge. Just remember to back convert for plotting & inferences!
tmp.gomp.log <- getInitial(log(1000*Y) ~ SSgompertz(X, Asym, b2, b3), data = dat_nonasymp) 
fit.gomp.log <- gnls(log(1000*Y) ~ SSgompertz(X, Asym, b2, b3), data = dat_nonasymp, weights= varExp(form = ~ fitted(.)))
out.gomp.log <- output.gomp.gnls(fit.gomp.log, Xes_nonasymp$X, CI = T, LOG = T)

#### NON-ASYMPTOTIC FITS ####

## Linear fit with intercept ####
tmp.lin  <- getInitial(Y~SSlin(X, a, M0), dat_nonasymp)
fit.lin  <- gnls(Y ~ SSlin(X, a, M0), data = dat_nonasymp) # predicts negative biomass, which is impossible.
out.lin <- output.lin.gnls(fit.lin, Xes_nonasymp$X, T)

## Linear fit NO intercept ####
tmp.linno <- getInitial(Y~ SSlinno(X, a), dat_nonasymp)
fit.linno <-  gnls(Y ~ SSlinno(X, a), data = dat_nonasymp)
out.linno <- output.lin.gnls(fit.linno, Xes_nonasymp$X, T)

## Exponential fit ####
tmp.exp  <- getInitial(Y~SS.exp(X,r,M0), dat_nonasymp)
fit.exp  <-  nls(Y ~ SS.exp(X, r, M0), data = dat_nonasymp) 
out.exp  <- output.exp.nls(fit.exp, Xes_nonasymp$X, T)

## 3-parameter Power-law fit, including initial mass ####
tmp.pow <- getInitial(Y ~ SS.pow(X, M0, r, beta), data = dat_nonasymp)
## These don't converge, even with initial parameters supplied.... 
#fit.pow <- gnls(Y ~ SS.pow(X, M0, r, beta), data = dat_nonasymp, list(M0 = 0.0043650, r = 0.0421200, beta = 0.8))
#fit.pow <- gnls(Y ~ SS.pow(X, M0, r, beta), start = list(M0 = 0.08, r = 0.04, beta = 0.9), control = nls.control(maxiter = 500), data = dat_nonasymp) 

## Alternative approach to fitting 3-parameter Power-law fit, including initial mass ####
# this routine determines parameter values by minimizing the RSS, while constraining all parameters to be positive

fit.power3 <- function(p, X, Y){
	M0   <- exp(p[1])
	beta <- exp(p[2])
	r    <- exp(p[3])
	y.pred <- (M0^(1-beta) + r*X*(1-beta))^(1/(1-beta))
	RSS    <- sum((Y-y.pred)*(Y-y.pred))/length(X)
	return(RSS)
	}
# Provide the routine with some starting guesses for parameter values
p <- log(c(M0 = 0.08, beta = 0.9, r = 0.04))

tmp <- optim(p, fit.power3, method="Nelder-Mead", X=dat_nonasymp$X, Y=dat_nonasymp$Y, control=list(maxit=20000,trace=0))
# Output parameter values
best.pow.optim <- exp(tmp$par)

## 2-parameter Power-law fit, initial mass fixed ####
seed.mass <- 0.00045 #(45 mg, Hautier et al 2010, J. Ecology) # change this to suit!
fmla.pow2 <- as.formula(paste("~(", seed.mass, "^(1-beta) + r*x*(1-beta))^(1/(1-beta))", sep = ""))
SS.pow2   <- selfStart(fmla.pow2, initial = Init.pow2, parameters = c("r", "beta"))
tmp.pow2 <- getInitial(Y ~ SS.pow2(X, r, beta), data = dat_nonasymp)
fit.pow2 <- gnls(Y ~ SS.pow2(X, r, beta), data = dat_nonasymp)
out.pow2 <- output.pow2.gnls(fit.pow2, Xes_nonasymp$X, CI = T, M0 = seed.mass)

## brute-force search to fit the 3-parameter Power-law fit, including initial mass ####
# STEP 1: Coarse search
##               M0        r     beta
length    <- c(   100,    100,   100)
limit_low <- c(   0.00,    0.00,  0.00)
limit_hi  <- c(   1.00,   0.20,  0.99)
grid.pow <- grid.search.pow(limit_low, limit_hi, length, dat_nonasymp$X, dat_nonasymp$Y) # This will run many times faster on a multi-core mac/unix than on a Windows machine.
write.csv(grid.pow, "coarse brute-force.csv", row.names = F)
grid.pow.coarse <- read.csv("coarse brute-force.csv")
best.pow.coarse <- grid.pow.coarse[which.max(grid.pow.coarse$LL),] #global maximum likelhood
best.pow.coarse

# STEP 2: Fine search
# Use coarse estimate from step one to refine parameter estimates
##               M0     r   beta
length    <- c( 500,  100,   100)
limit_low <- c(0.00,  0.00,  0.20)
limit_hi  <- c(0.002, 0.20,  0.80)
grid.pow <- grid.search.pow(limit_low, limit_hi, length, dat_nonasymp$X, dat_nonasymp$Y) # This will run many times faster on a multi-core mac/unix than on a Windows machine.write.csv(grid.pow,  "fine brute-force.csv", row.names = F)
grid.pow.fine <- read.csv("fine brute-force.csv")
best.pow.fine <- grid.pow.fine[which.max(grid.pow.fine$LL),] #global maximum likelhood
best.pow.fine

grid.pow.all <- rbind(grid.pow.coarse[grid.pow.coarse$M0 != 0,], grid.pow.fine)
best.pow.all <- grid.pow.all[which.max(grid.pow.all$LL),] #global maximum likelhood
best.pow.all
out.pow <- output.pow(dat_nonasymp$X, dat_nonasymp$Y, best.pow.brute, Xes_nonasymp$X) # Can't do CIs, because no var-cov matrix - so no way o draw good estimates of parameters. So, no error propagation. Also, this function has a different setup than the others, since it's coming from a grids earch, rather than from a nls fit.

# Compare fits between grid search and nls with M0 fixed at 0.00045
(best.pow.all$r      - coef(fit.pow2)[1])/best.pow.all $r     #  1.3% difference
(best.pow.all$beta   - coef(fit.pow2)[2])/best.pow.all$beta   # -0.9% difference
(best.pow.optim[3]   - coef(fit.pow2)[1])/best.pow.optim[3]   #  2.3% difference
(best.pow.optim[2]   - coef(fit.pow2)[2])/best.pow.optim[2]   #  1.5% difference
(best.pow.optim[3]   - best.pow.all$r)/best.pow.optim[3]      #  1.0% difference
(best.pow.optim[2]   - best.pow.all$beta)/best.pow.optim[2]   #  0.57% difference





#### FIGURES ####

## Figure 1: compare all biomass curves ####

pdf("Figure 1 biomass fits.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
par(las = 1, bty = "n", mar = c(4, 4, 1, 1), oma = c(0, 0, 2, 0), tcl = 0.2, xpd = F, xpd = NA, pty = "m", family = "mono", mgp = c(2.5, 0.5, 0))
layout(matrix(c(1:15), ncol = 3, byrow = T), heights = c(0.3, 1, 1, 1, 1))
#layout.show(15)

# Summary
plot(1:10, 31:40, type = "n", axes = F, ann = F)
text(0.5, 25, expression(paste("Model         ", R^2, "   ", Delta, "AIC", sep = "")), pos = 4)
text(0.5,  15, sprintf("Logistic      %1.2f %4.1f", out.logis$summary[1], 0 ), pos = 4, col = COL[7], font = 2)
text(0.5,   5, sprintf("Gompertz      %1.2f %4.1f", out.gomp$summary[1],  1 ),  pos = 4, col = COL[6])
text(0.5,  -5, sprintf("4-p Logistic  %1.2f %4.1f", out.fpl$summary[1],   1.9),   pos = 4, col = COL[8])
text(0.5, -15, sprintf("Monomolecular %1.2f %4.1f", out.mono$summary[1], 33.2),  pos = 4, col = COL[5])
mtext("Cerastium", font = 3, line = 0, family = "Helvetica")
mtext("Figure 1", adj = 0.0, cex = 1.2, line = 1.5, family = "Helvetica")

plot(1:10, 31:40, type = "n", axes = F, ann = F)
text(0.5, 25, expression(paste("Model         ", R^2, "   ", Delta, "AIC", sep = "")), pos = 4)
text(0.5,  15, sprintf("Power-law     %1.2f %4.1f", out.pow2$summary[1],  0),   pos = 4, col = COL[4], font = 2)
text(0.5,   5, sprintf("Exponential   %1.2f %4.1f", out.exp$summary[1],   12.5),   pos = 4, col = COL[3])
text(0.5,  -5, sprintf("Linear        %1.2f %4.1f", out.lin$summary[1],   25),   pos = 4, col = COL[1])
text(0.5, -15, sprintf("Linear no-int %1.2f %4.1f", out.linno$summary[1], 39.6), pos = 4, col = COL[2])
mtext("Holcus", font = 3, line = 0, at = 12, family = "Helvetica")

plot(1:10, 31:40, type = "n", axes = F, ann = F)
text(0.5, 25, expression(paste("Model         ", R^2, "    ", Delta, "AIC", sep = "")), pos = 4)
text(0.5,  15, sprintf("4-p Logistic  %1.3f %4.1f", out.fpl.log$summary[1],   0),   pos = 4, col = COL[8], font = 2)
text(0.5,   5, sprintf("Gompertz      %1.3f %4.1f", out.gomp.log$summary[1],  13.2),  pos = 4, col = COL[6])
text(0.5,  -5, sprintf("Monomolecular %1.3f %4.1f", out.mono.log$summary[1],  20.2),  pos = 4, col = COL[5])
text(0.5, -15, sprintf("Logistic      %1.3f %4.1f", out.logis.log$summary[1], 55.4), pos = 4, col = COL[7])
par(xpd = F, family = "Helvetica", pty = "s")

# Biomass
plot(dat_asymp$X, dat_asymp$Y, ylim = c(-0.2, 12), xlim = c(0, 200), xlab = "", ylab = "Biomass (g)")
lines(out.mono$rates$times,  out.mono$rates$M,  col = COL[5])     # monomolecular
lines(out.gomp$rates$times,  out.gomp$rates$M,  col = COL[6])     # gompertz
lines(out.logis$rates$times, out.logis$rates$M, col = COL[7])     # Logistic
lines(out.fpl$rates$times,   out.fpl$rates$M,   col = COL[8])     # 4-parameter logistic
abline(h = 0, col = "gray")
abline(h = out.mono$params[1],  col = COL[5], lty = 2)
abline(h = out.gomp$params[1],  col = COL[6], lty = 2)
abline(h = out.logis$params[1], col = COL[7], lty = 2)
abline(h = out.fpl$params[2],   col = COL[8], lty = 2)
mtext("a)", adj = 0.1, line = -1)

plot(dat_nonasymp$X, dat_nonasymp$Y, xlim = c(0, 83), ylim = c(-0.2, 16), xlab = "Days since sowing", ylab = "")
abline(h = 0, col = "gray")
lines(out.lin$rates$times,   out.lin$rates$M,   col = COL[1])    # Linear
lines(out.linno$rates$times, out.linno$rates$M, col = COL[2])    # Linear, no intercept
lines(out.exp$rates$times,   out.exp$rates$M,   col = COL[3])    # Exponential
lines(out.pow2$rates$times,  out.pow2$rates$M,  col = COL[4])    # Power-law fit 1
mtext("b)", adj = 0.1, line = -1)

plot(dat_nonasymp$X, dat_nonasymp$logY, xlim = c(0, 83), xlab = "", ylab = "", col = "black", axes = F)
axis(1); axis(2, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, 5, 15), las =1)
lines(out.mono.log$rates$times,  out.mono.log$rates$M,  col = COL[5])  
lines(out.gomp.log$rates$times,  log(exp(out.gomp.log$rates$M)/1000),  col = COL[6]) # gompertz requires special treatment, since it was fit in terms of milligrams, rather than grams.
lines(out.logis.log$rates$times, out.logis.log$rates$M, col = COL[7])
lines(out.fpl.log$rates$times,   out.fpl.log$rates$M,   col = COL[8])
mtext("c)", adj = 0.1, line = -1)

##AGR
plot(Xes_asymp$X,  out.mono$rates$AGR, xlab = "", ylab = expression(paste("AGR ", (g^-1%.%day^-1))), type = "n", xlim = c(0, 200), ylim = c(0, .15))
lines(Xes_asymp$X, out.mono$rates$AGR,  col = COL[5])  # monomolecular
lines(Xes_asymp$X, out.gomp$rates$AGR,  col = COL[6])  # Gompertz
lines(Xes_asymp$X, out.logis$rates$AGR, col = COL[7])  # Logistic
lines(Xes_asymp$X, out.fpl$rates$AGR,   col = COL[8])  # 4-parameter Logistic
mtext("d)", adj = 0.1, line = -1)

plot(Xes_nonasymp$X,  out.pow2$rates$AGR, xlab = "Days since sowing", ylab = "", type = "n", xlim = c(0, 83), ylim = c(0, 1))
lines(Xes_nonasymp$X, out.lin$rates$AGR,   col = COL[1]) # Linear 
lines(Xes_nonasymp$X, out.linno$rates$AGR, col = COL[2]) # Linear, no intercept
lines(Xes_nonasymp$X, out.exp$rates$AGR,   col = COL[3]) # Exponential
lines(Xes_nonasymp$X, out.pow2$rates$AGR,  col = COL[4]) # power=law
mtext("e)", adj = 0.1, line = -1)

plot(Xes_nonasymp$X,(out.mono.log$rates$AGR), xlab = "", ylab = "", type = "n", xlim = c(0, 83), ylim = c(0, 1))
lines(Xes_nonasymp$X, (out.mono.log$rates$AGR),  col = COL[5])  # monomolecular
lines(Xes_nonasymp$X, (out.gomp.log$rates$AGR)/1000,  col = COL[6])  # Gompertz divide by 1000 to bring it back to units of g, rather than mg
lines(Xes_nonasymp$X, (out.logis.log$rates$AGR), col = COL[7])  # Logistic
lines(Xes_nonasymp$X, (out.fpl.log$rates$AGR),   col = COL[8])  # 4-parameter Logistic
mtext("f)", adj = 0.1, line = -1)

##RGRt
plot(Xes_asymp$X,  out.mono$rates$RGRt, xlab = "", ylab = expression(paste("RGR ", (g%.%g^-1%.%day^-1))), type = "n", xlim = c(0, 200), ylim = c(0, 0.2))
lines(Xes_asymp$X[out.mono$rates$M > 0], out.mono$rates$RGRt[out.mono$rates$M > 0],  col = COL[5])  # monomolecular
lines(Xes_asymp$X, out.gomp$rates$RGRt,  col = COL[6])  # Gompertz
lines(Xes_asymp$X, out.logis$rates$RGRt, col = COL[7]) # Logistic
lines(Xes_asymp$X, out.fpl$rates$RGRt, col = COL[8])   # 4-parameter Logistic
mtext("g)", adj = 0.1, line = -1)

plot(Xes_nonasymp$X,  out.linno$rates$RGRt, xlab = "Days since sowing", ylab = "", type = "n", xlim = c(0, 83), ylim = c(0, 0.35))
lines(Xes_nonasymp$X[out.lin$rates$RGRt > 0], out.lin$rates$RGRt[out.lin$rates$RGRt > 0], col = COL[1])    # Linear # plot only where biomass is positive
lines(Xes_nonasymp$X, out.linno$rates$RGRt, col = COL[2]) # Linear, no intercept
lines(Xes_nonasymp$X, out.exp$rates$RGRt,   col = COL[3]) # Exponential
lines(Xes_nonasymp$X, out.pow2$rates$RGRt,  col = COL[4]) # power=law
mtext("h)", adj = 0.1, line = -1)

plot(Xes_nonasymp$X,  out.mono.log$rates$RGRt, xlab = "", ylab = "", type = "n", xlim = c(0, 83), ylim = c(0, 0.35))
lines(Xes_nonasymp$X, out.mono.log$rates$RGRt,  col = COL[5])  # monomolecular
lines(Xes_nonasymp$X, out.gomp.log$rates$RGRt,  col = COL[6])  # Gompertz
lines(Xes_nonasymp$X, out.logis.log$rates$RGRt, col = COL[7])  # Logistic
lines(Xes_nonasymp$X, out.fpl.log$rates$RGRt,   col = COL[8])  # 4-param Logistic
mtext("i)", adj = 0.1, line = -1)

#RGRm
plot(out.gomp$rates$M,   out.gomp$rates$RGRm, xlab = "", ylab = expression(paste("RGR ", (g%.%g^-1%.%day^-1))), type = "n", ylim = c(0, 0.2))
lines(out.mono$rates$M[out.mono$rates$M > 0], out.mono$rates$RGRm[out.mono$rates$M > 0],   col = COL[5])  # monomolecular
lines(out.gomp$rates$M,  out.gomp$rates$RGRm,   col = COL[6])  # Gompertz
lines(out.logis$rates$M, out.logis$rates$RGRm, col = COL[7]) # Logistic
lines(out.fpl$rates$M,   out.logis$rates$RGRm, col = COL[8]) # 4-param Logistic
mtext("j)", adj = 0.1, line = -1)

plot(out.pow2$rates$M,  out.pow2$rates$RGRm, xlab = "Biomass (g)", ylab = "", type = "n", ylim = c(0, 0.35), xlim = c(0, 16))
lines(out.lin$rates$M[out.lin$rates$M > 0], out.lin$rates$RGRm[out.lin$rates$M > 0], col = COL[1])    # Linear # plot up only where biomass is positive
lines(out.linno$rates$M, out.linno$rates$RGRm, col = COL[2]) # Linear, no intercept
lines(out.exp$rates$M,   out.exp$rates$RGRm,   col = COL[3]) # Exponential
lines(out.pow2$rates$M,  out.pow2$rates$RGRm,   col = COL[4]) # power=law
mtext("k)", adj = 0.1, line = -1)

plot(exp(out.mono.log$rates$M),  out.mono.log$rates$RGRm, xlab = "", ylab = "", type = "n", xlim = c(0, 16), ylim = c(0, 0.35), axes = T)
lines(exp(out.mono.log$rates$M), out.mono.log$rates$RGRm,   col = COL[5])  # monomolecular
lines(exp(out.gomp.log$rates$M)/1000,  out.gomp.log$rates$RGRm,  col = COL[6])  # Gompertz # back convert from miligrams to grams
lines(exp(out.logis.log$rates$M), out.logis.log$rates$RGRm, col = COL[7]) # Logistic
lines(exp(out.fpl.log$rates$M),   out.fpl.log$rates$RGRm,   col = COL[8]) # 4-param logis
mtext("l)", adj = 0.1, line = -1)
dev.off()


## Figure 2: Why classic RGR is bad ####
pdf("Figure 2 classic RGR is bad.pdf", paper = "a4", height = 7/2.54, width = 14/2.54)
par(las = 1, bty = "n", mfrow = c(1, 2), mar = c(3, 3, 0, 0), oma = c(0, 0, 1, 0), tcl = 0.2, pty = "s", family = "Helvetica", mgp = c(2.5, 0.5, 0))
temp <- aggregate(cbind(Y = dat_asymp$Y), by = list(X = dat_asymp$X), mean)
n.times <- nrow(temp)
plot(dat_asymp$X, dat_asymp$Y,pch = 1, ylab = "Biomass (g)", xlab = "", xlim = c(1, 189), ylim = c(0.01, 50), col = "gray", log = "y", axes = F)
axis(1)
axis(2, at = c(0.01, 0.1, 0.5, 1, 5, 10), labels = c(0.01, 0.1, 0.5, 1, 5, 10))
rgr <- log(temp$Y[n.times]/temp$Y[1])/(temp$X[n.times]-temp$X[1])
xes <- c(temp$X[1], temp$X[n.times])
yes <- c(temp$Y[1], temp$Y[n.times])
lines(xes, temp$Y[1] * exp(rgr * xes), lty = 2)
rgr <- diff(log(temp$Y))/diff(temp$X)
for(i in 1:(n.times-1)){
	xes <- c(temp$X[i] - 6, temp$X[i+1] + 6) - temp$X[i]
	lines(xes + temp$X[i], temp$Y[i] * exp(rgr[i] * xes), col = COL[i])
	}
title(xlab = "Days since sowing", outer = T, line = -1.5)
mtext("a)", adj = 0.1, line = -1)
mtext("Cerastium", font = 3, line = -1)
mtext("Figure 2", adj = 0.05, cex = 1.2)

temp <- aggregate(cbind(Y = dat_nonasymp$Y), by = list(X = dat_nonasymp$X), mean)
n.times <- nrow(temp)
options(scipen =5)
plot(dat_nonasymp$X, dat_nonasymp$Y, pch = 1, ylab = "", xlab = "", xlim = c(1, 89), ylim = c(0.0001, 50), col = "gray", log = "y", axes = F)
axis(1); axis(2, at = (c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, 5, 15), las =1)
rgr <- log(temp$Y[n.times]/temp$Y[1])/(temp$X[n.times]-temp$X[1])
xes <- c(temp$X[1], temp$X[n.times])
yes <- c(temp$Y[1], temp$Y[n.times])
lines(xes, temp$Y[1] * exp(rgr * xes), lty = 2)
rgr <- diff(log(temp$Y))/diff(temp$X)
for(i in 1:(n.times-1)){
	xes <- c(temp$X[i]-3, temp$X[i+1] + 3) - temp$X[i]
	lines(xes + temp$X[i], temp$Y[i] * exp(rgr[i] * xes), col = COL[i])
	}
mtext("b)", adj = 0.1, line = -1)
mtext("Holcus", font = 3, line = -1)
dev.off()

## Figure 3: compare timing of measuring RGRs ####
### Primary task illustrated in this figure: propagating error from parameters into estimates of RGR. 
#Before propagating error, we need to check that the parameter profiles are approximately V-shaped (and that the sampling intervals of the parameters are therefore approximately multivariate normal) before proceeding. 
#To do this, we need to fit a seperate nls fit for each species. We have to use nls(), rather than gnls() or nlsList, because profile methods do not exist for those functions.
fit.logis.CER   <- nls(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp_spp, subset = species == "CER")
fit.logis.GER   <- nls(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp_spp, subset = species == "GER")
par(mfrow = c(2, 3), oma = c(3, 2, 3, 1))
plot(profile(fit.logis.CER, alphamax = 0.1))
plot(profile(fit.logis.GER, alphamax = 0.1))
mtext("Confidence intervals based on the profile sum of squares", side = 3, outer = TRUE)
mtext("Cerastium (top row) and Geranium (Bottom row) - confidence levels of 50%, 80%, 90% and 95%", side = 1, outer = TRUE)
## We see that the parofiles are nicely v-shaped
## Now we proceed to propagating error.

# Fit logistic functions for the two species using nlsList
# Look at RGR for two species nutrients, using logistic function
fit.logis4.0   <- nlsList(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp_spp)
out.logis4     <- output.logis.nlsList(fit.logis4.0, times = Xes_annuals_spp$X, CI = T, LOG = F, alpha = 0.05)
fit.logis.nlme <- nlme(fit.logis4.0)
# WHY is the variance and covariance among estimated parameters so much greater in nlme than in nlsList? Because it includes the among-species variance, as well as the within-species variance. thus, the vcov matrix from nlme is inappropriate for comparing among species (or other treatment groups). So, use the nlsList fit.

# compute differences in timing and magnitude of AGR
time.max.AGR.N <- out.logis4$rates[[1]]$times[out.logis4$rates[[1]]$AGR == max(out.logis4$rates[[1]]$AGR)]
time.max.AGR.Y <- out.logis4$rates[[2]]$times[out.logis4$rates[[2]]$AGR == max(out.logis4$rates[[2]]$AGR)]
time.max.AGR.Y-time.max.AGR.N # difference of 45 days
#Magnitude of peak AGR
max.AGR.N <- max(out.logis4$rates[[1]]$AGR)
max.AGR.Y <- max(out.logis4$rates[[2]]$AGR)
(max.AGR.Y-max.AGR.N)/max.AGR.N # increase of 28%


pdf("Figure 3 timing of RGR comparisons.pdf", paper = "a4", width = 14/2.54, height = 14/2.54)
COL.CI <- "#55555530"
par(mfrow = c(2, 2), bty = "n", mar = c(4, 4, 1, 1), oma = c(0, 0, 1, 0), las = 1, family = "Helvetica", pty = "s", tcl = 0.2, mgp = c(2.5, 0.5, 0))
# Biomass over time
plot(dat_asymp_spp$X, dat_asymp_spp$Y, col = COL[7], pch = ifelse(dat_asymp_spp$species =="GER", 3, 1), xlab = "Days since sowing", ylab = "Biomass (g)", ylim = c(0, 25), xlim = c(0, 200))
for(i in 1:nrow(out.logis4$params)){
	polygon(x = c(out.logis4$rates[[i]]$times, rev(out.logis4$rates[[i]]$times)), y = c(out.logis4$rates[[i]]$M.lo, rev(out.logis4$rates[[i]]$M.hi)), col = COL.CI, border = NA)
	lines(out.logis4$rates[[i]]$times, out.logis4$rates[[i]]$M, col = COL[7], lty = i)    # Logistic
	}
legend(5, 23, legend = c("Geranium", "Cerastium"), col = COL[7], lty = c(2, 1), pch = c(3, 1), merge = F)
mtext("a)", adj = 0.1, line = -1)
mtext("Figure 3", adj = 0.05, cex = 1.2)

# AGR on a time basis
plot(out.logis4$rates[[2]]$times, out.logis4$rates[[2]]$AGR, xlab = "Days since sowing", ylab = "AGR (g/day)", type = "n", xlim = c(0, 200), ylim = c(0, 0.2))
for(i in 1:nrow(out.logis4$params)){
	polygon(x = c(out.logis4$rates[[i]]$times, rev(out.logis4$rates[[i]]$times)), y = c(out.logis4$rates[[i]]$AGR.lo, rev(out.logis4$rates[[i]]$AGR.hi)), col = COL.CI, border = NA)
	lines(out.logis4$rates[[i]]$times, out.logis4$rates[[i]]$AGR, col = COL[7], lty = i)    # Logistic
	}
mtext("b)", adj = 0.1, line = -1)

# RGR on a time basis
plot(out.logis4$rates[[2]]$times,  out.logis4$rates[[2]]$RGRt, xlab = "Days since sowing", ylab = "RGR (g/g/day)", type = "n", xlim = c(0, 200), ylim = c(0, 0.1))
for(i in 1:nrow(out.logis4$params)){
	polygon(x = c(out.logis4$rates[[i]]$times, rev(out.logis4$rates[[i]]$times)), y = c(out.logis4$rates[[i]]$RGRt.lo, rev(out.logis4$rates[[i]]$RGRt.hi)), col = COL.CI, border = NA)
	lines(out.logis4$rates[[i]]$times, out.logis4$rates[[i]]$RGRt, col = COL[7], lty = i)    # Logistic
	}
mtext("c)", adj = 0.1, line = -1)

# RGR on a mass basis
plot(out.logis4$rates[[2]]$M, out.logis4$rates[[2]]$RGRm, xlab = "Predicted biomass (g)", ylab = "", type = "n", axes = F, ylim = c(0, 0.1), xlim = c(0, 11))
axis(1)
axis(2, labels = NA)
for(i in 1:nrow(out.logis4$params)){
	polygon(x = c(out.logis4$rates[[i]]$M, rev(out.logis4$rates[[i]]$M)), y = c(out.logis4$rates[[i]]$RGRm.lo, rev(out.logis4$rates[[i]]$RGRm.hi)), col = COL.CI, border = NA)
	lines(out.logis4$rates[[i]]$M, out.logis4$rates[[i]]$RGRm, col = COL[7], lty = i)    # Logistic
	}
mtext("d)", adj = 0.1, line = -1)
dev.off()





## Figure 4: Brute-force grid search ####
pdf("Figure 4 brute-force search.pdf", paper = "a4", width = 7/2.54, height = 21/2.54)
par(family = "Helvetica", bty = "n", las = 1, pty = "s", tcl = 0.2, mfrow = c(3, 1), mar = c(4, 4, 1, 1), oma = c(0, 0, 1, 0), mgp = c(2.5, 0.5, 0))
options(scipen = 5)
Z <- grid.pow.all[grid.pow.all$M0 == best.pow.all$M0,]
X <- unique(Z$beta) # beta
Y <- unique(Z$r)    # r
Z <- matrix(Z$LL, nrow = length(X))
contour(X, Y, Z,  labcex = 1, method = "edge", xlim = c(0.2, 0.8), xlab = expression(beta), ylab = "r", col = "grey20", vfont = c("serif", "plain"))
points(best.pow.all$beta, best.pow.all$r, col = "red", pch = 4, cex = 2, lend = 1)
mtext("Figure 4", adj = 0.05, cex = 1.2, line = 0.7)
mtext("a)", adj = 0.1, line = -1)

Z <- grid.pow.all[grid.pow.all$r == best.pow.all$r,]
Z <- Z[order(Z$M0),]
X <- unique(Z$beta) # beta
Y <- unique(Z$M0)   # M0
Z <- matrix(Z$LL, nrow = length(X))
par(family = "Helvetica", bty = "n", las = 1, pty = "s", tcl = 0.2)
contour(X, log(Y+0.0001), Z, nlevels = 20,  labcex = 1, method = "flattest",  xlab = expression(beta), ylab = expression(M[0]), col = "grey20", vfont = c("serif", "plain"), axes = F)
axis(1)
axis(2, at = log(c(0.0001, 0.001, 0.01, 0.1, 0.5)), labels = c(0.0001, 0.001, 0.01, 0.1, 0.5))
points(best.pow.all$beta, log(best.pow.all$M0+0.0001), col = "red", pch = 4, cex = 2, lend = 1)
mtext("b)", adj = 0.1, line = -1)

Z <- grid.pow.all[grid.pow.all$beta == best.pow.all$beta,]
Z <- Z[order(Z$M0),]
X <- unique(Z$r) # r
Y <- unique(Z$M0)   # M0
Z <- matrix(Z$LL, nrow = length(X))
par(family = "Helvetica", bty = "n", las = 1, pty = "s", tcl = 0.2)
contour(X, log(Y+0.0001), Z, nlevels = 10,  labcex = 1, method = "flattest", xlab = "r", ylab = expression(M[0]), col = "grey20", vfont = c("serif", "plain"), axes = F)
axis(1)
axis(2, at = log(c(0.0001, 0.001, 0.01, 0.1, 0.5)), c(0.0001, 0.001, 0.01, 0.1, 0.5))
points(best.pow.all$r, log(best.pow.all$M0+0.0001), col = "red", pch = 4, cex = 2, lend = 1)
mtext("c)", adj = 0.1, line = -1)
dev.off()






###SUPPLEMENTAL FIGURES: Appendix 3
## Appendix Figure 1: CIs for differences in growth rates among groups ####
pdf("Appendix 3 Figure 1 CIs of differences in growth rates.pdf", paper = "a4", width = 14/2.54, height = 14/2.54)
COL.CI <- "#55555530"
Cg <- expression(italic("Cerastium")*plain(" greater"))
Gg <- expression(italic("Geranium")*plain(" greater"))
par(mfrow = c(2, 2), bty = "n", mar = c(4, 4, 1, 1), oma = c(0, 0, 1, 0), las = 1, family = "Helvetica", pty = "s", tcl = 0.2, mgp = c(2.5, 0.5, 0))
# plot diference in M and its CIs to show when sig diff
plot(dat_asymp_spp$X, dat_asymp_spp$Y, col = COL[7], xlab = "Days since sowing", ylab = "Difference in biomass (g)", ylim = c(-8, 2), xlim = c(0, 200), type = "n")
abline(h = 0, lty = 3)
for(i in 1:length(out.logis4$diffs)){
	polygon(x = c(out.logis4$diffs[[i]]$times, rev(out.logis4$diffs[[i]]$times)), y = c(out.logis4$diffs[[i]]$diffM.lo, rev(out.logis4$diffs[[i]]$diffM.hi)), col = COL.CI, border = NA)
	lines(out.logis4$diffs[[i]]$times, out.logis4$diffs[[i]]$diffM, col = COL[7], lty = i)    # Logistic
	}
text(200, c(0.4, -0.4), c(Cg, Gg), cex = 0.8, pos = 2)
mtext("a)", adj = 0.1, line = -1)
mtext("Appendix 3 Figure 1", adj = 0.05, cex = 1.2)

# plot diference in AGR and its CIs to show when sig diff
plot(out.logis4$rates[[2]]$times, out.logis4$diffs[[1]]$diffAGR,, xlab = "Days since sowing", ylab = expression(paste("Difference in AGR ", (g^-1%.%day^-1))), type = "n", xlim = c(0, 200), ylim = c(-0.2, 0.2))
abline(h = 0, lty = 3)
for(i in 1:length(out.logis4$diffs)){
	polygon(x = c(out.logis4$diffs[[i]]$times, rev(out.logis4$diffs[[i]]$times)), y = c(out.logis4$diffs[[i]]$diffAGR.lo, rev(out.logis4$diffs[[i]]$diffAGR.hi)), col = COL.CI, border = NA)
	lines(out.logis4$diffs[[i]]$times, out.logis4$diffs[[i]]$diffAGR, col = COL[7], lty = i)    # Logistic
	}
	
text(200, c(0.15, -0.15), c(Cg, Gg), cex = 0.8, pos = 2)
mtext("b)", adj = 0.1, line = -1)


# plot diference in RGRt and its CIs to show when sig diff
plot(out.logis4$rates[[2]]$times, out.logis4$diffs[[1]]$diffRGRt, xlab = "Days since sowing", ylab = expression(paste("Difference in RGR ", (g%.%g^-1%.%day^-1))), type = "n", xlim = c(0, 200), ylim = c(-0.04, 0.1))
abline(h = 0, lty = 3)
for(i in 1:length(out.logis4$diffs)){
	polygon(x = c(out.logis4$diffs[[i]]$times, rev(out.logis4$diffs[[i]]$times)), y = c(out.logis4$diffs[[i]]$diffRGRt.lo, rev(out.logis4$diffs[[i]]$diffRGRt.hi)), col = COL.CI, border = NA)
	lines(out.logis4$diffs[[i]]$times, out.logis4$diffs[[i]]$diffRGRt, col = COL[7], lty = i)    # Logistic
	}
text(200, c(0.05, -0.01), c(Cg, Gg), cex = 0.8, pos = 2)
mtext("c)", adj = 0.1, line = -1)

# DIFFERENCE in RGR on a mass basis
plot(out.logis4$rates[[2]]$M, out.logis4$rates[[2]]$RGRm, xlab = "Predicted biomass (g)", ylab = expression(paste("Difference in RGR ", (g%.%g^-1%.%day^-1))), type = "n", axes = F, ylim = c(-0.04, 0.1), xlim = c(0, 11))
axis(1)
axis(2)
abline(h = 0, lty = 3)
for(i in 1:length(out.logis4$diffs)){
	polygon(x = c(out.logis4$diffs[[i]]$Mseq, rev(out.logis4$diffs[[i]]$Mseq)), y = c(out.logis4$diffs[[i]]$diffRGRm.lo, rev(out.logis4$diffs[[i]]$diffRGRm.hi)), col = COL.CI, border = NA)
	lines(out.logis4$diffs[[i]]$Mseq, out.logis4$diffs[[i]]$diffRGRm, col = COL[7], lty = i)    # Logistic
	}
text(11, c(0.035, -0.035), c(Cg, Gg), cex = 0.8, pos = 2)
mtext("d)", adj = 0.1, line = -1)
dev.off()



## Appendix Figure 2: Diagnostic plots ####
pdf("Appendix 3 Figure 2 diagnostics.pdf", paper = "a4", width = 21/2.54, height = (7*2.3)/2.54)
par(pty = "s", las = 1, bty = "n", mar = c(4, 4, 1, 1), tcl = 0.2, oma = c(0, 0, 2, 0), tcl = 0.2, xpd = F, xpd = NA, pty = "m", family = "mono", mgp = c(2.5, 0.5, 0))
layout(matrix(c(1:9), ncol = 3, byrow = T), heights = c(0.3, 1, 1))
#layout.show(10)
par(xpd = NA, pty = "m", family = "mono")
plot(1:10, 31:40, type = "n", axes = F, ann = F)
text(0.5, 25, sprintf("              RMSE"), pos = 4)
text(0.5, 15, sprintf("Monomolecular %1.2f", out.mono$summary[3]),  pos = 4, col = COL[5])
text(0.5,  5, sprintf("Logistic      %1.2f", out.logis$summary[3]), pos = 4, col = COL[7])
text(0.5, -5, sprintf("4-p Logistic  %1.2f", out.fpl$summary[3]),   pos = 4, col = COL[8])
text(0.5,-15, sprintf("Gompertz      %1.2f", out.gomp$summary[3]),  pos = 4, col = COL[6], font = 2)
mtext("Appendix 3 Figure 2", adj = 0.05, cex = 1.2, line = 0.5, family = "Helvetica")
mtext("Cerastium", font = 3, line = -0.5, family = "Helvetica")

plot(1:10, 31:40, type = "n", axes = F, ann = F)
text(0.5, 25, sprintf("              RMSE"), pos = 4)
text(0.5, 15, sprintf("Linear        %1.2f", out.lin$summary[3]), pos = 4, col = COL[1])
text(0.5,  5, sprintf("Linear no-int %1.2f", out.linno$summary[3]), pos = 4, col = COL[2])
text(0.5, -5, sprintf("Exponential   %1.2f", out.exp$summary[3]), pos = 4, col = COL[3])
text(0.5,-15, sprintf("Power-law     %1.2f", out.pow2$summary[3]), pos = 4, col = COL[4], font = 2)
mtext("Holcus", font = 3, line = -0.5, at = 12, family = "Helvetica")

plot(1:10, 31:40, type = "n", axes = F, ann = F)
text(0.5, 25,  sprintf("              RMSE"), pos = 4)
text(0.5, 15,  sprintf("Monomolecular %1.2f", out.mono.log$summary[3]),  pos = 4, col = COL[5])
text(0.5,  5,  sprintf("Logistic      %1.2f", out.logis.log$summary[3]), pos = 4, col = COL[7])
text(0.5,  -5, sprintf("4-p Logistic  %1.2f", out.fpl.log$summary[3]),   pos = 4, col = COL[8], font = 2)
text(0.5, -15, sprintf("Gompertz      %1.2f", out.gomp.log$summary[3]),  pos = 4, col = COL[6])
par(xpd = F, family = "Helvetica", pty = "s")

# predicted and observed masses
par(family = "Helvetica")
plot(out.mono$data$fitted, dat_asymp$Y, ylab = "Observed biomass (g)", xlab = "",  type = "n", asp = 1, xlim = c(-0.5, 8))
abline(0, 1, col = "gray")
points(out.mono$data$fitted, dat_asymp$Y, col = COL[5])    # monomolecular 
points(out.gomp$data$fitted, dat_asymp$Y, col = COL[6])    # Gompertz
points(out.logis$data$fitted,dat_asymp$Y,     col = COL[7])    # Logistic
mtext("a)", adj = 0.1, line = -1)

plot(out.lin$data$fitted, dat_nonasymp$Y, ylab = "", xlab = "", asp = 1, type = "n", xlim = c(0, 16))
abline(0, 1, col = "gray")
points(out.lin$data$fitted,   dat_nonasymp$Y, col = COL[1])        # Linear
points(out.linno$data$fitted, dat_nonasymp$Y, col = COL[2])  # linear, no intercept
points(out.exp$data$fitted,   dat_nonasymp$Y, col = COL[3])  # exponential
points(out.pow2$data$fitted,   dat_nonasymp$Y, col = COL[4])  # Power-law
mtext("b)", adj = 0.1, line = -1)

plot(out.mono.log$data$fitted, dat_nonasymp$logY, ylab = "", xlab = "", asp = 1, type = "n", axes = F)
axis(1, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, NA, 15), las =1)
axis(2, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, 5,  15), las =1)
abline(0, 1, col = "gray")
points(out.mono.log$data$fitted,  dat_nonasymp$logY, col = COL[5])        # 
points(log(exp(out.gomp.log$data$fitted)/1000), dat_nonasymp$logY, col = COL[6])  
points(out.logis.log$data$fitted, dat_nonasymp$logY, col = COL[7])  # 
points(out.fpl.log$data$fitted,   dat_nonasymp$logY, col = COL[8])  # Power-law
mtext("c)", adj = 0.1, line = -1)

# Residual plot
plot(out.mono$data$fitted, out.mono$data$resid, xlab = "", ylab = "Residuals", type = "n", xlim = c(-0.5, 8))
abline(h = 0, col = "gray")
points(out.mono$data$fitted, out.mono$data$resid,  col = COL[5]) # monomolecular
points(out.gomp$data$fitted, out.gomp$data$resid,  col = COL[6]) # Gompertz
points(out.logis$data$fitted,out.logis$data$resid, col = COL[7]) # Logistic
points(out.fpl$data$fitted,  out.fpl$data$resid,   col = COL[8]) # Logistic
mtext("d)", adj = 0.1, line = -1)

# Residual plot
plot(out.lin$data$fitted, out.lin$data$resid, xlab = "Predicted biomass (g)", ylab = "", type = "n", xlim = c(0, 16))
abline(h = 0, col = "gray")
points(out.lin$data$fitted, out.lin$data$resid,   col = COL[1])                 # Linear
points(out.linno$data$fitted, out.linno$data$resid, col = COL[2]) # linear, no intercept
points(out.exp$data$fitted, out.exp$data$resid, col = COL[3])    # exponential
points(out.pow2$data$fitted, out.pow2$data$resid, col = COL[4])    # Power-law
mtext("e)", adj = 0.1, line = -1)

plot(out.mono.log$data$fitted, out.mono.log $data$resid, xlab = "", ylab = "", type = "n", ylim = c(-3, 3), axes = F)
axis(1, at = log(c(0.001, 0.01, 0.1, 1, 5, 15)), labels = c(0.001, 0.01, 0.1, 1, NA, 15), las =1); axis(2)
abline(h = 0, col = "gray")
points(out.mono.log$data$fitted,   out.mono.log$data$resid,   col = COL[5])                 # Linear
points(log(exp(out.gomp.log$data$fitted)/1000), out.gomp.log$data$resid, col = COL[6]) # linear, no intercept
points(out.logis.log$data$fitted,  out.logis.log$data$resid, col = COL[7])    # exponential
points(out.fpl.log $data$fitted,   out.fpl.log $data$resid, col = COL[8])    # Power-law
mtext("f)", adj = 0.1, line = -1)
dev.off()
