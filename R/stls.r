funcval.STLS <- function(bet,x,y) 		
{
	funcvalue <- sum((y-pmax((1/2)*y, x%*%bet))^2)
	return(funcvalue)
}

stls.fit<-function (mf, x, y, point, direction, bet, ...) 
{
    	dots <- list(...)
    	if (is.null(dots$control)) 
        	control <- list(maxit = 2000)
    	else control <- dots$control
    	z <- optim(par = bet, fn = funcval.STLS, control = control, 
        	x = x, y = y)
	    z$counts <- z$counts[1]
    	if (z$counts < 10) 
        	warning("Convergence reached after extremely few iterations. Make sure the specifications \n in the function call are correct (point, direction, starting values etc.).")
   	  b <- z$par
    	
    	z$residuals <- y - x %*% b
    	
    	if (direction == "right")
      { 
        	z$par <- -z$par
        	bet <- -bet
     	}
    
    	if (point != 0)
      { 
        	z$par[1, 1] <- z$par[1, 1] + point
          bet[1, 1] <- bet[1, 1] + point
      }
    
    	z$startcoef <- bet
    	rownames(z$startcoef) <- c("(Intercept)", names(mf[-1]))
    	colnames(z$startcoef) <- c("")
    	dimnames(z$par) <- dimnames(z$startcoef)
    	z$coefficients <- t(z$par)
    	z$df.residual <- length(mf[, 1]) - ncol(x)
    	z$fitted.values <- x %*% z$par
    	return(z)
}



stls <- function (formula, data, point = 0, direction = "left", beta = "ols", covar = FALSE, na.action, ...) 
{
    	cll <- match.call()
    	mf <- match.call(expand.dots = FALSE)
    	m <- match(c("formula", "data", "na.action"), names(mf), 0L)
    	mf <- mf[c(1L, m)]
    	mf$drop.unused.levels <- TRUE
    	mf[[1L]] <- as.name("model.frame")
    	mf <- eval(mf, parent.frame())
    	if (point != 0)
        	mf[, 1] <- mf[, 1] - point
    	
    	if (direction == "right") 
        	mf[, 1] <- -mf[, 1]
   	 
    	y <- model.response(mf, "numeric")
    	y <- matrix(y)
    	x <- model.matrix(formula, mf)
	if (is.numeric(beta)) 
	{
	 bet <- matrix(beta, ncol(x), 1)
 	  if (point != 0) 
    {
   	  bet[1, 1] <- bet[1, 1] - point
   	}
    if (direction == "right") 
    {
   	  bet <- -bet
   	}
	}
	else 
	{
		if (beta == "ols") 
		{
			bet <- lm(data = mf)$coef
			bet <- matrix(bet)
		}
		if (beta == "ml") 
		{
			mlcof <- mlcoef(mf)
			bet <- matrix(mlcof[-length(mlcof)])
  	}
  	if (beta != "ols" && beta != "ml") 
  	 stop("'beta' must be numeric or a valid method.")
	}

	z <- stls.fit(mf, x, y, point, direction, bet, ...)
	
	if (covar) 
	{
		dots <- list(...)
		if (is.null(dots$R)) 
			R <- 200
		else R <- dots$R
		if (is.null(dots$control)) control <- list(maxit=2000) else control <- dots$control
		bootobj <- covar.bootSTLS(mf, x, funcSTLS, R = R, beta, bet, point, direction, control)
		z$covariance <- bootobj[[1]]
		rownames(z$covariance) <- c("(Intercept)", names(mf[-1]))
		colnames(z$covariance) <- c("(Intercept)", names(mf[-1]))
		z$bootrepl <- bootobj[[2]]
	}
	class(z) <- c("stls")
	z$call <- cll
z
}
