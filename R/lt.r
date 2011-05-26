funcval.LT <- function(bet,x,y,cl,cu) 		
{
	ind <- ifelse(x%*%bet>cl,(ifelse(-cl<=(y-x%*%bet)&(y-x%*%bet)<=cu,0.5*(y-x%*%bet)^2,0)+ifelse((y-x%*%bet)<(-cl),0.5*cl^2,0)
		+ifelse((y-x%*%bet)>cu,0.5*cu^2,0)),(ifelse(-cl<=(y-cl)&(y-cl)<=cu,0.5*(y-cl)^2,0)+ifelse((y-cl)<(-cl),0.5*cl^2,0)
		+ifelse((y-cl)>cu,0.5*cu^2,0)))
	funcvalue <- sum(ind)
	return(funcvalue)
}


lt.fit<-function (mf, x, y, point, direction, bet, cl, cu, ...) 
{
    	dots <- list(...)
    	if (is.null(dots$control)) 
        	control <- list(maxit = 2000)
    	else control <- dots$control
    	z <- optim(par = bet, fn = funcval.LT, control = control, 
        	x = x, y = y, cl = cl, cu=cu)
   	  z$counts <- z$counts[1] 
      if(z$counts < 10) 
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
        	z$par[1, 1] <- z$par[1, 1] + (point)
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



lt <- function (formula, data, point = 0, direction = "left", clower = "ols", const=1,	 
    cupper = 2, beta = "ols", covar = FALSE, na.action, ...) 
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

	if(!is.numeric(cupper))
		stop("'cupper' must be numeric.")

	if(length(cupper)!=1)
		stop("'cupper' must be a number or a numeric vector of length 1.")
	
	if (is.numeric(clower)) 
	{
		if(length(clower)!=1)
			stop("'clower' must be a number, a numeric vector of length 1, or a valid method.")
		cl <- clower
		cu <- cupper		
	}
	else 
	{
		if (clower == "ols") 
		{
			lmOLS <- lm(data = mf)
			cl <- sqrt(deviance(lmOLS)/df.residual(lmOLS))
			cl <- const*cl
			cu <- cupper*cl
		}
		if (clower == "ml") 
		{
			mlcof <- mlcoef(mf)
			cl <- mlcof[length(mlcof)]
			cl <- const * cl
			cu <- cupper*cl
  	}
  	if (clower != "ols" && clower != "ml") 
  	 stop("'clower' must be numeric or a valid method.")
	}
	z <- lt.fit(mf, x, y, point, direction, bet, cl, cu, ...)
	if (covar) 
	{
		dots <- list(...)
		if (is.null(dots$R)) 
			R <- 200
		else R <- dots$R
		if (is.null(dots$control)) 
        	control <- list(maxit = 2000)
    	else control <- dots$control
		bootobj <- covar.bootLT(mf, x, funcLT, R = R, beta, bet, clower, const, cupper, point, direction, control)
		z$covariance <- bootobj[[1]]
		rownames(z$covariance) <- c("(Intercept)", names(mf[-1]))
		colnames(z$covariance) <- c("(Intercept)", names(mf[-1]))
		z$bootrepl <- bootobj[[2]]
	}
	class(z) <- c("lt")
	z$call <- cll
	z$cvalues <- cbind(cl,cu)
	colnames(z$cvalues) <- c("cl", "cu")
z
}
