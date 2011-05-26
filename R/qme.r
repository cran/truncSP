funcval.QME <- function(bet,x,y,cv)
{
	max.xtb.cval <- pmax(x%*%bet, cv)
	ind <- ifelse(abs(y-max.xtb.cval)<cv,((y-max.xtb.cval)^2-cv^2),0)
	funcvalue <- sum(ind)
	return(funcvalue)
}

mlcoef <- function(mf)
{
	mlc<-truncreg(as.formula(mf), data=mf, point=0, direction="left")$coef
	return(mlc)
}



qme.fit <- function(mf,x,y,point,direction,bet,cv,...)
{
	dots <- list(...)

	if (is.null(dots$control)) control <- list(maxit=2000) else control <- dots$control

	z <- optim(par=bet,fn=funcval.QME,control=control,x=x,y=y,cv=cv)

  z$counts <- z$counts[1] 
	if(z$counts < 10)
		warning("Convergence reached after extremely few iterations. Make sure the specifications \n in the function call are correct (point, direction, starting values etc.).")
	
	b <- z$par
	
	z$residuals <- y-x%*%b
	
	if(direction=="right")
	{
		z$par <- -z$par
    bet <- -bet
  }
  
	if(point!=0)
	{
		z$par[1,1] <- z$par[1,1] + point
		bet[1,1] <- bet[1,1] + point
	}
	
	z$startcoef <- bet
	rownames(z$startcoef) <- c("(Intercept)",names(mf[-1]))
	colnames(z$startcoef) <- c("")
	dimnames(z$par) <- dimnames(z$startcoef)
	z$coefficients <- t(z$par)
	z$df.residual <- length(mf[,1])-ncol(x)	
	z$fitted.values <- x%*%z$par
	return(z)
}


qme <- function(formula,data,point=0,direction="left",cval="ols", const=1,beta="ols",covar=FALSE,na.action,...)
{
	cll <- match.call()

	mf <- match.call(expand.dots = FALSE)	
	m <- match(c("formula", "data", "na.action"), names(mf), 0L)
	mf <- mf[c(1L, m)]				
	mf$drop.unused.levels <- TRUE			
	mf[[1L]] <- as.name("model.frame") 		
	mf <- eval(mf, parent.frame())
		
	if(point!=0)
		mf[,1] <- mf[,1]-point

	if(direction=="right")
		mf[,1] <- -mf[,1]

	y <- model.response(mf, "numeric") 	
	y <- matrix(y)
	x <- model.matrix(formula,mf)
	
	if(is.numeric(beta))
	{
		bet <- matrix(beta, ncol(x),1)
		if(point!=0)
		{
			bet[1,1] <- bet[1,1]-point
		}

		if(direction=="right")
		{
			bet <- -bet
		}
	}
	else
	{
		if(beta=="ols")
		{
			bet <- lm(data=mf)$coef
			bet <- matrix(bet)
		}

		if(beta=="ml")
		{
			mlcof<-mlcoef(mf)
			bet <- matrix(mlcof[-length(mlcof)])
		}
		
		if(beta!="ols" && beta!="ml")		
			stop("'beta' must be numeric or a valid method.")
	}

	
	if(is.numeric(cval))
	{
		cv <- cval
		if(length(cval)!= 1)
			stop("'cval' must be one number or a vector of length 1.")
	}
	
	else
	{
		if(cval=="ols")
		{
			lmOLS <- lm(data=mf)
			cv <- sqrt(deviance(lmOLS)/df.residual(lmOLS))
			cv <- const*cv
		}

		if(cval=="ml")
		{
			mlcof<-mlcoef(mf)
			cv <- mlcof[length(mlcof)]
			cv <- const*cv	
		}

		if(cval!="ols" && cval!="ml")
			stop("'cval' must be numeric or a valid method.")
	}

	z <- qme.fit(mf,x,y,point,direction,bet,cv,...)
	
	if(covar)
	{
		dots <- list(...)
		if(is.null(dots$R)) R <- 200	
		else R <- dots$R
		if (is.null(dots$control)) control <- list(maxit=2000) else control <- dots$control
		bootobj <- covar.bootQME(mf,x,funcQME,R=R,beta,bet,cval,const,point,direction,control)
		z$covariance <- bootobj[[1]]
		rownames(z$covariance) <- c("(Intercept)",names(mf[-1]))
		colnames(z$covariance) <- c("(Intercept)",names(mf[-1]))
		z$bootrepl <- bootobj[[2]]
	}

	class(z) <- c("qme")
	z$call <- cll
	z$cval <- cv
	names(z$cval) <- c("cval")
z
}