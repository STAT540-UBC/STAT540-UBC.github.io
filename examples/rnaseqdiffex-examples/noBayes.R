# Populate slots in limma fit without doing shrinkage.

noBayes <- function(fit)
{
  eb <- nobayes(fit=fit)
  
  fit$t <- eb$t
  fit$df.total <- eb$df.total
  fit$df.prior <- 0 # needed
  fit$p.value <- eb$p.value
  fit$s2.post <- eb$s2.post
  
  if(!is.null(fit$design) && is.fullrank(fit$design)) {
    F.stat <- classifyTestsF(fit,fstat.only=TRUE)
    fit$F <- as.vector(F.stat)
    df1 <- attr(F.stat,"df1")
    df2 <- attr(F.stat,"df2")
    if(df2[1] > 1e6) # Work around bug in R 2.1
      fit$F.p.value <- pchisq(df1*fit$F,df1,lower.tail=FALSE)
    else
      fit$F.p.value <- pf(fit$F,df1,df2,lower.tail=FALSE)
  }
  fit
}

nobayes <- function(fit)
{
  coefficients <- fit$coefficients
  stdev.unscaled <- fit$stdev.unscaled
  sigma <- fit$sigma
  df.residual <- fit$df.residual
  if(is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) stop("No data, or argument is not a valid lmFit object")
  if(all(df.residual==0)) stop("No residual degrees of freedom in linear model fits")
  if(all(!is.finite(sigma))) stop("No finite residual standard deviations")
  
  out <-  list(df.prior=0,var.prior=fit$scale)
  out$s2.prior <- 0
  out$s2.post <- sigma^2
  out$var.prior <- out$var.post <- NULL
  out$t <- coefficients / stdev.unscaled / fit$sigma
  df.total <- df.residual
  df.pooled <- sum(df.residual,na.rm=TRUE)
  df.total <- pmin(df.total,df.pooled)
  out$df.total <- df.total
  out$p.value <- 2*pt(-abs(out$t),df=df.total)
  
  out
}


