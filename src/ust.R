UST <- new.env();

## x ---- design matrix for covariates
## y ---- response variable
dg2 <- function(y, x=NULL, w, ...)
{
    ## number of samples
    M <- length(y);
    
    ## prepare regression matrix based on covariants
    ## r = x(x'x)^x'
    x <- cbind(1, x); # append intercept term
    r <- tcrossprod(x %*% solve(crossprod(x)), x);

    ## standardize y to m=0, s=1
    y <- rank(y);
    y <- (y-mean(y))/sd(y);

    ## exclude liner effect of coveriant on y
    y <- y - r %*% y;

    ## the U kernel is the pair wise similarity between
    ## all y.i and y.j
    f <- tcrossprod(y);

    ## get product of all weight terms.
    w <- Reduce(f='*', x=list(w, ...));
    diag(w) <- 0; # ??
    
    ## compute U score
    u <- sum(a*f);

    ## exclude coveriant and intercept effect on both
    ## dimensions of w
    w <- tcrossprod(r %*% w, r);

    ## calculate p-value of u.
    coef <- eigen(w, symmetric=T, only.values=T)$values;
    p <- davies(u, coef, acc=0.000001)$Qq;
    p
}
