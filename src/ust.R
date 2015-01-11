UST <- new.env();

## x ---- design matrix for covariates
## y ---- response variable
dg2 <- function(y, x=NULL, w1, ...)
{
    ## number of samples
    M <- length(y);
    
    ## prepare regression residual solution matrix based
    ## on covariants, R = I - X(X'X)^X'
    x <- cbind(1, x); # append intercept term
    r <- diag(M) - tcrossprod(x %*% solve(crossprod(x)), x);
    
    ## standardize y to m=0, s=1
    y <- rank(y);
    y <- (y-mean(y))/sd(y);
    
    ## exclude liner effect of coveriant on y
    ## residual of Y = RY
    y <- r %*% y;
    y <- y/sqrt(sum(y^2)/(M-ncol(x)));
 
    ## the U kernel is the pair wise similarity between
    ## all y.i and y.j
    f <- tcrossprod(y);
    
    ## get product of all weight terms.
    w <- Reduce(f='*', x=list(w1, ...));
    diag(w) <- 0; # ??
    
    ## compute U score
    u <- sum(w*f);

    ## exclude coveriant effect on both dimensions of w.
    ## square residual of W is RWR', since R is symmetric,
    ## RWR' equals RWR
    w <- r %*% w
    w <- w %*% r;

    ## calculate p-value of u.
    coef <- eigen(w, symmetric=T, only.values=T)$values;
    p <- davies(u, coef, acc=0.000001)$Qq;
    p
}
