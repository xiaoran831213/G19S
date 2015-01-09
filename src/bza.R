library(mgcv)
library(nlme)

## bezier function analysis
BZA <- new.env();

## check and flip genotype codeing to minimize genotype value
## turbulance between two variants
## v1 ---- the reference variant to be compared against
## v2 ---- the target variant to be check and flipped.
## NA is tolarated by ignoring these indivudials' contribution
## to terbulence.
## level ---- the top genotype value, use 2 for dosage genotype
## coded as 0, 1, and 2.
BZA$flp <- function(v1, v2, level = 2L)
{
    ## compliment of v2 <=> flipped v2
    v2.c <- level - v2;

    ## compare the before and after flipping terbulence between
    ## two variants across all individual, if NA turns up at an
    ## individual of either varaint, this individual contribute
    ## no turbulence.
    ## d is in monotune with sum(|v1-v2|-|v1-v2.c|)
    d <- sum((level-v1*2L)*(level-v2.c*2L), na.rm=TRUE);
    
    if(d > 0L)
    {
        return (v2.c);
    }
    else
    {
        return (v2);
    }
}

## flip genotype codeing to reduce turbulance between variants,
## so the fiited curve is smoother.
## gmx -- genotype matrix, row variant, column individual
BZA$gfp <- function(gmx)
{
    ## allocate the output - flipped genotype matrix;
    out <- matrix(data=0L, nrow=nrow(gmx), ncol=ncol(gmx));
    
    ## the first variant doesn't require flipping 
    out[1L,] <- gmx[1L,];

    ## deal with the rest, use previous(i-1 th.) variant as a
    ## reference, check and flip current variant if necessary
    for(i in 2L:nrow(gmx))
    {
        out[i,] <- BZA$flp(out[i-1L,], gmx[i,]);
    }

    out;
}

## Elevate degree of fitting, aka. increase the number
## of guild point by R, while keeping the curve intact.
## Return new guild points with R more members
## than the original set.
BZA$elv <- function(R, P)
{
    ## R --- degree elevation
    ## P --- original guild points

    ## original degree
    D <- length(P);

    ## allocate vector for final P
    P <- c(P, rep_len(0.0, R));

    ## Q[i]=P[i-1]
    Q <- double(D+R);

    ## final degree - 1
    N <- D-1;
    E <- D+R-1;
    while(N < E)
    {
        N <- N+1L;
        ## update coeficient
        C <-(0:E)/N;

        ## Q[i]=P[i-1] (i=2,...,N+1), Q[1]=P[0]=ANY_THING
        Q[-1L] <- P[1L:E];
        P <- C*Q + (1-C)*P;
    }
    P
}

## Get bezier coeficient
## N --- degree of fit, number of guild points
## T --- time to sample the fitted curve, must span
## accross [0,1].
BZA$cof <- function(N, T)
{
    ## number of sample frames
    M <- length(T);
    
    ## B[t,i] = C(N-1,i-1)*t^(i-1)*(1-t)^(N-i)
    B <- matrix(data=double(), nrow=M, ncol=N);
    
    D <- log(1.0-T);
    T <- log(T);

    ## B[t,1] = (1-t)^(N-1)
    ## log(B[t,1]) = (N-1)log(1-t)
    B[,1L] <- (N-1L)*D;

    ## B[t,i+1] = B[t,i]*(N-i)/i*t/(1-t);
    ## log(B[t,i+1]) = log(B[t,i])+log(N-i)-log(i)+log(T)-log(1-T);
    for(i in 1L:(N-1L))
    {
       B[,i+1L] <- B[,i]+log(N-i)-log(i)+T-D;
    }
    B <- exp(B);

    ## When t=1, B[1,i]=NA, since log(1-1)=-Inf, directly
    ## calculate B[1,i] instead,
    ## When i<N, B[1,i]=C(N-1,i-1)*1^(i-1)*(1-1)^(N-i)=0
    ## when i=N, B[1,N]=C(N-1,N-1)*1^(N-1)*(1-1)^(N-N)=1
    if(T[M]==0.0) # t=1 <=> log(t)==0
    {
        B[M, ] <- 0.0;
        B[M,N] <- 1.0;
    }
    B;
}

## eliminate NA by elevate valid points to uniformed
## dimension
## dat --- data to cleanup NAs, in column major, row
##     count decide the final dimension
BZA$clr <- function(dat)
{
    ## points in each data unit
    N <- nrow(dat);
    out <- apply(dat, MARGIN=2L, FUN=function(d)
    {
        d <- d[!is.na(d)];
        BZA$elv(R=N-length(d), P=d);
    });
    out;
}

## fitting bezier curves to genotypes
## gmx ---- genotype matrix.
## pos ---- genotype position.
## mrg ---- margin of matrix.
##    1: row major -- row variant, col individual
##    2: col major -- col variant, row individual
## flp: flop genotype coding to reduce turbulance.
## frm: sample frames
##    a single number mean the number of frames
##    a vector give precise frame posion, which can be
##    non-uniformly intervaled.
## res: sample resolution, only this many frames are
##    sampled from fitted funcionx.
## return function matrix of resoluation res
BZA$fit <- function(gmx, pos, flp=1, frm=0, div=1, res=0)
{
    ## t: time frame to sample the fitted curve, which by
    ## default uses standarized genotype positon
    ## M: number of sample frames, by default equal to the
    ## number of variants.
    if(length(frm)==1L)
    {
        if(frm==0)
        {
            frm <- pos;
        }
        else
        {
            frm <- 1L:frm;
        }
    }
    M <- length(frm);
    t <- (frm-frm[1L])/(frm[M]-frm[1L]);
    
    ## consider sample resolution
    if(res>0L)
    {
        res <- round((0L:(res-1L))*(M-1L)/(res-1L))+1L;
        t <- t[res];
        M <- length(t);
    }
    else
    {
        res <- M;
    }
    
    ## number of guide points
    N <- nrow(gmx);

    ## flip genotype variant coding for smoother curves
    if(flp>0L)
        gmx <- BZA$gfp(gmx);

    ## standardize genotypes and positions to [0,1]
    gmx <- gmx/2;
    pos <- pos-head(pos,1L);
    pos <- pos/tail(pos,1L);
    
    ## expand position vector to matrix, each individual
    ## now has its own variant map.
    ## position information is dropped if a variant has
    ## missing genotype value
    pos <- matrix(pos, nrow=nrow(gmx), ncol=ncol(gmx));
    pos[is.na(gmx)] <- NA;
    
    ## Deal with missing values, get fianl guide points
    gmx <- BZA$clr(gmx);
    pos <- BZA$clr(pos);

    ## fit genotype curve now!
    g <- matrix(double(), nrow=M, ncol=ncol(gmx));
    p <- matrix(double(), nrow=M, ncol=ncol(pos));
    ixM <- 0L;
    for(i in 1L:div)
    {
        ix1 <- ixM+1L;
        ixM <- round(M*i/div);
        b <- BZA$cof(N,t[ix1:ixM]);
        g[ix1:ixM,] <- b %*% gmx;
        p[ix1:ixM,] <- b %*% pos;
    }
    list(gmx=g, pos=p, time=t,
         mrg=1L, frm=frm, res=res, div=div);
}

BZA$pic <- function(r, out=NULL, res=1.0, ...)
{
    if(r$mrg==2L)
    {
        gmx <- t(r$gmx);
        pos <- t(r$pos);
    }
    else
    {
        gmx <- r$gmx;
        pos <- r$pos;
    }

    ## determine plot resolution.
    M <- length(r$time);
    if (res > 1.0)
    {
        res <- ceiling(res);
    }
    else if(res > 0.0)
    {
        res <- ceiling(M*res);
    }
    else if(res > -1.0)
    {
        ## including res = 0.0
        res <- ceiling(M*(1.0+res));
    }
    else
    {
        res <- ceiling(M + res);
    }
    res <- seq(from=0.0, to=1.0, length.out=res);
    res <- as.integer(round(res*(M-1L)))+1L;

    ## pick out points according to resolution.
    gmx <- gmx[res,];
    pos <- pos[res,];

    ## save plot settings
    pardefault <- par(no.readonly = T)

    ## connect to ploting device
    if(!is.null(out))
        png(out);

    ## dummy plot to layout the background
    plot(x=0.0, y=0.0, xlim=c(0,1), ylim=c(0,1));

    ## plot every individual
    for(i in 1L:ncol(gmx))
    {
        ## sample from fitted bezier curve
        par(new=TRUE);
        plot(pos[,i], gmx[,i], type="l", ylab="", xlab="", axes=F);
    }

    ## close ploting device.
    if(!is.null(out))
        dev.off();
    
    ## restore plot settings
    pardefault <- par(pardefault)
}

## get cumulative pairwise distance for a function set
## fmx --- function matrix,
##    row ---- variant
##    col ---- function
##    for a genotype matrix, one column is one sample's genotype
## ... additional function matrix, must be indentical in dimension
## time --- time frame to take sample from the functions.
BZA$dff <- function(time, fmx, ..., norm=NULL)
{
    ## list of fitted functions
    lsf <- list(fmx, ...);
    
    ## number of individuals
    S <- ncol(fmx);

    ## number of sample time
    M <- length(time);
    
    ## all pair combinations
    C <- combn(S,2L);

    ## helper indices
    p <- 2L:M;
    q <- p-1L;

    # M-1 sample intervals
    dt <- time[p]-time[q];

    ## iterate individual pairs, fill up pairwise distance matrix
    dm <- matrix(0, S, S);
    for(k in 1L:ncol(C))
    {
        i <- C[1L,k];
        j <- C[2L,k];

        ## get square distance between individual i and j
        df <- 0;
        for(f in lsf)
        {
            df <- df + (f[,i]-f[,j])^2;
        }
        
        ## integrate over all sample intervals
        dm[i,j] <- sum((df[p]+df[q])*dt)/2;
        dm[j,i] <- dm[i,j];
    }

    if(!is.null(norm))
        dm <- norm(dm);
    dm;
}

## scale to [0,1] normalization
BZA$norm.01 <- function(x)
{
    r <- range(x);
    (x-r[1L])/(r[2L]-r[1L]);
}

## rank normal quantile normalization
BZA$norm.rq <- function(x)
{
    qnorm((rank(x)-0.5)/length(x));
}
