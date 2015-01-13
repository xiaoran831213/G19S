source('src/hwu.R')
source('src/hlp.R')
source('src/bza.R')

ANL <- new.env();

## f ---- U kernel
## r ---- residual matrix
## w1, ... ---- weight terms.
ANL$dg2 <- function(f, r, w1, ...)
{
    ## get product of all weight terms.
    w <- Reduce(f='*', x=list(w1, ...));
    diag(w) <- 0; # ??
    
    ## compute U score
    u <- sum(w*f);

    ## exclude coveriant effect on both dimensions of w. Residual of W is RWR',
    ## since R is symmetric, RWR' equals RWR
    w <- r %*% w
    w <- w %*% r;

    ## calculate p-value of u.
    r <- try(eigen(w, symmetric=T, only.values=T));
    if (inherits(r, "try-error"))
    {
        r <- -1;
    }
    else
    {
        r <- r$values;
        r <- davies(u, r, acc=0.000001)$Qq;
    }
    r;
}

## gmx --- genomic matrix
##     row: individual
##     col: genome variant
## emx --- expression matrix
##     row: individual
##     col: transcriptome variant
## rsp --- response variable, phenotype
## cvr --- covariants
##     row: individual
##     col: varaint
ANL$bza <- function(gmx, pos, emx, f, r, ...)
{
    out <- list();
    bz <- BZA$fit(gmx, pos);
    ##dm <- BZA$dff(t=bz$time, f=bz$gmx, bz$pos, sqr=T)
    
    w1 <- HWU$weight.gaussian(t(bz$gmx));
    w2 <- HWU$weight.gaussian(t(bz$pos));
    out$G <- ANL$dg2(f, r, w1, w2);

    ## g <- HWU$collapse.burden(t(gmx));
    ## wg <- HWU$weight.gaussian(g);
    ## out$G <- dg2(rsp, cvr, wg);
    
    ## we <- HWU$weight.gaussian(t(emx));
    ## out$J <- dg2(rsp, cvr, wg, we);

    ## we <- we-mean(we);
    ## out$T <- dg2(rsp, cvr, wg, we);

    ## out$MIN <- Reduce(min, out);
    out;
}
    
ANL$go <- function(exp, rng, phe, rsp, cvr, pcs=NULL, imp=F, ...)
{
    ## get individual list from the first gene record.
    rng <- rng[, list(SEQ, CHR, BP1, BP2, GEN, PRB)];
    idv <- HLP$gld(rng[1L,])$idv;

    ## extract phenotype and covariants.
    phe<-phe[, c("IID", rsp, cvr)];
    phe<-HLP$clrMss(phe);

    ## align tables by indidivual id
    iid <- sort(Reduce(intersect, list(idv, exp$idv, phe$IID)));
    gdx <- match(iid, idv);
    edx <- match(iid, exp$idv);
    pdx <- match(iid, phe$IID);
    
    ## integrity check
    stopifnot(identical(idv[gdx], phe$IID[pdx]));
    stopifnot(identical(idv[gdx], exp$idv[edx]));
    
    ## response and covariate
    y<-as.matrix(phe[pdx, rsp], rownames.force = F);   # response variable
    M <- length(y);

    ## standardize y to m=0, s=1
    y <- rank(y);
    y <- (y-mean(y))/sd(y);
    
    ## covariants
    x<-as.matrix(phe[pdx, cvr], rownames.force = F);
    x <- cbind(1, x);       # intercept
    if(!is.null(pcs))       # PCAs
    {
        x <- cbind(x, pcs[gdx,]);
    }
    
    ## regression residual matrix, R = I - X(X'X)^X'
    r <- diag(1, M, M) - tcrossprod(x %*% solve(crossprod(x)), x);
    
    ## exclude liner covariant effect on y, leave residual of Y
    y <- r %*% y;
    y <- y/sqrt(sum(y^2)/(M-ncol(x)));
 
    ## the U kernel is the pair wise similarity between phenotypes
    f <- tcrossprod(y);
    
    ## prepare output containers
    nrg <- nrow(rng);
    out <- as.list(rep(NA, times=nrg));
    nvr <- rep.int(0L, nrg); # number of variants in gene ranges
    err <- rep.int(NA, nrg); # error recorde
    
    ## iterate genome variants
    require(data.table);
    for(i in 1L:nrg)
    {
        ## get expression matrix
        prb<-rng[i,PRB];
        emx<-exp$emx[exp$prb==prb, edx, drop=F]; # row->probe, col->indvidual
        if(length(emx)==0L)
        {
            err[i] <- 'EXP=N'; ## no expression data
            next;
        }
        
        ## get testing range i and pick variants within
        gen <- HLP$gld(rng[i,]);
        gmx <- gen$gmx;
        pos <- gen$map$POS;
        
        ## record number of variants in the i.th gene range.
        if(is.null(gmx))
        {
            err[i]<-'NVR=0' # no variants in gene i.
            next;
        }
        nvr[i] <- nrow(gmx);

        ## save population mean genotype value for later imputation
        avg<-apply(gmx, 1L, mean, na.rm=T); 
        
        ## exclude unavailable individual
        gmx<-gmx[,gdx, drop=F];  # row->variant, col->individual
        
        ## exclude degenerated variants
        idx<-HLP$gmxClr(gmx);
        if(length(idx)==0L)
        {
            err[i]<-'NES=0' ## number of effective sample is zero
            next;
        }
        gmx <- gmx[idx,,drop=F];
        pos <- pos[idx];
        
        ## asign population mean genotype value to NA if imputation
        ## is requested.
        if(imp)
        {
            for(j in 1L:nrow(gmx))
                gmx[j, is.na(gmx[j,])]<-avg[j];
        }
        
        o <- ANL$bza(gmx, pos, emx, f, r);
            ##FUN(gmx=gmx, emx=emx, rsp=rsp, cvr=cov, pos=pos, ...);
        if (inherits(o, "try-error"))
        {
            err[i] <- errMsg(o);
            next;
        }
        
        ## record result.
        out[[i]] <- o;

        ## report progress.
        HLP$shwPrg(nrg, i, 100L);
    }
    data.table(rng, HLP$mktab(out), NVR=nvr, ERR=err, key=key(rng));
}

ANL$qvl <- function(out)
{
    i <- which(is.na(out$ERR));
    out$GT[i] <- p.adjust(out$GT[i], 'fdr');
    out$G[i] <- p.adjust(out$G[i], 'fdr');
    out$T[i] <- p.adjust(out$T[i], 'fdr');
    out$MIN <- pmin(out$GT, out$G, out$T);
    out;
}

