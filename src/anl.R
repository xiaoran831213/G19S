source('src/hwu.R')
source('src/hlp.R')
source('src/bza.R')

ANL <- new.env();

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
ANL$bza <- function(gmx, pos, emx, rsp, cvr, ...)
{
    out <- list();
    f <- BZA$fit(gmx, pos, frm=0L);
     
    w1 <- HWU$weight.gaussian(t(f$gmx));
    w2 <- HWU$weight.gaussian(t(f$pos));
    
    out$bgp <- dg2(rsp, cvr, w1*w2);

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
    
ANL$go <- function(exp, rng, phe, rsp, cov, pcs=NULL, imp=T, FUN=ANL$hwu, ...)
{
    ## get individual list from the first gene record.
    rng <- rng[, list(SEQ, CHR, BP1, BP2, GEN, PRB)];
    idv <- HLP$gld(rng[1L,])$idv;

    ## extract phenotype and covariants.
    phe<-phe[, c("IID", rsp, cov)];
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
    rsp<-as.matrix(phe[pdx, rsp], rownames.force = F);   # response variable
    cov<-as.matrix(phe[pdx, cov], rownames.force = F);   # covariants

    ## incroporate PCAs into covariants matrix.
    if(!is.null(pcs))
    {
        cov <- cbind(cov, pcs[gdx,]);
    }
    
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
        
        r<-FUN(gmx=gmx, emx=emx, rsp=rsp, cvr=cov, pos=pos, ...);
        if (inherits(r, "try-error"))
        {
            err[i] <- errMsg(p);
            next;
        }
        
        ## record result.
        out[[i]] <- r;

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

