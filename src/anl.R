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
ANL$one <- function(gmx, emx, rsp, cvr, gsm='dist')
{
    ## return list
    out <- list();
    
    ## join genotype into one dimensional variable
    g <- HWU$collapse.burden(gmx);

    ## combined test, weight matrix is transcriptome similarity
    k<-HWU$weight.gaussian(emx);
    r<-HWU$run(Tn=rsp, geno=g, X=cvr, K=k, gsim = gsm, appx = 'davis');
    out$GT <- r$p;

    ## only test genotype, weight matrix is I
    i<-matrix(1, nrow(gmx), nrow(gmx));
    r<-HWU$run(Tn=rsp, geno=g, X=cvr, K=i, gsim = gsm, appx = 'davis');
    out$G <- r$p;

    ## test for existance of transcriptome(heterogeneity) effect, normalize
    ## previous transcriptome similarity matrix to zero center
    k<-k-mean(k);
    r<-HWU$run(Tn=rsp, geno=g, X=cvr, K=k, gsim = gsm, appx = 'davis');
    out$T <- r$p;

    out;
}

ANL$bza <- function(gmx, pos, emx, rsp, cvr)
{
    out <- list();
    bz <- BZA$fit(gmx, pos);
    df <- BZA$dff(f$gmx, f$pos);

    wg <- exp(-df);
    we <- HWU$weight.gaussian(emx);
    out$GT <- dg2(rsp, cvr, wg, we);

    we <- we-mean(we);
    out$T <- dg2(rsp, cvr, wg, we);

    out$G <- dg2(rsp, cvr, wg);

    out;
}
    
ANL$run1<-function(gno, exp, rng, phe, rsp, cov, pcs=NULL, imp=T, FUN=ANL$one, ...)
{
    phe<-phe[, c("IID", rsp, cov)];
    phe<-HLP$clrMss(phe);

    ## align tables by indidivual id
    tmp<-HLP$algIdv(gno, exp, phe);
    gdx<-tmp$gdx;
    edx<-tmp$edx;
    pdx<-tmp$pdx;
    
    ## integrity check
    stopifnot(identical(phe$IID[pdx], gno$idv[gdx]));
    stopifnot(identical(phe$IID[pdx], exp$idv[edx]));
    
    ## response and covariate
    rsp<-as.matrix(phe[pdx, rsp], rownames.force = F);   # response variable
    cov<-as.matrix(phe[pdx, cov], rownames.force = F);   # covariants
    if(!is.null(pcs))
    {
        cov <- cbind(cov, pcs[gdx,]);
    }
    
    ## prepare output containers
    nrg <- nrow(rng);
    out <- matrix(data = NA, nrow = nrg, ncol = 3L); # 3 types of test
    nvr <- rep.int(NA, nrg); # number of variants in gene ranges
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
        chr<-rng[i,CHR];
        bp1<-rng[i,BP1]-5001L;
        bp2<-rng[i,BP2]+5001L;
        idx<-gno$map[CHR==chr][bp1<POS&POS<bp2, IDX];

        ## record number of variants in the i.th gene range.
        nvr[i] <- length(idx); 
        if(nvr[i]==0)
        {
            err[i]<-'NVR=0' # no variants in gene i.
            next;
        }

        ## row->variant, col->individual
        gmx<-gno$gmx[idx,, drop=F]; 
        stopifnot(identical(gmx, HLP$gld(rng[i,])$gmx));

        ## save population mean genotype value for later imputation
        avg<-apply(gmx, 1L, mean, na.rm=T); 
        
        ## exclude unavailable individual
        gmx<-gmx[,gdx, drop=F];  # row->variant, col->individual
        
        ## exclude degenerated variants
        idx<-HLP$gmxClr(gmx);
        if(length(idx)==0L)
        {
            err[i]<-'NES=0' ## number of effect sample is zero
            next;
        }
        gmx<-gmx[idx,,drop=F];
        
        ## asign population mean genotype value to NA if imputation
        ## is requested.
        if(imp)
        {
            for(j in 1L:nrow(gmx))
                gmx[j, is.na(gmx[j,])]<-avg[j];
        }
        
        ## call U sta
        gmx<-t(gmx); # row->individual, col->genomic variant
        emx<-t(emx); # row->individual, col->RNA probe
        r<-FUN(gmx, emx, rsp, cov, ...);
        if (inherits(r, "try-error"))
        {
            err[i] <- errMsg(p);
            next;
        }
        
        ## record result.
        out[i,] <- c(r$GT.p, r$G.p, r$T.p);

        ## report progress.
        HLP$shwPrg(nrg, i, 256L);
    }
    out<-data.table(
        rng, 
        GT=out[,1L],
        G=out[,2L],
        T=out[,3L],
        MIN=pmin(out[,1L], out[,2L], out[,3L]),
        NVR=nvr,
        ERR=err, key='CHR,BP1,BP2');
    out;
}

ANL$run2<-function(exp, rng, phe, rsp, cov, pcs=NULL, imp=T, FUN=ANL$one, ...)
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
    if(!is.null(pcs))
    {
        cov <- cbind(cov, pcs[gdx,]);
    }
    
    ## prepare output containers
    nrg <- nrow(rng);
    out <- matrix(data = NA, nrow = nrg, ncol = 3L); # 3 types of test
    out <- list();
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
        gmx<-gmx[idx,,drop=F];
        
        ## asign population mean genotype value to NA if imputation
        ## is requested.
        if(imp)
        {
            for(j in 1L:nrow(gmx))
                gmx[j, is.na(gmx[j,])]<-avg[j];
        }
        
        ## call U sta
        gmx<-t(gmx); # row->individual, col->genomic variant
        emx<-t(emx); # row->individual, col->RNA probe
        r<-FUN(gmx, emx, rsp, cov, ...);
        if (inherits(r, "try-error"))
        {
            err[i] <- errMsg(p);
            next;
        }
        
        ## record result.
        ## out[i,] <- c(r$GT.p, r$G.p, r$T.p);
        out[[length(out)+1L]] <- r;

        ## report progress.
        HLP$shwPrg(nrg, i, 256L);
    }
    out <- HLP$mktab(out);
    out<-data.table(
        rng,
        GT=out[,1L], G=out[,2L], T=out[,3L],
        MIN=pmin(out[,1L], out[,2L], out[,3L]),
        NVR=nvr, ERR=err, key='CHR,BP1,BP2');
    out;
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

