source('src/hwu.core.R')
source('src/helper.R')

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
ANL$one<-function(gmx, emx, rsp, cvr)
{
    ## return list
    out <- list();
    
    ## join genotype into one dimensional variable
    g <- HWU$collapse.burden(gmx);

    ## only test genotype, weight matrix is I
    k<-matrix(1, nrow(gmx), nrow(gmx));
    r<-HWU$run(Tn=rsp, geno=g, X=cvr, K=k, gsim = 'add', appx = 'davis');
    out$G.u <- r$u;
    out$G.p <- r$p;

    ## combined test, weight matrix is transcriptome similarity
    k<-HWU$weight.gaussian(emx);
    r<-HWU$run(Tn=rsp, geno=g, X=cvr, K=k, gsim = 'add', appx = 'davis');
    out$GT.u <- r$u;
    out$GT.p <- r$p;

    ## test for existance of transcriptome(heterogeneity) effect, normalize
    ## previous transcriptome similarity matrix to zero center
    k<-k-mean(k);
    r<-HWU$run(Tn=rsp, geno=g, X=cvr, K=k, gsim = 'add', appx = 'davis');
    out$T.u <- r$u;
    out$T.p <- r$p;

    out;
}

ANL$run<-function(gno, exp, rng, phe, rsp, cov)
{
    phe<-phe[, c("IID", rsp, cov)];
    phe<-hlp$clrMss(phe);

    ## align tables by indidivual id
    tmp<-hlp$algIdv(gno, exp, phe);
    gdx<-tmp$gdx;
    edx<-tmp$edx;
    pdx<-tmp$pdx;
    
    ## integrity check
    stopifnot(identical(phe$IID[pdx], gno$idv[gdx]));
    stopifnot(identical(phe$IID[pdx], exp$idv[edx]));
    
    ## response and covariate
    rsp<-as.matrix(phe[pdx, rsp], rownames.force = F);   # response variable
    cov<-as.matrix(phe[pdx, cov], rownames.force = F);   # covariants
    
    ## prepare output
    nrg<-nrow(rng);
    out<-matrix(data = NA, nrow = nrg, ncol = 6L);
    err<-rep.int(NA, nrg);
    
    ## iterate genome variants
    for(i in 1L:nrg)
    {
        if(i %% 0xFF == 0x01L)
        {
            cat('\r', round(i*100L/nrg, 2L), '%');
        }
        
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
        if(length(idx)==0)
        {
            err[i]<-'NVR=0' # number of variants is zero for range i
            next;
        }
        gmx<-gno$gmx[idx,, drop=F]; # row->variant, col->individual
        
        ## save population mean genotype value for later imputation
        avg<-apply(gmx, 1L, mean, na.rm=T); 
        
        ## exclude unavailable individual
        gmx<-gmx[,gdx, drop=F];  # row->variant, col->individual
        
        ## exclude degenerated variants
        idx<-hlp$gmxClr(gmx);
        if(length(idx)==0L)
        {
            err[i]<-'NES=0' ## number of effect sample is zero
            next;
        }
        gmx<-gmx[idx,,drop=F];
        
        ## asign population mean genotype value to NA
        for(j in 1L:nrow(gmx))
            gmx[j, is.na(gmx[j,])]<-avg[j];
        
        ## call U sta
        gmx<-t(gmx); # row->individual, col->genomic variant
        emx<-t(emx); # row->individual, col->RNA probe
        r<-ANL$one(gmx, emx, rsp, cov);
        if (inherits(r, "try-error"))
        {
            err[i] <- errMsg(p);
            next;
        }
        
        ## record result.
        out[i,] <- c(r$GT.u, r$GT.p, r$G.u, r$G.p, r$T.u, r$T.p);
    }
    out<-data.table(
        rng, 
        GT=out[,2L], 
        G=out[,4L],
        T=out[,6L],
        ERR=err, key='CHR,BP1,BP2');
    out;
}
