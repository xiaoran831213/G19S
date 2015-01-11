HLP<-new.env();

## remove individuals with any missing phenotype
HLP$clrMss<-function(phe, fun = function(x) x<0 | is.na(x))
{
    msk<-rep(x=F, times=nrow(phe));
    for(i in 1 : ncol(phe))
        msk = msk | fun(phe[,i]);
    phe<-phe[!msk,];
}

## align gnotype, expression and phenotype by individual ID
HLP$algIdv<-function(gno, exp, phe)
{
    # collect common individual id
    iid <- Reduce(intersect, list(gno$idv, exp$idv, phe$IID));
    iid <- sort(iid);
    
    # index common individual in phenotype
    pdx <- match(iid, phe$IID);
    
    # index common individual in genotype
    gdx <- match(iid, gno$idv);
    
    # index common individual in expression data
    edx <- match(iid, exp$idv);
    
    # return
    list(gdx=gdx, edx=edx, pdx=pdx);
}

## non-physically remove degenerated variants by index
HLP$gmxClr<-function(gmx, idx=1L:nrow(gmx))
{
    ## number of variants
    n <- length(idx);
    
    ## number of individuals
    m <- ncol(gmx);

    ## get variations of variants
    nv <- rep.int(0L,n);
    for(i in 1L:n)
    {
        g <- gmx[idx[i],];
        n0<-sum(g == 0L, na.rm = T);     # 1.Homo-major
        n1<-sum(g == 1L, na.rm = T);     # 2.Hete
        n2<-sum(g == 2L, na.rm = T);     # 3.Homo-minor
        nn<-m-sum(n0, n1, n2);           # 4.missing
        nv[i] <- m-nn-max(n0,n1,n2);     # variation > 0?
    }
    
    ## exclude degenerated variants
    idx[nv>0L];
}

## get error message from try-error object
HLP$errMsg<-function(try_error)
{
    sub('\\(converted from warning\\) ', '', attr(try_error, 'condition')[['message']])
}

HLP$getGvr<-function(gno, chr=0L, bp1=0L, bp2=0L, wnd=0L)
{
    b1<-b1-1L-wnd;
    b2<-b2+1L+wnd;

    idx<-gno$map[CHR==chr][bp1<POS&POS<bp2, IDX];
    list(gmx=gno$gmx[idx,, drop=F], map=gno$map[idx,, drop=F], idv=gno$idv);
}

HLP$getGvr <- function(gno, rng, wnd=0L)
{
    map <- gno$map;
    gmx <- gno$gmx;
    lsg <- list();
    for(i in 1L:nrow(rng))
    {
        chr <- rng[i, CHR];
        bp1 <- rng[i, BP1] -1 - wnd;
        bp2 <- rng[i, BP2] +1 + wnd;
        idx <- map[CHR==chr][bp1<POS & POS<bp2, IDX];
        tag <- paste('G', rng[i, GEN], sep='.');
        gen <- list();
        gen$gmx <- gmx[idx,,drop=F];
        gen$map <- map[idx,,drop=F];
        gen$idv <- gno$idv;
        lsg[[tag]] <- gen;
    }
    lsg$rng <- rng;
    lsg$wnd <- wnd;
    lsg;
}

HLP$getGvn<-function(gno, chr, bp1, bp2, wnd=0L)
{
    bp1<-bp1-1L-wnd;
    bp2<-bp2+1L+wnd;
    gvn<-gno$map[CHR==chr][bp1<POS&POS<bp2, .N];
    gvn;
}

## show progress on screen
HLP$shwPrg <- function(n, i, f=1L)
{
    if(i==n)
    {
        cat('\r, 100%\n');
    }
    else
    {
        if(f==1L | i%%f ==1L)
            cat('\r', round(i*100L/n, 2L), '%');            
    }
}

## save genes to files according to given range list and genotype
HLP$gsv <- function(gno, rng, wnd=5000L, rut='dat/gen')
{
    require(data.table);

    ## number of gene regions
    n <- nrow(rng);
    
    for(i in 1L:n)
    {
        r <- as.list(rng[i,]);
        
        ## get genome range for gene region i.
        gmx <- NULL;
        map <- gno$map[CHR==r$CHR][r$BP1-1L-wnd<POS & POS<r$BP2+1L+wnd,];
        if(nrow(map)>0L)
            gmx <- gno$gmx[map$IDX, ,drop=F];
        
        ## get unique record name and file name
        whr <- sprintf('%s/G%04X', rut, r$SEQ);
        
        ## pack and save
        gen <- list(
            seq=r$SEQ, chr=r$CHR, bp1=r$BP1, bp2=r$BP2, wnd=wnd,
            gen=r$GEN, prb=r$PRB,
            gmx=gmx, map=map, idv=gno$idv);
        save(gen, file=whr);

        ## report progress
        HLP$shwPrg(n, i, 100L);
    }
}

## load genes from filees according to given rang list, they must
## be first extracted from genotype data and saved.
HLP$gld <- function(rng, rut='dat/gen', rdc=T)
{
    require(data.table);

    out <- new.env()
    for(i in 1L:nrow(rng))
    {
        seq <- sprintf('G%04X', rng[i, SEQ]);
        load(sprintf('%s/%s', rut, seq));

        tag <- rng[i, GEN];
        key <- sprintf('%s.%s', seq, tag);
        out[[key]] <- gen;
    }
    if(length(out)==1L & rdc)
        out <- gen;
    out;
}

HLP$mktab <- function(lst)
{
    ## number of row, number of columns
    n <- length(lst);                   # row

    ## decide column count
    m <- max(sapply(lst, length));      # col

    ## column heads
    h <- Reduce(union, lapply(lst, names));
    
    ## column types
    t <- as.list(rep(NA, times=m));
    names(t) <- h;

    ## decide column classes
    for(l in lst)
    {
        i <- which(is.na(t));
        if(length(i)<1L)
            break;
        j <- which(!is.na(l));
        if(length(j)<1L)
            next;
        x <- intersect(names(t[i]), names(l[j]));
        if(length(x)<1L)
            next;
        t[x] <- sapply(l[x], class);
    }
    t <- as.character(t);

    ## temporary chrarcter matrix and empty data.frame.
    d <- matrix(NA, n, m, dimnames=list(NULL, h));
    f <- as.data.frame(d);

    ## fill temporary matrix
    for(i in 1L:n)
    {
        l <- lst[[i]];
        if(length(l) > 0L)
        {
            d[i,] <- as.character(l[h]);
        }
    }

    ## fill data.frame, complete type cast.
    for(i in 1L:m)
    {
        f[,i] <- as(d[,i], t[i]);
    }
    f;
}
