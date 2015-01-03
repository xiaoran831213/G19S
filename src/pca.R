source('src/helper.R')
PCA <- new.env();

## standardize gnotype matrix
## gmx ---- genotype matrix, row variant, column individual
PCA$std <- function(gno, min_maf=0.05, min_dst=1000L, rst=F)
{
    bin <- 'dat/pca.std.bin';

    ## check cached file
    if(rst==F & file.exists(bin))
    {
        load(file=bin);
        return(out);
    }

    ## get genotype map, constraint on MAF
    idx <- which(gno$map$MAF>=min_maf);

    ## constraint on interval
    pos <- gno$map$POS[idx];
    tmp <- pos[1L];
    for(i in 2L:length(idx))
    {
        d <- pos[i]-tmp;
        if(d >= min_dst | d < 0L)
        {
            tmp <- pos[i];
        }
        else
        {
            idx[i] <- 0L; # mark removal
        }
    }
    idx <- idx[idx>0L];
    
    ## constraint on degeneration
    idx <- HLP$gmxClr(gno$gmx, idx);

    ## extract and standardize genotypes
    map <- gno$map[idx,];
    dnm <- sqrt(map$MAF*(1-map$MAF));
    
    gmx <- matrix(double(), nrow=length(idx), ncol=ncol(gno$gmx));
    for(i in 1L:length(idx))
    {
        g <- gno$gmx[idx[i],];
        g <- (g-mean(g, na.rm=T))/dnm[i];
        g[is.na(g)] <- 0;
        gmx[i,] <- g;
    }

    ## repack, build cache, return.
    out <- list(gmx=gmx, map=map, idv=gno$idv);
    save(out, file=bin);
    return(out);
}

## adjust genotype against LD
PCA$adj <- function(gmx, rst=F, max_LD=1)
{
    ## check cached result
    bin <- 'dat/pca.adj.bin';
    if(rst==F & file.exists(bin))
    {
        load(file=bin);
        return (out);
    }

    n <- nrow(gmx);
    m <- ncol(gmx);

    ## perform simple regression, treat the i.th variant as x and the i+1.th as y,
    ## see how much one variant can determine its next.
    ## Y=A+BX+R
    ## since geotype is standardized, A is 0, the model is simplifys to
    ## Y=BX+R
    ## the residual R is supposedly the portion of genotype free of LD, that is,
    ## the portion which cannot be determined by its predecessor.
    ## The 2nd to the last variant will be replaced by residual.
    y <- gmx[1L,];
    rsd <- matrix(double(), nrow=n, ncol=ncol(gmx))
    bta <- double(n);
    r2s <- double(n);
    rsd[1L,] <- y;
    bta[1L] <- NA;
    r2s[1L] <- NA;

    for(i in 2L:n)
    {
        x <- y;
        y <- gmx[i,];

        ## beta = sample_mean(xy)/sample_mean(x2)
        b <- sum(x*y)/sum(x*x);
        bta[i] <- b;
        
        ## r = y - beta*x
        r <- y-b*x;
        rsd[i,] <- r;
        
        ## sum of square residule
        r2s[i] <- sum(r*r);
    }
    r2s <- r2s/(m*m);

    ## exclude high LD
    msk <- which(abs(bta) <= max_LD);
    if(length(msk>0))
    {
        bta <- bta[msk];
        r2s <- r2s[msk];
        rsd <- rsd[msk,];
    }

    ## build cache, return
    out <- list(bta=bta, adj=rsd, r2s=r2s);
    save(out, file=bin);
    out;
}

PCA$cov <- function(dat, rst=F)
{
    ## check cached binery
    bin='dat/pca.cov.bin';
    if(rst==F & file.exists(bin))
    {
        load(file=bin);
        return(out);
    }

    ## Calculate covariance between variables -- the m individuals.
    ## Here genotype variants are treated as n observations.
    out <- cov(dat);

    ## build cache, then return.
    save(out, file=bin);
    out;
}

PCA$run <- function(gno, rst=F)
{
    cat('standardize variants\n');
    t <- system.time(gmx <- PCA$std(gno, rst=rst)$gmx);
    print(t);

    cat('adjust for LD\n');
    t <- system.time(gmx <- PCA$adj(gmx, rst=rst)$adj);
    print(t);
    
    cat('get between individual covariance\n');
    t <- system.time(cov <- PCA$cov(gmx, rst=rst));
    print(t);

    cat('solve principle component\n');
    system.time(egn <- eigen(cov, symmetric=T));
    print(t);

    list(val=egn$values, vec=egn$vectors, cov=cov);
}

PCA$clr <- function()
{
    lsf <- dir('dat', 'pca.*.bin');
    if(length(lsf)<1L)
    {
        ret <- 'cache is already clean';
    }
    else
    {
        lsf <- paste('rm dat/', lsf, sep='');
        ret <- sapply(lsf, system);
        cbind(lsf, ret);
    }
    ret
}
