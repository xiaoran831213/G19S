
PCA <- new.env();

## standardize gnotype matrix
## gmx ---- genotype matrix, row variant, column individual
PCA$std <- function(gmx, maf, rst=F)
{
    whr <- 'dat/pca.std.bin';
    if(rst==T | !file.exists(whr))
    {
        idx <- which(maf>0.05);
        gmx <- gmx[idx,];
        maf <- maf[idx];
        out <- matrix(double(), nrow(gmx), ncol(gmx));
        for(i in 1L:nrow(gmx))
        {
            g <- gmx[i,];
            p <- maf[i];
            g <- (g-mean(g, na.rm=T))/sqrt(p*(1-p))
            g[is.na(g)] <- 0;
            out[i,] <- g;
        }
        save(out, file=whr);
        gc();
    }
    else
    {
        load(file=whr);
    }
    out;
}

PCA$cov <- function(gmx, rst=F)
{
    whr='dat/pca.cov.bin';
    if(rst==T | !file.exists(whr))
    {
        out <- cov(gmx);
        save(out, file=whr);
    }
    else
    {
        load(file=whr);
    }
    out;
}

PCA$run <- function(gno, rst=F)
{
    gmx <- PCA$std(gno$gmx, gno$map$MAF, rst);
    cov <- PCA$cov(gmx, rst);
    egn <- eigen(cov, symmetric=T);
    list(val=egn$values, vec=egn$vectors, cov=cov);
}
