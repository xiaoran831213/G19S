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

HLP$getGvr<-function(gno, chr, bp1, bp2, wnd=0L)
{
	bp1<-bp1-1L-wnd;
	bp2<-bp2+1L+wnd;
	idx<-gno$map[,IDX:=.I][CHR==chr][bp1<POS&POS<bp2, IDX];
	list(gmx=gno$gmx[idx,, drop=F], map=gno$map[idx,, drop=F], idv=gno$idv);
}

HLP$getGvn<-function(gno, chr, bp1, bp2, wnd=0L)
{
	bp1<-bp1-1L-wnd;
	bp2<-bp2+1L+wnd;
	gvn<-gno$map[CHR==chr][bp1<POS&POS<bp2, .N];
	gvn;
}
