hlp<-new.env();

#remove individuals with any missing phenotype
hlp$clrMss<-function(phe, fun = function(x) x<0 | is.na(x))
{
	msk<-rep(x=F, times=nrow(phe));
	for(i in 1 : ncol(phe))
		msk = msk | fun(phe[,i]);
	phe<-phe[!msk,];
}

# align gnotype, expression and phenotype by individual ID
hlp$algIdv<-function(gno, exp, phe)
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

# non-physically remove degenerated variants by index
hlp$gmxClr<-function(gmx)
{
	# calculate genotype statistics is necessary
	ngv<-nrow(gmx); # number of genome variants
	ndv<-ncol(gmx); # number of individuals
	# get value statistics for all variants
	col<-c("N0", "N1", "N2", "NN", "NV")
	stt<-matrix(data = 0L, nrow = ngv, ncol = 5L, dimnames = list(NULL, col));
	for(i in 1L:ngv)	
	{
		n0<-sum(gmx[i,] == 0L, na.rm = T);      # 1.Homo-major
		n1<-sum(gmx[i,] == 1L, na.rm = T);      # 2.Hete
		n2<-sum(gmx[i,] == 2L, na.rm = T);      # 3.Homo-minor
		nn<-ndv-sum(n0, n1, n2);     # 4.missing
		nv<-ndv-nn-max(n0, n1, n2);  # 5.variation
		stt[i,]<-c(n0, n1, n2, nn, nv);
	}
	
	# mask degenerated variants.
	msk<-which(stt[,'NV']>0); # variation count
	msk
}

# get error message from try-error object
hlp$errMsg<-function(try_error)
{
	sub('\\(converted from warning\\) ', '', attr(try_error, 'condition')[['message']])
}

hlp$getGvr<-function(gno, chr, bp1, bp2, wnd=0L)
{
	bp1<-bp1-1L-wnd;
	bp2<-bp2+1L+wnd;
	idx<-gno$map[,IDX:=.I][CHR==chr][bp1<POS&POS<bp2, IDX];
	list(gmx=gno$gmx[idx,, drop=F], map=gno$map[idx,, drop=F], idv=gno$idv);
}

hlp$getGvn<-function(gno, chr, bp1, bp2, wnd=0L)
{
	bp1<-bp1-1L-wnd;
	bp2<-bp2+1L+wnd;
	gvn<-gno$map[CHR==chr][bp1<POS&POS<bp2, .N];
	gvn;
}