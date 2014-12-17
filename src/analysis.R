source('src/hwu.core.R')
source('src/helper.R')
# gmx --- genomic matrix, row<-idv, col<-var
# emx --- expression
# rsp --- response variable
# cvr --- covariants, row<-idv, col<-var
run_inr<-function(gmx, emx, rsp, cvr)
{
	g <- hwu$collapse.burden(gmx); # join genotype into one dimensional variable
	
	k<-matrix(1, nrow(gmx), nrow(gmx));
	n<-hwu$HWU(Tn=rsp, geno=g, X=cvr, K=k, gsim = 'add', appx = 'davis');

	k <- hwu$weight.gaussian(emx); # get weight matrix k
	h<-hwu$HWU(Tn=rsp, geno=g, X=cvr, K=k, gsim = 'add', appx = 'davis');
	
	k <- k-mean(k);
	p<-hwu$HWU(Tn=rsp, geno=g, X=cvr, K=k, gsim = 'add', appx = 'davis');

	list(HWU=list(U=h$U, P=h$p), NWU=list(U=n$U, P=n$p), PWU=list(U=p$U, P=p$p));
}

run_std<-function(gno, exp, rng, phe, rsp, cov)
{
	phe<-phe[, c("IID", rsp, cov)];
	phe<-hlp$clrMss(phe);
	
	tmp<-algIdv(gno, exp, phe);
	gdx<-tmp$gdx;
	edx<-tmp$edx;
	pdx<-tmp$pdx;
	rm(tmp);
	
	# integrity check
	stopifnot(identical(phe$IID[pdx], gno$idv[gdx]));
	stopifnot(identical(phe$IID[pdx], exp$idv[edx]));
	
	# response and covariate
	rsp<-as.matrix(phe[pdx, rsp], rownames.force = F);   # response variable
	cov<-as.matrix(phe[pdx, cov], rownames.force = F);   # covariants
	
	#prepare output
	nrg<-nrow(rng);
	out<-matrix(data = NA, nrow = nrg, ncol = 6L);
	err<-rep.int(NA, nrg);
	# iterate genome variants
	for(i in 1L:nrg)
	{
		if(bitwAnd(i, 0x00FF) == 1L)
		{
			cat('\r', round(i*100L/nrg, 2L), '%');
		}
		
		# get expression matrix
		prb<-rng[i,PRB];
		emx<-exp$emx[exp$prb==prb, edx, drop=F]; # row->probe, col->idvidual
		if(length(emx)==0L)
		{
			err[i] <- 'EXP=N'; # no expression data
			next;
		}
				
		# get testing range i and pick variants within
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
		
		# save population mean genotype value for later imputation
		avg<-apply(gmx, 1L, mean, na.rm=T); 
		
		# exclude unavailable individual
		gmx<-gmx[,gdx, drop=F];  # row->variant, col->individual
				
		# exclude degenerated variants
		idx<-hlp$gmxClr(gmx);
		if(length(idx)==0L)
		{
			err[i]<-'NES=0' # number of effect sample is zero
			next;
		}
		gmx<-gmx[idx,,drop=F];
		
		# asign population mean genotype value to NA
		for(j in 1L:nrow(gmx))
			gmx[j, is.na(gmx[j,])]<-avg[j];
		
		# call U sta
		gmx<-t(gmx); # row->individual, col->genomic variant
		emx<-t(emx); # row->individual, col->RNA probe
		r<-run_inr(gmx, emx, rsp, cov);
		if (inherits(r, "try-error"))
		{
			err[i] <- errMsg(p);
			next;
		}
		
		# record result.
		out[i,] <- c(r$HWU$U, r$HWU$P, r$NWU$U, r$NWU$P, r$PWU$U, r$PWU$P);
	}
	out<-data.table(
		rng, 
 		HWU=out[,1L], 
		HWP=out[,2L],
 		NWU=out[,3L],
		NWP=out[,4L],
 		PWU=out[,5L],
		PWP=out[,6L],
		ERR=err, key='CHR,BP1,BP2');
	out;
}