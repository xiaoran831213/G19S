#include stuff.
require(data.table);
source('src/analysis.R');

# data file root name.
rut<-"dat/";

if(!exists('phe', inherits = F))
{
	phe<-read.table(file= 'dat/phe.csv', header=TRUE, sep=',', as.is=T);
	ped<-read.table(file= 'dat/ped.csv', header=TRUE, sep=',', as.is=T);
	
	# numerical part of individual id
	iid<-as.integer(sub('T2DG', '', phe$t2dgID));
	
	# re-arrange phenotype to pick only the necessary and sex from pedigree table
	phe<-data.frame(
		IID=iid,          # individual id
		AGE=phe$AGE_1,    # covariant: age
		SEX=ped$SEX,      # covariant: sex
		MED=phe$BPMEDS_1, # covariant: medication
		SMK=phe$SMOKE_1,  # covariatn: smoking
		HTN=phe$HTN_1,    # phenotype: hypertension
		SBP=phe$SBP_1,    # phenotype: systolic blood pressure
		DBP=phe$DBP_1,    # phenotype: distolic blood pressure
		row.names = NULL, stringsAsFactors = F);
	rm(ped, iid);
}

# read genotype
if(!exists('gno', inherits = F))
{
	if(file.exists('dat/gno.bin'))
	{
		load('dat/gno.bin');
	}
	else
	{
		# read genotype individual list
		idv<-read.table('dat/gno.fam', col.names = c("FID", "IID", "PID", "MID", "SEX", "PHE") , as.is = T)
		idv<-as.integer(sub('T2DG', '', idv$IID));
		
		# read genotype map & MAF
		# the MAF
		cls<-c('NULL', 'NULL', 'NULL', 'NULL', 'double', 'NULL');
		frq<-read.table(file = 'dat/gno.frq', header = T, colClasses = cls);
		# the MAP
		cls<-c("integer", "character", "NULL", "integer", "character", "character");
		map<-read.table(file = 'dat/gno.bim', header = F, colClasses = cls);
		map<-data.table(
			IDX=1L:nrow(map),
			CHR=map[,1L], POS=map[,3L],	SNP=map[,2L], A1=map[,4L], A2=map[,5L], MAF=frq$MAF,
			key='CHR,POS');
		rm(frq, cls);
		
		# read genotype matrix
		gmx<-sprintf("tail -n +2 %s | cut -f7- -d ' '", 'dat/gno.raw');
		tmp<-pipe(gmx, "r");
		gmx<-scan(file = tmp, what = integer(0L));
		close(tmp);
		rm(tmp);
		dim(gmx)<-c(nrow(map), length(idv));
		
		# enlist & cleanup
		gno<-list(gmx=gmx, map=map, idv=idv);
		rm(gmx, map, idv);
		save(gno,file = 'dat/gno.bin'); # create binary image
	}
	gc();
}

# read gene expression
if(!exists('exp', inherits = F))
{
	# scan expression file for individual ID
	idv<-sprintf('cut -f1 %s -d \',\'|tail -n+2', 'dat/exp.csv');
	tmp<-pipe(idv, 'r');
	idv<-scan(tmp, what="");
	close(tmp);
	rm(tmp);
	idv<-as.integer(sub('T2DG', '', idv));
	
	# scan for probe ID
	prb<-'dat/exp.csv';
	prb<- scan(file=prb, what="", sep = ',', nlines = 1L);
	prb<- prb[-1L];
	
	# scan for expression matrix
	emx<-sprintf('tail -n+2 %s|cut -f2- -d \',\'', 'dat/exp.csv');
	tmp<-pipe(emx, 'r');
	emx<-scan(tmp, what=double(0), sep=',');
	close(tmp);
	rm(tmp);
	dim(emx)<-c(length(prb), length(idv));	
	
	# packup & cleanup
	exp=list(emx=emx, prb=prb, idv=idv);
	rm(emx, prb, idv);
}

# read gene ranges
if(!exists('rng', inherits = F))
{
	# scan gene-probe-position table
	rng<- 'dat/exp.rng';
	rng<- read.table(rng, header=T, as.is = T, na.strings = 'NULL')
	tmp<- grep('^[0-9]*$', rng$chr);
	chr<- as.integer(rng$chr[tmp]);
	rng<- data.table(
		key='CHR,BP1,BP2',
		CHR=chr, BP1=rng$bp1[tmp], BP2=rng$bp2[tmp],
		GEN=rng$gen[tmp], PRB=rng$prb[tmp]);
	rm(tmp, chr);
}

# prepare file surfix and covariant
cov<-c('AGE', 'SEX', 'MED','SMK');
sfx<-paste(cov, collapse = '.')

if(!exists('tst', inherits = F))
{
	tst<-list('HTN', 'SBP', 'DBP');
	lsp<-list();
}

while(length(tst)>0)
{
	rsp<-tst[[1]];
	lsp[[rsp]]<-system.time(out<-run_std(gno, exp, rng, phe, rsp, cov));
	whr<-sprintf('dat/%s.%s.%s', rsp, sfx,'out');
	write.table(out, whr, row.names = F, quote = F, sep = '\t');
	
	rpt<-data.frame(
		CHR=out$CHR, BP1=out$BP1, BP2=out$BP2, GEN=out$GEN, PRB=out$PRB,
		HWP=format(out$HWP, digits = 3L, scientific = T, trim=T),
		NWP=format(out$NWP, digits = 3L, scientific = T, trim=T),
		PWP=format(out$PWP, digits = 3L, scientific = T, trim=T),
		row.names = NULL, stringsAsFactors = F);
	whr<-sprintf('dat/%s.%s.%s', rsp, sfx, 'rpt');
	write.table(rpt, whr, row.names = F, quote = F, sep = '\t');
	
	rm(whr, rsp, out, rpt);
	
	tst[[1]]<-NULL;
	gc();
}

if(exists('rut', inherits = F)) rm(rut);
if(exists('sex', inherits = F)) rm(sex);
if(exists('cov', inherits = F)) rm(cov);
if(exists('sfx', inherits = F)) rm(sfx);
