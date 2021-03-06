#include stuff.
source("hwu.core.R");
source("func.R");

# read expressions (e-seq phenotype)
# no missing value handler for no missing expressions.
exp<-read.table(file= paste("exp/", cfg.prb, sep=""), header=FALSE, stringsAsFactors=FALSE);
idv.exp<-read.table(file="exp/idv", header=FALSE, stringsAsFactors=FALSE)[,1]; # individual index in expression files.
row.names(exp)<-idv.exp;

# read covariants
cov<-read.table(file= "phe.csv", header=TRUE, row.names=1, sep=",", stringsAsFactors=FALSE)[cfg.cov];
cov<-miss.rmv(cov, function(x){return(x<0)})

idv.gno<-read.table(file="gno/idv", header=FALSE, stringsAsFactors=FALSE)[,1]; # individual index in genotype files
idv.cov<-row.names(cov);


# find common individual among genotype, covariant and expression.
idv.ins<-intersect(x=idv.gno, y= idv.exp);
idv.ins<-intersect(x=idv.ins, y= idv.cov);
idv.ins<-as.data.frame(idv.ins);

#align expression and covariants
C<-merge(x=cov, y=idv.ins, by.x=0, by.y=1); #Covariant
C<-as.matrix(C[,2:length(C)]);
E<-merge(x=exp, y=idv.ins, by.x=0, by.y=1); #Expression
E<-as.matrix(E[,2:length(E)]);

rng<-read.table(file="rng.txt", header=TRUE, stringsAsFactors=FALSE);
rng.skp<-data.frame(matrix(ncol = 3, nrow=0), stringsAsFactors=FALSE);
rst.nwu<-data.frame(matrix(ncol = 3, nrow=0), stringsAsFactors=FALSE); # non-heterogeneity test
colnames(rng.skp)<-c("seq", "prb", "dsc");

for (i in 1:nrow(rng))
{
	seq_i<-as.character(rng$seq[i]);
	prb_i<-as.character(rng$prb[i]);
	
	# basic file existence check.
	gno_whr<-paste("gno/", seq_i, sep="");
	if(!file.exists(gno_whr))
	{
		rng.skp<-rbind(rng.skp, data.frame(seq=seq_i, prb=prb_i, dsc= "GNO not found."));
		next;
	}
	
	# read i th. genotype and expression
	capture<-tryCatch(
	{
		gno_i<-read.table(file=gno_whr, header=FALSE); # fetch the i th. genotype
	},
	warning = function(war)	{return(war);},
	error = function(err)	{return(err);}, 
	finally = {}); # END tryCatch
	if(inherits(capture, "error") || inherits(capture, "warning"))
	{
		rng.skp<-rbind(rng.skp, data.frame(seq=seq_i, prb=prb_i, dsc=conditionMessage(capture)));
		next;
	}
	row.names(gno_i)<-idv.gno;               # assign genotype row names (individual id)
	gno_i<-miss.avg(gno_i,function(x){x>2}); # imputate missing genotype values
	
	
	print(paste("processing range: ", seq_i, prb_i));
	# align genotype to common idividual ID list. 
	{
		G<-merge(x=gno_i, y=idv.ins, by.x=0, by.y=1); #Genotype
		G<-as.matrix(G[,2:length(G)]); # drop first column - individual IDs
	}
	
	# invoke U core
	capture<-tryCatch(
	{
		G_1<-collapse.burden(G); # join genotype into one dimensional variable
		nwu<-NHWU(Tn=E[,1], geno=G_1, X=C);
	},
	warning=function(war) {return (war);},
	error=function(err)   {return (err);},
	finally =	{});
	if(inherits(capture, "error") || inherits(capture, "warning"))
	{
		rng.skp<-rbind(rng.skp, data.frame(seq=seq_i, prb=prb_i, dsc=conditionMessage(capture)));
		next;
	}
	
	# record U core output
	{
		rst.nwu<-rbind(rst.nwu, data.frame(seq=seq_i, U=signif(nwu$U,4), P=signif(nwu$p,4)));
		#write.table(x=signif(nwu$coef,4),file=paste("coef_nwu/",seq_i,sep=""),quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE);
	}
	
	#clean up i loop
	{
		rm(capture);
		rm(G);
		rm(nwu);
		rm(G_1);
	}
}
#clean up globle
rm(prb_i);
rm(gno_i);
rm(seq_i);
rm(gno_whr);
rm(i);
rm(idv.cov);
rm(idv.gno);
rm(idv.exp);
rm(C);
rm(E);

write.table(x=rst.nwu, file="rst.nwu", quote=FALSE, sep='\t', row.names=FALSE);
write.table(x=rng.skp, file="rng.skp", quote=FALSE, sep='\t', row.names=FALSE);
