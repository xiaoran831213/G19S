#include stuff.
source("hwu.core.R");
source("func.R");

phe<-read.table(file= "phe.csv", header=TRUE, row.names=1, sep=",", stringsAsFactors=FALSE)[cfg.phe];
phe<-miss.rmv(phe, function(x){return(x<0)})
cov<-read.table(file= "phe.csv", header=TRUE, row.names=1, sep=",", stringsAsFactors=FALSE)[cfg.cov];
cov<-miss.rmv(cov, function(x){return(x<0)})

idv.gno<-read.table(file="gno/idv", header=FALSE, stringsAsFactors=FALSE)[,1]; # individual index in genotype files
idv.phe<-row.names(phe);
idv.cov<-row.names(cov);

# find common individual among phenotype, genotype, covariant and expression.
idv.ins<-idv.gno;
idv.ins<-intersect(x=idv.ins, y= idv.phe);
idv.ins<-intersect(x=idv.ins, y= idv.cov);
idv.ins<-as.data.frame(idv.ins);

#align phenotype and covariants
P<-merge(x=phe, y=idv.ins, by.x=0, by.y=1); #Phenotype
C<-merge(x=cov, y=idv.ins, by.x=0, by.y=1); #Covariant
P<-as.matrix(P[,2:length(P)]);
C<-as.matrix(C[,2:length(C)]);

rng<-read.table(file="rng.txt", header=TRUE, stringsAsFactors=FALSE);
rng.skp<-data.frame(matrix(ncol = 3, nrow=0), stringsAsFactors=FALSE);
colnames(rng.skp)<-c("seq", "gen", "dsc");
rst.nwu<-data.frame(matrix(ncol = 3, nrow=0), stringsAsFactors=FALSE); # non-heterogeneity test

for (i in 1:nrow(rng))
{
	seq_i<-as.character(rng$seq[i]);
	gen_i<-as.character(rng$gen[i]);
	# basic file existence check.
	gno_whr<-paste("gno/", seq_i, sep="");
	if(!file.exists(gno_whr))
	{
		rng.skp<-rbind(rng.skp, data.frame(seq=seq_i, gen=gen_i, dsc= "GNO not found."));
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
		rng.skp<-rbind(rng.skp, data.frame(seq=seq_i, gen=gen_i, dsc=conditionMessage(capture)));
		next;
	}
	row.names(gno_i)<-idv.gno;               # assign genotype row names (individual id)
	gno_i<-miss.avg(gno_i,function(x){x>2}); # imputate missing genotype values
	
	print(paste("processing range: ", seq_i, gen_i));
	
	# align genotype by common idividual ID. 
	{
		G<-merge(x=gno_i, y=idv.ins, by.x=0, by.y=1); #Genotype
		G<-as.matrix(G[,2:length(G)]); # drop first column - individual IDs
	}
	
	# invoke U core
	capture<-tryCatch(
	{
		G_1<-collapse.burden(G); # join genotype into one dimensional variable
		nwu<-NHWU(Tn=P[,1], geno=G_1, X=C);
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
	rst.nwu<-rbind(rst.nwu, data.frame(seq=seq_i, U=signif(nwu$U,4), P=signif(nwu$p,4)));
}
#clean up globle
rm(i);
rm(capture);
rm(nwu);
rm(G_1);
rm(gno_i);
rm(seq_i);
rm(gen_i);
rm(gno_whr);
rm(idv.phe);
rm(idv.cov);
rm(idv.gno);
rm(G);
rm(C);
rm(P);

write.table(x=rst.nwu, file="rst.nwu", quote=FALSE, sep='\t', row.names=FALSE);
write.table(x=rng.skp, file="rng.skp", quote=FALSE, sep='\t', row.names=FALSE);
