source('src/helper.R')
require(data.table);


if(!exists('rpt.DBP', inherits = F))
{
	rpt.DBP<-read.table(file = 'out/DBP.AGE.SEX.MED.SMK.rpt', header = T)
	rpt.DBP<-data.table(rpt.DBP, key = 'CHR,BP1,BP2');
	rpt.DBP[, MNP:=min(HWP,NWP,PWP, na.rm = T), by=.(CHR,BP1,GEN,PRB)]->rpt.DBP;
	rpt.DBP[HWP>1, HWP:=1];
}
if(!exists('rpt.SBP', inherits = F))
{
	rpt.SBP<-read.table(file = 'out/SBP.AGE.SEX.MED.SMK.rpt', header = T)
	rpt.SBP<-data.table(rpt.SBP, key = 'CHR,BP1,BP2');
	rpt.SBP[, MNP:=min(HWP,NWP,PWP, na.rm = T), by=.(CHR,BP1,GEN,PRB)]->rpt.SBP;
	rpt.SBP[HWP>1, HWP:=1];
}
if(!exists('rpt.HTN', inherits = F))
{
	rpt.HTN<-read.table(file = 'out/HTN.AGE.SEX.MED.SMK.rpt', header = T)
	rpt.HTN<-data.table(rpt.HTN, key = 'CHR,BP1,BP2');
	rpt.HTN[, MNP:=min(HWP,NWP,PWP, na.rm = T), by=.(CHR,BP1,GEN,PRB)]->rpt.HTN;
	rpt.HTN[HWP>1, HWP:=1];
}
if(!exists('rpt.T15', inherits = F))
{
	rpt.T15<-read.table(file = 'rpt/T15.txt', header = T);
	gvn<-sapply(1L:nrow(rpt.T15), FUN = function(i)
	{
		rpt<-rpt.T15[i,];
		hlp$getGvn(gno, rpt$CHR, rpt$BP1, rpt$BP2, 5000L);
	});
	rpt.T15=data.table(rpt.T15, GVN=gvn);
}

qunf=seq(from=0,to=1,length=nrow(rng)); # quantiles of uniform distrubution on [0, 1]

png(paste('rpt/', 'ALL_G+T','.png', sep=''), width=3, height=1, units="in",res=900, pointsize=4);
pardefault <- par(no.readonly = T) # save plot settings

par(mfrow=c(1, 3));

qqplot(y=rpt.SBP$HWP, x=qunf, xlab="Quntile(U(0,1))", main="SBP", ylab="Quntile(p-value)");
qqline(y=rpt.SBP$HWP, distribution=function(p) qunif(p,0,1), prob=c(0,1), col=2);

qqplot(y=rpt.DBP$HWP, x=qunf, xlab="Quntile(U(0,1))", main="DBP", ylab="");
qqline(y=rpt.DBP$HWP, distribution=function(p) qunif(p,0,1), prob=c(0,1), col=2);

qqplot(y=rpt.HTN$HWP, x=qunf, xlab="Quntile(U(0,1))", main="HTN", ylab="");
qqline(y=rpt.HTN$HWP, distribution=function(p) qunif(p,0,1), prob=c(0,1), col=2);

dev.off();
pardefault <- par(pardefault)    # restore plot settings
