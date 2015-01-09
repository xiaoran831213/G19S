require(data.table);
source('src/anl.R');

CTL <- new.env();

## read phenotype
CTL$phe <- function()
{
    ## phenotypes and covariants
    phe<-read.table(file= 'dat/phe.csv', header=TRUE, sep=',', as.is=T);

    ## only need sex from pedigree file
    ped<-read.table(file= 'dat/ped.csv', header=TRUE, sep=',', as.is=T);

    ## integrity check
    stopifnot(identical(phe$t2dgID, ped$ID));
    
    ## numerical part of individual id
    iid<-as.integer(sub('T2DG', '', phe$t2dgID));
    
    ## re-arrange phenotype to pick only the necessary and sex from pedigree table
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
}

## read genotype
CTL$gno <- function(rst=F)
{
    bin <- 'dat/gno.bin'
    if(file.exists(bin) & !rst)
    {
        load(bin);
    }
    else
    {
        ## read genotype individual list
        col <- c("FID", "IID", "PID", "MID", "SEX", "PHE");
        idv<-read.table('dat/gno.fam', col.names = col, as.is = T)
        idv<-as.integer(sub('T2DG', '', idv$IID));
        
        ## read genotype map & MAF
        ## the MAF
        cls<-c('NULL', 'NULL', 'NULL', 'NULL', 'double', 'NULL');
        frq<-read.table(file = 'dat/gno.frq', header = T, colClasses = cls);
        ## the MAP
        cls<-c("integer", "character", "NULL", "integer", "character", "character");
        map<-read.table(file = 'dat/gno.bim', header = F, colClasses = cls);
        map<-data.table(
            IDX=1L:nrow(map),
            CHR=map[,1L], POS=map[,3L],
            SNP=map[,2L], A1=map[,4L], A2=map[,5L], MAF=frq$MAF,
            key='CHR,POS');
        
        ## read genotype matrix
        gmx<-sprintf("tail -n +2 %s | cut -f7- -d ' '", 'dat/gno.raw');
        tmp<-pipe(gmx, "r");
        gmx<-scan(file = tmp, what = integer(0L));
        close(tmp);
        dim(gmx)<-c(nrow(map), length(idv));
        
        ## enlist & cleanup
        out<-list(gmx=gmx, map=map, idv=idv);
        save(out, file = bin); # create binary image
    }
    out;
}

## read gene outression
CTL$exp <- function(rst=F)
{
    bin <- 'dat/exp.bin';
    whr <- 'dat/exp.csv';
    if(file.exists(bin) & rst==F)
    {
        load(bin);
    }
    else
    {
        ## scan expression file for individual ID
        idv<-sprintf('cut -f1 %s -d \',\'|tail -n+2', whr);
        tmp<-pipe(idv, 'r');
        idv<-scan(tmp, what="");
        close(tmp);
        idv<-as.integer(sub('T2DG', '', idv));
        
        ## scan for probe ID
        prb<- scan(file=whr, what="", sep = ',', nlines = 1L);
        prb<- prb[-1L];
        
        ## scan for outression matrix
        emx<-sprintf('tail -n+2 %s|cut -f2- -d \',\'', whr);
        tmp<-pipe(emx, 'r');
        emx<-scan(tmp, what=double(0), sep=',');
        close(tmp);
        dim(emx)<-c(length(prb), length(idv));    
        
        ## packup and save binery image
        out=list(emx=emx, prb=prb, idv=idv);
        save(out, file=bin);
    }
    out;
}

## read gene ranges
CTL$rng <- function(rst=F)
{
    bin <- 'dat/rng.bin';
    whr <- 'dat/rng.txt';
    if(file.exists(bin) & rst==F)
    {
        load(bin);
    }
    else
    {
        ## scan gene-probe-position table
        out<- read.table(whr, header=T, as.is = T, na.strings = 'NULL')

        ## only take numbered chromosomes, exclude sex and mitochondria
        tmp<- grep('^[0-9]*$', out$chr);
        chr<- as.integer(out$chr[tmp]);
        out<- data.table(
            key='CHR,BP1,BP2',
            CHR=chr, BP1=out$bp1[tmp], BP2=out$bp2[tmp],
            GEN=out$gen[tmp], PRB=out$prb[tmp]);

        ## insert sequence number as the first column.
        out <- data.table(SEQ=1L:nrow(out), out);
        setkey(out, CHR, BP1, BP2);
        save(out, file=bin);
    }
    out;
}

CTL$pcs <- function()
{
    pcs <- read.table(file='dat/pcs.ldv', header=T, sep='\t', as.is=T);
    pcs$IID <- NULL;
    as.matrix(pcs);
}

CTL$run <- function(
    pck=NULL,
    pcs=NULL,
    cvr=c('AGE', 'SEX', 'MED','SMK'),
    rsp=c('DBP', 'SBP', 'HTN'))
{
    bin=format.Date(Sys.Date(), 'dat/run.%y-%m-%d %H:%M:%S');
    
    ## prepare file surfix and covariant
    cat('fatching materials\n');
    gno <- CTL$gno();
    exp <- CTL$exp();
    rng <- CTL$rng();
    phe <- CTL$phe();

    ## rng <- rng[1L:100L,]; #debug perposep
    if(!is.null(pck))
        rng <- rng[pck,];
    
    tmr <- list();
    out <- list();
    for(y in rsp)
    {
        cat('\ranalysing phenotype: ', y, '\n');
        t <- system.time(o<-ANL$run(gno, exp, rng, phe, y, cvr, pcs));
        out[[y]] <- o;
        tmr[[y]] <- t;
    }
    out[['PCS']] <- pcs;
    out[['TMR']] <- tmr;
    save(out, file=bin);

    print(tmr);
    out;
}

CTL$rpt <- function(out, dst=NULL, order=F)
{
    if(order)
        out <- CTL$odr(out, lowest=F, simple=F);

    out <- format(out, trim=T, digits=3L, scientific=T);
    if(!is.null(dst))
    {
        write.table(out, dst, quote=F, sep='\t');
        out <- dst;
    }
    out;
}

CTL$odr <- function(out, lowest=T, simple=T)
{
    odr <- out[order(MIN),];
    if(lowest)
    {
        odr <- odr[odr[,list(IDX=min(.I)), by=GEN][,IDX],];
    }
    if(simple)
    {
        odr <- odr[, list(GEN, GT, G, T, MIN)];
    }
    odr;
}

CTL$qqp <- function(OUT, dst=NULL)
{
    ## quantiles of uniform distrubution on [0, 1]
    qunf=seq(from=0,to=1,length=nrow(OUT$SBP));

    ## create file
    if(!is.null(dst))
    {
        png(dst, width=3, height=1, units="in",res=900, pointsize=4);
    }
    
    ## save plot settings
    pardefault <- par(no.readonly = T)

    par(mfrow=c(1, 3));

    qqplot(y=OUT$SBP$GT, x=qunf, xlab="U(0,1)", main="SBP", ylab="p-value");
    qqline(y=OUT$SBP$GT, distribution=function(p) qunif(p,0,1), prob=c(0,1), col=2);

    qqplot(y=OUT$DBP$GT, x=qunf, xlab="U(0,1)", main="DBP", ylab="");
    qqline(y=OUT$DBP$GT, distribution=function(p) qunif(p,0,1), prob=c(0,1), col=2);

    qqplot(y=OUT$HTN$GT, x=qunf, xlab="U(0,1)", main="HTN", ylab="");
    qqline(y=OUT$HTN$GT, distribution=function(p) qunif(p,0,1), prob=c(0,1), col=2);

    ## restore plot settings
    pardefault <- par(pardefault)

    ## release file
    if(!is.null(dst))
    {
        dev.off();
    }
}

CTL$top <- function(OUT, n=5, dst=0L, lwp=T)
{
    sbp <- CTL$odr(OUT$SBP, lowest=lwp, simple=F)[1:n,];
    sbp <- cbind(PHE='SBP', sbp);

    dbp <- CTL$odr(OUT$DBP, lowest=lwp, simple=F)[1:n,];
    dbp <- cbind(PHE='DBP', dbp);

    htn <- CTL$odr(OUT$HTN, lowest=lwp, simple=F)[1:n,];
    htn <- cbind(PHE='HTN', htn);

    t15 <- rbind(sbp, dbp, htn);
    t15$ERR <- NULL;

    if(is.character(dst))
    {
        t15 <- format(t15, trim=T, digits=3, scientific=T)
        write.table(t15, dst, quote=F, sep='\t', row.names=F);
        t15 <- dst;
    }
    if(dst==0L)
    {
        t15 <- t15[,list(PHE,GEN,GT,G,T,MIN)];
    }
    t15;
}
