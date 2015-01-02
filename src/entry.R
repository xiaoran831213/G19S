require(data.table);
source('src/analysis.R');

## read phenotype
getPhe <- function()
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
getGno <- function(rst=F)
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
getExp <- function(rst=F)
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
getRng <- function(rst=F)
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
        save(out, file=bin);
    }
    out;
}

start <- function()
{
    ## prepare file surfix and covariant
    cvr<-c('AGE', 'SEX', 'MED','SMK');
    sfx<-paste(cvr, collapse = '.')

    lsp<-list();
    for(rsp in c('HTN', 'SBP', 'DBP'))
    {
        lsp[[rsp]] <- system.time(out<-run_std(gno, exp, rng, phe, rsp, cvr));
        whr<-sprintf('dat/%s.%s.%s', rsp, sfx, 'out');
        write.table(out, whr, row.names = F, quote = F, sep = '\t');
        
        rpt<-data.frame(
            CHR=out$CHR, BP1=out$BP1, BP2=out$BP2, GEN=out$GEN, PRB=out$PRB,
            HWP=format(out$HWP, digits = 3L, scientific = T, trim=T),
            NWP=format(out$NWP, digits = 3L, scientific = T, trim=T),
            PWP=format(out$PWP, digits = 3L, scientific = T, trim=T),
            row.names = NULL, stringsAsFactors = F);
        whr<-sprintf('dat/%s.%s.%s', rsp, sfx, 'rpt');
        write.table(rpt, whr, row.names = F, quote = F, sep = '\t');
        gc();
    }
}
