miss.avg<-function(T, f=is.na)
{
	for(i in 1:ncol(T))
	{
		avg=mean(T[,i][!f(T[,i])]);
		T[,i][f(T[,i])]<-avg;
	}
	return(T);
}
miss.rmv<-function(tab, f=is.na)
{
	flag<-rep(FALSE,nrow(tab));
	for(i in 1:ncol(tab))
	{
		flag[f(tab[,i])]<-TRUE;
	}
	rs<-tab[!flag,];
	rs<-as.data.frame(rs);
	row.names(rs)<-row.names(tab)[!flag];
	colnames(rs)<-colnames(tab);
	return(rs);
}
join.inner<-function(t1,t2)		#inner join two table by row.names
{
	rs=merge(x=t1,y=t2,by.x=0,by.y=0);
	row.names(rs)<-rs[[1]];
	rs<-rs[2:length(rs)];
	return(rs);
}
