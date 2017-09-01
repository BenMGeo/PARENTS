# princomp vs prcomp

loc.DF2raster=function(DF,Raster){
  
  # DF: data frame of raster layers fitting to non-NA in Raster
  # Raster: projected Raster
  
  blanko=Raster
  
  list.R=list()
  
  for (i in 1:ncol(DF)){
    blanko[][!is.na(blanko[])]=DF[,i]  
    list.R[i]=blanko
  }
  
  Res=stack(list.R)
  names(Res)=colnames(DF)
  
  return(Res)
  
}


D=na.omit(as.data.frame(Rstack))


D=as.data.frame(MS)
matplot(D,type="l")


PR=prcomp(D)

PRIN=princomp(D)


P.PR=predict(PR)

P.PRIN=predict(PRIN)


matplot(P.PR,type="l")

matplot(P.PRIN,type="l")

matplot(abs(P.PRIN)-abs(P.PR),type="l")

varimax(loadings(PRIN)) #works

varimax(loadings(PR)) #doesn't work

HKs=fa(D,nfactors=6,rotate="none",fm="pa",warnings=F)

P.HKs=predict(HKs,D)
matplot(P.HKs,type="l")

P.HKs=loc.DF2raster(P.HKs,mean(Rstack))
P.PR=loc.DF2raster(P.PR,mean(Rstack))
P.PRIN=loc.DF2raster(P.PRIN,mean(Rstack))

target.rot(loadings(HKs))

plot(P.PR)
plot(P.PRIN)
plot(P.HKs)