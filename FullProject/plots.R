overview.figure=function(RES, col.map, col.scat, labels=c()){
  
  normMat=function(Mat,min,max){Mat=(Mat-min)/(max-min);return(Mat)}
  
  num=length(RES$Corrs)
  
  m0=matrix(c(1,2,2,3,3,4,4),ncol=1)
  m=m0
  count=1
  
  while (ncol(m)<num){
    m=cbind(m,m0+4*count)
    count=count+1
  }
         
  nf=layout(m)
  par(oma=c(0.5,4,1.5,1)+0.1)
  
  for (i in 1:num){
    
    name=names(RES$Corrs)[i]
    
    par(mar=c(0.5,1,8,1))
    
    ran=range(RES$RotOpt[[name]][[name]][],na.rm=T)
    
    color.bar.comp(lut=col.map,floor(ran[1]*10)/10,ceiling(ran[2]*10)/10,nticks=4, at=round(seq(floor(ran[1]*10)/10,ceiling(ran[2]*10)/10,length.out = 4),2),labels=round(seq(floor(ran[1]*10)/10,ceiling(ran[2]*10)/10,length.out = 4),2),shiftlab=10, adj=0.5,turn=T,side="top",cexlab=1.4)
    mtext(labels[i],line=4,cex=1.5)
    par(mar=c(0,0,1,0))
    
    DefaultR=as.raster(matrix(NA,ncol=dim(RES$RotOpt[[name]])[[2]],nrow=dim(RES$RotOpt[[name]])[[1]]))
    
    plotpat=lapply(lapply(lapply(lapply(unstack(RES$RotOpt[[name]]),as.matrix),normMat,floor(ran[1]*10)/10,ceiling(ran[2]*10)/10),"*",length(col.map)-1),"+",1)
    # 
    plotpat=lapply(lapply(plotpat,mat2col,col.map),colmat2ras,DefaultR)
    
    if (i==1){
      plot(0,type="n",xlab="",ylab="",xlim=c(xmin(RES$RotOpt[[name]]),xmax(RES$RotOpt[[name]])),ylim=c(ymin(RES$RotOpt[[name]]),ymax(RES$RotOpt[[name]])),asp=1.5,xaxt="n",cex.axis=1.3)
    }else{
      plot(0,type="n",xlab="",ylab="",xlim=c(xmin(RES$RotOpt[[name]]),xmax(RES$RotOpt[[name]])),ylim=c(ymin(RES$RotOpt[[name]]),ymax(RES$RotOpt[[name]])),asp=1.5,xaxt="n",yaxt="n")
    }
    rasterImage(plotpat[[2]],xmin(RES$RotOpt[[name]]),ymin(RES$RotOpt[[name]]),xmax(RES$RotOpt[[name]]),ymax(RES$RotOpt[[name]]),asp=1.5)
    if (i==1){
      mtext("original",side=2, line=2.5,cex=1.3)
    }
    
    par(mar=c(1,0,0,0))
    if (i==1){
      plot(0,type="n",xlab="",ylab="",xlim=c(xmin(RES$RotOpt[[name]]),xmax(RES$RotOpt[[name]])),ylim=c(ymin(RES$RotOpt[[name]]),ymax(RES$RotOpt[[name]])),asp=1.5,cex.axis=1.3)
    }else{
      plot(0,type="n",xlab="",ylab="",xlim=c(xmin(RES$RotOpt[[name]]),xmax(RES$RotOpt[[name]])),ylim=c(ymin(RES$RotOpt[[name]]),ymax(RES$RotOpt[[name]])),asp=1.5,yaxt="n",cex.axis=1.3)
    }
      rasterImage(plotpat[[1]],xmin(RES$RotOpt[[name]]),ymin(RES$RotOpt[[name]]),xmax(RES$RotOpt[[name]]),ymax(RES$RotOpt[[name]]),asp=1.5)
    if (i==1){
      mtext("rotated",side=2, line=2.5,cex=1.3)
    }
    
    par(mar=c(3.5,3.5,2,2))
    smoothScatter(RES$RotOpt[[name]][[name]][], RES$RotOpt[[name]][[paste(name,"_mod",sep="")]][], asp=1, colramp=col.scat, xlab="", ylab="",xaxt="n",yaxt="n")
    axis(1, mgp=c(3, .6, 0),cex.axis=1.2)
    axis(2, mgp=c(3, .6, 0),cex.axis=1.2)
    abline(0,1)
    abline(lm(RES$RotOpt[[name]][[paste(name,"_mod",sep="")]][]~RES$RotOpt[[name]][[name]][]),lty=2,col="red",lwd=2)
    usr <- par( "usr" )
    mtext(paste("r =",round(RES$Corrs[[name]],2)),side=1,line = -2.0,adj = 0.9,cex=1.2)
    mtext("original", side=1,line = 2.0,cex=1.3)
    mtext("rotated", side=2,line = 2.0,cex=1.3)
    
  }
  
}

plotpos = "R:\\iP.PCArotator\\"

png(paste(plotpos,"fig.png",sep=""),width=length(RES$Corrs)*7.5/cm(1),height=4*6/cm(1),units="in",res=300)

col.map=colorRampPalette(c("blue", "red","yellow"))(101)
col.scat=colorRampPalette(c("white", "blue"))
labels=c(
  "height a.s.l.  [ m ]",
  "hill shade [ - ]",
  bquote(atop("thermal inertia","["~10^-3~cal~cm^-2~degree*C~s^0.5~"]")),
  expression(paste("LAI  [ ",m^2," ",m^-2," ]", sep = "")),
  "uniform noise  [ - ]"
)
overview.figure(RES, col.map, col.scat,labels)

dev.off()