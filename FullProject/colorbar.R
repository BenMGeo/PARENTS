
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='',at=NA,labels=letters[1:nticks],side="right",turn=F,cexlab=1) {
  scale = (length(lut)-1)/(max-min)
  if(!turn){
    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    if (side=="left") {a.side=2} else {if (side=="right"){a.side=4} else {return("wrong combination of side and turn")}}
    axis(a.side, ticks, labels=F,las=0)
    if(all(is.na(at))){
      axis(a.side, ticks, labels=labels,tick=F,las=1,cex=cexlab)
    }else{
      axis(a.side, at, labels=labels,tick=F,las=1,cex=cexlab)
    }
    
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }  
  }else{
    plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    if (side=="top") {a.side=3} else {if (side=="bottom"){a.side=1} else {return("wrong combination of side and turn")}}
    axis(a.side, ticks, labels=F,las=0)
    if(all(is.na(at))){
      axis(a.side, ticks, labels=labels,tick=F,las=1,cex=cexlab)
    }else{
      axis(a.side, at, labels=labels,tick=F,las=1,cex=cexlab)
    }
    
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(y,0,y+1/scale,10, col=lut[i], border=NA)
    }  
  }
}

color.bar.comp <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='',at=NA,labels=letters[1:nticks],side="right",turn=F,srt=0,adj=0,shiftlab=1,cexlab=1) {
  scale = (length(lut)-1)/(max-min)
  if(!turn){
    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    if (side=="left") {a.side=2} else {if (side=="right"){a.side=4} else {print("wrong combination of side and turn");print("side should be right/left or turn should be T");return(NULL)}}
    axis(a.side, ticks, labels=F,las=0)  
    if (a.side==2){
      if(all(is.na(at))){
        axis(a.side, ticks, labels=F,tick=F,las=1)
        text(y=ticks, x =par("usr")[1]-shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }else{
        axis(a.side, at, labels=F,tick=F,las=1)
        text(y=at, x =par("usr")[1]-shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }
    }else{
      if(all(is.na(at))){
        axis(a.side, ticks, labels=F,tick=F,las=1)
        text(y=ticks, x = par("usr")[2]+shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }else{
        axis(a.side, at, labels=F,tick=F,las=1)  
        text(y=at, x = par("usr")[2]+shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }
    }
    
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }  
  }else{
    plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    if (side=="top") {a.side=3} else {if (side=="bottom"){a.side=1} else {print("wrong combination of side and turn");print("side should be top/bottom or turn should be F");return(NULL)}}
    axis(a.side, ticks, labels=F,las=0)
    
    if (a.side==1){
      if(all(is.na(at))){
        axis(a.side, ticks, labels=F,tick=F,las=1)
        text(x=ticks, y =par("usr")[3]-shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }else{
        axis(a.side, at, labels=F,tick=F,las=1)
        text(x=at, y =par("usr")[3]-shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }
    }else{
      if(all(is.na(at))){
        axis(a.side, ticks, labels=F,tick=F,las=1)
        text(x=ticks, y = par("usr")[4]+shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }else{
        axis(a.side, at, labels=F,tick=F,las=1)  
        text(x=at, y = par("usr")[4]+shiftlab, srt = srt, adj = adj, labels = labels, xpd = TRUE,cex=cexlab)
      }
    }
    
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(y,0,y+1/scale,10, col=lut[i], border=NA)
    }  
  }
}

mat2col=function(M,colmap){return(colmap[M])}

colmat2ras=function(mat,ras){ras[!is.na(mat)]=mat[!is.na(mat)];return(ras)}