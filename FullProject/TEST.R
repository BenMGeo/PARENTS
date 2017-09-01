source('R:/iP.PCArotator/FullProject/Main_Reduction.R')
require(plyr)
require(raster)
require(rgdal)
require(fBasics)

#Test Main_Reduction.R

P1=matrix(c(rep(10,9),10,rep(8,7),10,10,8,rep(5,5),8,10,10,8,5,rep(2,3),5,8,10,10,8,5,2,0,2,5,8,10,10,8,5,rep(2,3),5,8,10,10,8,rep(5,5),8,10,10,rep(8,7),10,rep(10,9)),9,9)
P2=matrix(9:1,9,9)
P3=t(P2)
P4=matrix(c(1,0),9,9)

Patts=data.frame(as.vector(P1),as.vector(P2),as.vector(P3))

Patts=predict(prcomp(Patts))

P1=matrix(Patts[,1],9,9)
P2=matrix(Patts[,2],9,9)
P3=matrix(Patts[,3],9,9)

M1=runif(1,0,3)*P1+runif(1,1,3)*P2+runif(1,0,2)*P3+runif(9*9,-1,1)
M1=(M1-min(M1))/(max(M1)-min(M1))*8
M2=runif(1,0,3)*P1+runif(1,1,3)*P2+runif(1,0,2)*P3+runif(9*9,-1,1)
M2=(M2-min(M2))/(max(M2)-min(M2))*8
M3=runif(1,0,3)*P1+runif(1,1,3)*P2+runif(1,0,2)*P3+runif(9*9,-1,1)
M3=(M3-min(M3))/(max(M3)-min(M3))*8
M4=runif(1,0,3)*P1+runif(1,1,3)*P2+runif(1,0,2)*P3+runif(9*9,-1,1)
M4=(M4-min(M4))/(max(M4)-min(M4))*8
M5=runif(1,0,3)*P1+runif(1,1,3)*P2+runif(1,0,2)*P3+runif(9*9,-1,1)
M5=(M5-min(M5))/(max(M5)-min(M5))*8
M6=runif(1,0,3)*P1+runif(1,1,3)*P2+runif(1,0,2)*P3+runif(9*9,-1,1)
M6=(M6-min(M6))/(max(M6)-min(M6))*8

MS=stack(raster(M1),raster(M2),raster(M3),raster(M4),raster(M5),raster(M6))
names(MS)=LETTERS[1:6]

Patts=stack(raster(P1),raster(P2),raster(P4))

projection(MS)=CRS("+init=epsg:32632")
projection(Patts)=CRS("+init=epsg:32632")

names(Patts)=paste("P",1:3,sep="_")

plot(MS)
plot(Patts)

Patts=unstack(Patts)
names(Patts)=paste("P",1:3,sep="_")

RES=SignalReduction(MS,Patts[1:3],MaxEval=1000,printint=250,minPlot=0.001,maxFrac=0.35,minFrac=0.05,numFrac=11)

plot(MS)

plot(stack(RES$RedInput),maxnl=3*6)

plot(RES$RotOpt$P_3)

plot(RES$RotOpt$P_2)

plot(RES$RotOpt$P_1)

RES$Corrs

#####

CCI=stack("R:/Work/CCI_SM_20151018.nc",varname="sm")
LST=stack("R:/Work/lst_cut.nc")
LST[which(values(LST)==max(values(LST),na.rm=T))]=NA
LAI=stack("R:/Work/lai_cut.nc")
LAI[which(values(LAI)==max(values(LAI),na.rm=T))]=NA



#####

Rstack=stack(file.choose())

Rstack[][Rstack[]==0]=NA

names(Rstack)=paste("Pic",LETTERS[1:nlayers(Rstack)],sep="_")

DEM=raster(file.choose())

DEM.der=terrain(DEM,opt=c('slope', 'aspect', 'TPI', 'TRI', 'roughness', 'flowdir'))

HS=hillShade(DEM.der$slope,DEM.der$aspect,angle=70,direction=180)

P=list(HS=HS,DEM=DEM)

RES=SignalReduction(Rstack,P)

plot(stack(RES$RedInput),maxnl=30)

plot(RES$RotOpt$DEM)

plot(RES$RotOpt$HS)

RES$Corrs

