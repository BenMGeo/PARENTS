SignalReduction=function(Input,Patterns,...){
  
  # Input: raster tack of original Data
  # Patterns: list of patterns to be extracted, best named
  # ...: prameters, see below
  
  # RedInput: Input reduced by Patterns
  # RotOpt: list of representations of the modeled patterns
  # Corr: list of correlation hits
  
  ##### Packages #####
  
  require(plyr)
  require(raster)
  require(rgdal)
  require(fBasics)
  
  ##### Parameters #####
  
  ## Defaults
  pca.scale=F # PCA scaling
  plotter=T # plot of progress
  MaxEval=200 # maximum of evaluation processes with optimization
  maxFrac=0.9 # maximum fraction
  minFrac=0.1 # minimum fraction
  numFrac=9 # number of fractions
  printint=50 # print interval
  minPlot=0.01 # minimum for optimization plot
  maxStart=1 # maximum starting point for DDS
  maxSearch=10 # maximum number of searches for starting point 
  
  ## from ...
  
  if( length(list(...)) ){
    # read list
    Lst <- list(...)
    
    if( !is.null(Lst$maxFrac) ){
      maxFrac <- Lst$maxFrac
    }
    
    if( !is.null(Lst$minFrac) ){
      minFrac <- Lst$minFrac
    }
    
    if( !is.null(Lst$numFrac) ){
      numFrac <- Lst$numFrac
    }
    
    if( !is.null(Lst$pca.scale) ){
      pca.scale <- Lst$pca.scale
    }
    
    if( !is.null(Lst$plotter) ){
      plotter <- Lst$plotter
    }
    
    if( !is.null(Lst$MaxEval) ){
      MaxEval <- Lst$MaxEval
    }
    
    if( !is.null(Lst$printint) ){
      printint <- Lst$printint
    }
    
    if( !is.null(Lst$minPlot) ){
      minPlot <- Lst$minPlot
    }
    
    if( !is.null(Lst$maxStart) ){
      maxStart <- Lst$maxStart
    }   
    
    if( !is.null(Lst$maxSearch) ){
      maxSearch <- Lst$maxSearch
    }
  }
  
  ##### used functions #####
  
  loc.myDDS=function(X,LB,UB,FUN,...,pars=list(R=0.2,MaxEval=1000,len=10,pos=1,printint=100,plotter=F)){
    
    # Fun:  Handle to the objective Function, where "..." is additional input for function
    # X:    Innitial guess for the optimum
    # LB:   Lower Boundary of the Parameter Space
    # UB:   Upper Boundary of the Parameter Space
    # Xopt: Empirical optimum
    # Fopt: Value of the Objective Function in the optimum
    # pars: list of DDS parameters
    
    
    # Step 1: DDS parameters
    R       = pars$R;                  # neighborhood pertubation size (0.2 default by authors)
    MaxEval = pars$MaxEval;            # maximum number of function evaluations
    len = pars$len                     # length of colours
    pos = pars$pos                     # which position to plot
    printint = pars$printint;          # interval for printing actual optimum
    plotter = pars$plotter;            # if results shall be plotted
    
    # Step 2: initial evaluation of the function
    count   = 1;
    F       = FUN(X,...);  
    Fopt    = F;
    Xopt    = X;
    quit    = 0;
    
    #optional plot
    if(plotter){fopts=c(Fopt)}
    
    while(quit==0){
      
      # Step 3: Randomly select J of the N variables
      Pi = 1-log(count)/log(MaxEval);
      Pn = runif(length(X));
      J=which(Pn<Pi);
      nj=length(J);
      if (nj==0){
        J=floor(runif(1)*length(X)+1);
      }
      
      # Step 4: calculate updated xnew
      sig    = R*(UB-LB)*rnorm(length(X));
      Xnew   = Xopt;
      Xnew[J]= Xnew[J]+sig[J];
      
      i      = which(Xnew<LB);
      Xnew[i]= LB[i]+(LB[i]-Xnew[i]);
      i1     = which(Xnew[i]>UB[i]);
      Xnew[i1]=LB[i1];
      
      z       = which(Xnew>UB);
      Xnew[z] = UB[z] - (Xnew[z]-UB[z]);
      z1      = which(Xnew[z]<LB[z]);
      Xnew[z1]= UB[z1];
      
      # Step 5: 
      count   = count+1;
      F       = FUN(Xnew,...);  
      if (F<Fopt){
        Fopt    = F;
        Xopt    = Xnew;
      }
      
      #optional
      if(plotter){fopts=c(fopts,Fopt)}
      
      if (count%%printint==0){
        print(paste("step:",count))
        print(paste("Fopt:",Fopt))
        
        #optional plot
        if(plotter){
          points((pos-1)*MaxEval+1:length(fopts),fopts,pch=19,cex=0.5,col=colorRampPalette(c("red","green"))(len)[pos])}
      }
      
      # Step 6: Stop Criteria
      if(count==MaxEval){
        quit=1
        print("MaxEval reached!")
      }  
    }
    return(list(Xopt=Xopt,Fopt=Fopt))
  }
  
  loc.RotTransMat=function(alpha){
    
    # alpha: vector of radiant angles of rotation between n-dimensions in the order of (1<>2, 1<>3, 1<>4, 1<>5,..., 1<>n, 2<>3, 2<>4, 2<>5,..., 2<>n,...,n-1<>n)
    # Rfull: rotation matrix belonging to the angles alpha
    
    Rrot=function(r,c,theta,n){
      
      # Help function for getting the submatrix for the actual rotation between dimension r and c out of n for angle theta 
      
      Rrot=diag(1,n)
      Rrot[r,r]=Rrot[c,c]=cos(theta)
      Rrot[r,c]=-sin(theta)
      Rrot[c,r]=sin(theta)
      return(Rrot)
    }
    
    # Step 1: get angles in matrix form
    alphaM=diag(0,(1+sqrt(1+8*length(alpha)))/2)
    alphaM[upper.tri(alphaM)]=alpha
    
    
    if (ncol(alphaM)!=nrow(alphaM)){
      return("Input is not compatible!") #not expected
    }else{
      # Step 2: get a list of all rotation matrices for the angles 
      bl=ncol(alphaM)
      R=list()
      for (i in 1:(bl-1)){
        for (j in (i+1):bl){
          R[[paste(i,j,sep="->")]]=Rrot(i,j,alphaM[i,j],bl)
        }
      }
      
      # Step 3: multiply the list elements
      Rfull=Reduce('%*%',R)
      
      return(Rfull)
      
    }
  }
  
  loc.PCAstack=function(Stack,scale=T){
    
    # Stack: raster stack for PCA calculation
    # scale: logical, should raster stack be scaled
    
    # CompStack: raster stack of components from PCA
    # s: summary of PCA (sufficient for most problems)
    # full: full information from PCA
    
    # Step 1: get basic data frame from raster stack
    stackframe=na.omit(as.data.frame(Stack))
    
    # Step 2: calculate and predict PCA based on data frame
    stackframe.pca=prcomp(stackframe,scale=scale)
    s=summary(stackframe.pca)
    PCA=predict(stackframe.pca)
    
    # Step 3: bring back to raster stack form
    EmptyLay=mean(Stack)
    CompStack=list()
    
    for (i in 1:dim(stackframe)[2]){
      EmptyLay[]=EmptyLay[]*NA
      PCi=PCA[,i]
      EmptyLay[as.numeric(names(PCi))]=PCi
      CompStack[[i]]=EmptyLay    
    }
    
    CompStack=stack(CompStack)
    
    return(list(PCA=CompStack,summary=s,full=stackframe.pca))
  }
  
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
  
  loc.rotoptfun=function(alpha,PCA.DF,P.AIM){
    
    # alpha: angle parameters to be optimised
    # PCA.DF: component arrays in data frame
    # P.AIM: pattern that the first component is to be aimed for
    
    # Step 1: rotate in actual direction and extract first array
    RMat=loc.RotTransMat(alpha)
    dir1=(as.matrix(PCA.DF)%*%RMat)[,1]
    
    # Step 2a: scale pattern to first array's domain and calculate correlation (minimizable)
    aimmin=loc.normDF(t(t(P.AIM)))
    aimmin=aimmin*diff(range(dir1))+min(dir1) 
    resmin=(1-abs(cor(unlist(aimmin),unlist(dir1))))
#     resmin=resmin * (sqrt(mean((unlist(aimmin)-unlist(dir1))^2)))
    
    # Step 2b: scale inverted pattern to first array's domain and calculate correlation (minimizable)
    aimmax=-(loc.normDF(t(t(P.AIM)))-1)
    aimmax=aimmax*diff(range(dir1))+min(dir1)
    resmax=(1-abs(cor(unlist(aimmax),unlist(dir1))))
#     resmax=resmax * (sqrt(mean((unlist(aimmax)-unlist(dir1))^2)))
    
    # Step 3: return minimum of both possibilities
    return(min(c(resmin,resmax)))
  }
  
  loc.normDF=function(DF,rowwise=F,interval="0-1"){
    
    # DF: data frame to be normalized
    # rowwise: logical, columnwise if false
    # interval: 2 options to choose; "0-1" or "-1-1" for different minima
    
    if (interval=="0-1"){
      normV=function(V){
        V=(V-min(V,na.rm=T))/(max(V,na.rm=T)-min(V,na.rm=T))
        return(V)
      }
    }else{
      if (interval=="-1-1"){
        normV=function(V){
          V=V/max(abs(V),na.rm=T)
          return(V)
        }
      }else{
        print("Wrong Setup for intervals!")
        return(NULL)
      }
    }
    
    if (!rowwise) {return(colwise(normV)(as.data.frame(DF)))}
    else {return(t(colwise(normV)(as.data.frame(t(DF)))))}
  }
  
  loc.ProjCheck=function(R1,R2){
    
    # R1: raster to be changed
    # R2: raster to be checked
    
    if (!(projection(R1)==projection(R2))){
      R1=projectRaster(R1,R2)
    }
    
    if (!(extent(R1)==extent(R2)) || !(res(R1)==res(R2))){
      R1=resample(R1,R2)
    }
    return(R1)
  }
  
  ##### actual function #####
  
  # Step 0: Initialization of Output
  RotOpt=list()
  Corrs=list()
  Outputs=list()
  Angles=list()
  
  # Step 1: prepare original derivated data
  Input.save=Input
  Patterns=lapply(Patterns,loc.ProjCheck,Input)
  tridi=((nlayers(Input))^2-(nlayers(Input)))/2 #number of entries for angles
  
  print("Input prepared, calculation starts!")
  
  # Step 2 (looped through patterns): Calculation
  for (P in names(Patterns)){
    
    # Step 2a: calculate PCA with actual rotation and produce data frame with pattern
    I.PCA=loc.PCAstack(Input,scale=pca.scale)
    I.PCA.c=loc.DF2raster(as.data.frame(as.matrix(na.omit(as.data.frame(Input)))%*%I.PCA$full$rotation),mean(Input))
    FullData=na.omit(as.data.frame(stack(I.PCA.c,Patterns[[P]])))
    

    # Step 2b: initialize angles 
    ncount=0
    cInv.Opt=Inf
    while(cInv.Opt > maxStart && ncount< maxSearch){
      ALPHA.Opt=list(runif(tridi,-pi,pi))
      cInv.Opt=loc.rotoptfun(ALPHA.Opt[[1]],FullData[,1:nlayers(Input)],FullData[,nlayers(Input)+1])
      ncount=ncount+1
    }
    
    # Step 2c: optimization
    if(plotter){
      par(mfrow=c(1,1))
      plot(0,cInv.Opt,
           main=paste("current optimization process for pattern",P),type="p",log="y",
           ylim=c(minPlot,cInv.Opt),xlim=c(0,MaxEval*numFrac),
           ylab="fopt",col="black",pch=19,cex=2)
    }
    
    FracSeq=seq(maxFrac,minFrac,l=numFrac)
    for (R in FracSeq){
      print(paste("Freedom:",round(R,3)))
      ALPHA.Opt=loc.myDDS(ALPHA.Opt[[1]],rep(-pi,tridi),rep(pi,tridi),loc.rotoptfun,FullData[,1:nlayers(Input)],FullData[,nlayers(Input)+1],pars=list(R=R,MaxEval=MaxEval,len=length(FracSeq),pos=which(FracSeq==R),printint=printint,plotter=plotter))
    }
    
    # Step 2d: calculation of rotated results
    RotMat=loc.RotTransMat(ALPHA.Opt[[1]]) #[1:tridi]?
    fullRotAngleU1=acos(RotMat[1,1]/sqrt(sum(RotMat[,1]^2)))*180/pi
    PCA.c.rot=as.data.frame(as.matrix(na.omit(as.data.frame(stack(Input,Patterns[[P]])))[,-nlayers(Input)+1])%*%I.PCA$full$rotation%*%RotMat)
    #rotPCs=PCA.c.rot
    PCA.c.rot=loc.DF2raster(PCA.c.rot,mean(stack(I.PCA$PCA,Patterns[[P]])))
    names(PCA.c.rot)=paste(P,"PC",1:nlayers(PCA.c.rot),sep=".")
    
    # Step 2e: check data for which orientation is best
    Pmod1<-Pmod2<-PCA.c.rot[[1]]
    Pmod1[]=loc.normDF(t(t(Pmod1[])))[,1]*diff(range(Patterns[[P]][],na.rm=T))+min(Patterns[[P]][],na.rm=T)
    Pmod2[]=-(loc.normDF(t(t(Pmod2[])))[,1]-1)*diff(range(Patterns[[P]][],na.rm=T))+min(Patterns[[P]][],na.rm=T)
    
    choose=abs(cor(na.omit(as.data.frame(stack(Patterns[[P]],Pmod1,Pmod2)))))
    choose.P=which.max(choose[1,2:3])
    
    if(choose.P[1]==1){
      Pmod=Pmod1
    }else{
      if(choose.P[1]==2){
        Pmod=Pmod2
      }else{
        print("This should never happen!")
      }
    }
    
    loc.RotOpt=stack(Pmod,Patterns[[P]])
    names(loc.RotOpt)=c(paste(P,"mod",sep="_"),P)
    
    RotOpt[[P]]=loc.RotOpt
    Corrs[[P]]=choose[1,1+choose.P[1]]
    Angles[[P]]=fullRotAngleU1
    
    # Step 2f: calculate reduced Input for next step
    
    Outputs[[P]]=PCA.c.rot
    
#     # Step 2f: calculate reduced Input for next step
#     # !!! Is there a better way???? !!! This seems wrong
#     loads=I.PCA$full$rotation
#     
#     NDF=list()
#     
#     for(i in 1:nlayers(Input)){
#       
#       NDF[[names(Input)[i]]]=0
#       
#       for(j in 2:ncol(rotPCs)){
#         
#         NDF[[names(Input)[i]]]=NDF[[names(Input)[i]]]+rotPCs[,j]*loads[i,j]
#         
#       }
#     }
#     
#     NDF=as.data.frame(NDF)
#     
#     Input=loc.DF2raster(NDF,mean(stack(I.PCA$PCA,Patterns[[P]])))

    
#     # Step 2f: calculate reduced Input for next step 
#     # !!! Is there a better way???? !!! This seems wrong
#     FullDF=na.omit(as.data.frame(stack(Input,stack(unstack(PCA.c.rot)[-1]))))
#     
#     NDF=list()
#     
#     for(i in 1:nlayers(Input)){
#       LM=lm(as.formula(paste(names(Input)[i],"~", paste("V",c(2:nlayers(Input)),sep=""),sep="")),data=FullDF)
#       D=predict(LM,data=FullDF)
#       NDF[[colnames(FullDF)[i]]]=D
#     }
#     NDF=as.data.frame(NDF)
#     
#     Input=loc.DF2raster(NDF,mean(stack(I.PCA$PCA,Patterns[[P]])))
    
    print(paste("Pattern",P,"finished!"))
  }
  
  # Step 3: Setup output and return
  Corrs=as.data.frame(Corrs)
  #RedInput=Input
  RedInput=Outputs
  print("Calculation finished!")
  
  return(list(RedInput=RedInput,RotOpt=RotOpt,Corrs=Corrs,Angles=Angles))
}



