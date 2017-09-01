
myDDS=function(X,LB,UB,FUN,...,pars=list(R=0.2,MaxEval=1000,printint=100,plotter=F)){
  
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
  printint= pars$printint;           # interval for printing actual optimum
  plotter= pars$plotter;             # if results shall be plotted
  
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
      print(length(fopts))
      print(range(((1-R)*10-1)*MaxEval+(1:length(fopts))))
      points(((1-R)*10-1)*MaxEval+(1:length(fopts)),fopts,pch=19,cex=0.5,col=colorRampPalette(c("green","red"))(9)[R*10])}
    }
    
    # Step 6: Stop Criteria
    if(count==MaxEval){
      quit=1
      print("MaxEval reached!")
    }  
  }
  return(list(Xopt=Xopt,Fopt=Fopt))
}