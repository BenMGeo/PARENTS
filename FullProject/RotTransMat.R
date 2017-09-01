RotTransMat=function(alpha){
  
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
        R[[paste(i,j,sep="->")]]=Rrot(i,j,alphaM[i,j],bl+1)
      }
    }
    
    # Step 3: multiply the list elements
    Rfull=Reduce('%*%',R)
    
    return(Rfull)
    
  }
}