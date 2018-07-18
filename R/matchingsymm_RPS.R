#' This function obtains the individual resistant-symmetric shape for 2D
#' matching-symmetry data. The input is an array A of size
#' n (landmarks) x p (dimensions) x 2k (objects: the left-right sides for each)
#' Configurations are ordered in this way: left side Object 1, right side Object 1,
#' left side Object 2, right side Object 2, etc
#'
#' @param A  an array of size n (landmarks) x 2 (in 2D) x 2k (left/right sides for k configurations)
#' @param ctr Centering options: "gmedian" (the spatial or gemetric median, default choice), "median" (the componentwise median), "mean" (the average)
#' @param legend.loc The location of the legend for the plot.result function
#'
#' @usage
#' matchingsymm_RPS(A,ctr="gmedian",legend.loc="topleft")
#'
#' @author Federico Lotto, Sebastian Torcida
#'
#' @export
matchingsymm_RPS<-function(A,ctr="gmedian",legend.loc="topleft"){

  nl<-length(A[,1,1])#numero de landmarks
  n<-length(A[,1,1]) # The number of landmarks
  k<-dim(A)[3]/2     # The number of configurations (objects)

  # Initial centering----------------
  Mc<-center(A,cent=ctr)

  for(i in 1:k){
    #Side arrays are created----------
    Mcl<-Mc[,,(2*i-1),drop=FALSE]
    Mcr<-Mc[,,(2*i),drop=FALSE]
  }

  dv<-Mcl-Mcr  # Difference vectors between corresponding left-right sided lanmarks

  mod<-as.matrix(sqrt((dv[,1,]^2)+(dv[,2,]^2)))#Norm of diff. vectors
  N<-aperm(array(mod,c(dim(dv)[1],n,2)),c(1,3,2))#Creates a copy
  ndv<-dv/N #Normalises diff. vectors
  ndvp<-aperm(apply(ndv,c(1,3),function(x)if(x[2]<0) {x<-x*-1}else{x<-x}),c(2,1,3))#Makes vectors point in x>0 direction
  par.mat<-as.matrix(acos(ndvp[,1,]))#The angle with x axis
  #Leaves angles in the range [0,pi]

  dr<-apply(par.mat,2,median) # The median angle of the diff. vectors
  vr<-cbind(cos(dr),sin(dr))  # Computes the spherical median direction

  # Computes the Househölder reflection matrix using the spherical median as mirror
  Ur<-array(0,dim(Mcl))
  for(i in 1:k){
    R<-diag(2) - 2*vr[i,]%*%t(vr[i,])
    Ur[,,i] <- Mcl[,,i]%*%R  # Performs the corresponding reflection
  }

  S<-array(t(c(Mcr,Ur)),c(nl,2,(n*2)))
  print("please wait...",quote = F)
  capture.output(U<- RPS::robgit_RPS(S))
  #A resistant Procrustes superposition to filter out remaining differences
  #in size and/or orientation between sides
  Ud<-U[,,1:n,drop=FALSE]
  Ui<-U[,,(n+1):(n*2),drop=FALSE]

  TT<-(Ud+Ui)/2 #Computes the symmetric shape

  # Computes each landmark contribution percentage (%) to total asymmetry
  distances<-dist.contrib(mc=Ui,mre=Ud)

  # Plots
  plot.result(mc=Ud, mre =Ui, mt=TT, nconf = k, object=T, legloc = legend.loc)
  return(list(distances,Ui-Ud))
}

center<-function(A,cent){
    ni<-dim(A)[3]
    switch(cent,
         mean={Mc<-apply(A,c(2,3), function(x) x-mean(x))},
         median={Mc<-apply(A,c(2,3), function(x) x-median(x))},
         gmedian={Mc<-array(0,c(dim(A)))
         gmed<-t(apply(A,c(3), function(x) Weiszfeld(x)))
         for(i in 1:ni){
           Mc[,,i]<-t(t(A[,,i])-c(gmed[[i]][[1]]))
         }
         }
  )
  return(Mc)
}


dist.contrib<-function(mc=NULL,mre=NULL){

  eudist<-sqrt(((mc[,1,]-mre[,1,])^2)+((mc[,2,]-mre[,2,])^2))

  suma<-colSums(eudist) # Computes total shape difference (asymmetry)
  perc<-round((eudist/suma)*100,digits=4) # Landmark contribution percentage

  results<-array(0,c(dim(mc)))
  results[,1,]<-round(eudist,digits=6)
  results[,2,]<-perc

  print("Distance and contribution % to total asymmetry for each configuration:",quote = F)
  print(results)
  return(results)
}


plot.result<-function(mc=NULL,mre=NULL,mt=NULL,nconf,object=TRUE,legloc="topleft"){
  conf<-1:nconf
  if(object[1]){
    for(c in conf){
      plot(mc[,,c],asp=1,xlab = "",ylab = "",main = paste("Config",c)) # Side 1
      points(0,0,pch=3) # Center in (0,0)
      points(mre[,,c],pch=20) # Reflected
      points((mt[,,c]+mre[,,c]),col="red",pch=20,cex=0.8) # Symmetric (residual + reflected)
      text(mc[,1,c],mc[,2,c],c(seq(1:(dim(mc)[1]))),pos=2,offset=0.5,cex=0.55) # Landmark number
      legend(x=legloc,c("Side 1","Reflected Side 2","Symmetric"),pch=c(1,20,20) ,col = c("black","black","red"), cex=0.7)
    }

  }else{

    for(c in conf){
      plot(mc[,,c],asp=1,xlab = "",ylab = "",main = paste("Config",c))#original
      points(0,0,pch=3)#center in (0,0)
      points(mre[,,c],col="red",pch=20,cex=0.8)#reflected
      text(mc[,1,c],mc[,2,c],c(seq(1:dim(mc)[1])),pos=2,offset=0.5,cex=0.55)
      legend(x=legloc,c("Left","Right"),pch=c(20,1) ,col = c("black","red"), cex=0.7)
    }
  }

}

