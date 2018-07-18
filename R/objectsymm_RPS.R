#' This function obtains the individual resistant-symmetric shape for 2D
#' object-symmetry data. The input is an array A of size
#' n (landmarks) x p (dimensions) x k (objects)
#' Landmarks must be in this order: saggital (or unpaired) landmarks first,
#' then left paired landmarks and finally right paired landmarks
#' Configurations are ordered in this way: L side Object 1 and R side Object 1,
#' L side Object 2 and R side Object 2, etc
#'
#' @param A Input data: an array or matrix of size n (landmarks) x 2 (in 2D) x k (objects)
#' @param ctr Centering options: "gmedian" (the spatial or gemetric median, default choice), "median" (the componentwise median), "mean" (the average)
#' @param prs.file This is a .txt file indicating the L+R paired landmarks as rows: e.g. 7 15; 8 16; etc.
#' @param proj.met The choice to compute the saggital axis: sum or median of projections
#' @param legend.loc The location of the legend for the plot.result function
#'
#' @return w
#'
#' @author Federico Lotto, Sebastian Torcida
#'
#' @usage
#'  objectsymm_RPS(A,ctr="gmedian",prs.file,proj.met="msum",legend.loc="topleft")
#'
#'@importFrom "graphics"  "legend"  "plot" "points" "text"
#'@importFrom "stats" "median" "runif"
#'@importFrom "utils" "capture.output" "combn" "head" "read.table"
#'@importFrom "Gmedian" "Weiszfeld"
#'
#' @export
#'
objectsymm_RPS<-function(A,ctr="gmedian",prs.file,proj.met="msum",legend.loc="topleft"){

  pares<-read.table(prs.file,sep=" ",header=FALSE)
  sag<-length(A[,1,1])-(length(pares[,1])*2) # number of saggital or unpaired landmarks
  ppaired<-length(A[,1,1])-sag # number of paired landmarks
  n<-dim(A)[3]

  #Initial centering----------------
  Ac<-center(A,cent=ctr)  # centered data

  #Saggital, left-paired and right-paired landmarks in that order are stored in different arrays
  m0<-Ac[1:sag,,,drop=FALSE]
  m1<-Ac[(sag+1):((sag+1)+(ppaired/2)-1),,,drop=FALSE]
  m2<-Ac[((sag+1)+(ppaired/2)):((sag)+(ppaired)),,,drop=FALSE]

  #Saggital-------------------------
  pares.sag<-t(combn(1:sag,2)) # combinations of size 2 of saggital landmarks
  vs<-m0[pares.sag[,1],,,drop=FALSE]-m0[pares.sag[,2],,,drop=FALSE]# Saggital difference vectors
  mod<-as.matrix(sqrt((vs[,1,]^2)+(vs[,2,]^2))) # Saggital diff. vectors length
  N<-aperm(array(mod,c(dim(pares.sag)[1],n,2)),c(1,3,2))
  vsn<-vs/N # Saggital difference vectors are made of length 1


  #Paired landmarks-----------------
  vp<-m1[pares[,1]-sag,,,drop=FALSE]-m2[pares[,2]-sag-(ppaired/2),,,drop=FALSE] #Paired difference vectors

  # Saggital differences will be ranked by the sum (or median)
  # of the projections onto them of the paired diff. vectors

  P<-array(0,dim = c(dim(vs)[1],(ppaired/2),n))

  for(i in 1:n){
    P[,,i]<-(vsn[,,i]%*%t(vp[,,i]))
  }

  P<-abs(P) # Projections are stored

  p.sum<-apply(P,c(1,3), sum) # Sum of projections; each column corresponds to an object
  p.median<-apply(P,c(1,3),median) # Median of projections

  # Identifies the saggital diff. vector w/minimum sum (median) of project
  switch(proj.met,
    msum={i.min<-apply(p.sum,c(2),which.min)},
    mmedian={i.min<-apply(p.median,c(2),which.min)}
  )

  vr<-NULL # An auxiliar zero matrix

  #Select those sagg. diff. vectors with minimum sum (median) of projections
  for(j in 1:n){
    vaux<-vsn[i.min[j],,j]
    vr<-rbind(vr,vaux)
  }

  #Reflection----------------
  #Computes the unitary direction tha is orthogonal to the vector achieving minimum sum
  e<-matrix(0,n,2)
  e[,2]<-vr[,1]
  e[,1]<- -(vr[,2])

  #Computes the corresponding Househölder reflection
  Ur<-array(0,dim(A))
  for(i in 1:n){
    R<-diag(2) - 2*e[i,]%*%t(e[i,])
    Ur[,,i] <- Ac[,,i]%*%R  # Performs the reflection
  }

  #Auxiliar relabelling of paired landmarks
  Ure<-Ur
  for (j in 1:(ppaired/2)) {
    Ure[pares[j, 1],,] <- Ur[pares[j,2],,]
    Ure[pares[j, 2],,] <- Ur[pares[j,1],,]
  }

  TT<-(Ac+Ure)/2 #Computes the symmetric shape

  #Computes each landmark contribution percentage (%) to total asymmetry
  distances<-dist.contrib(mc=Ac,mre=Ure)

  #Plots
  plot.result(mc=Ac,mre=Ure,mt=TT,nconf=n,object=T,legloc=legend.loc)

  #Returns the contributions to asymmetry and the symmetric shape
  return(list(distances,(Ac-Ure)/2))

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

  suma<-colSums(eudist)#Computes total shape difference (asymmetry)
  perc<-round((eudist/suma)*100,digits=4)#landmark contribution percentage

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
      plot(mc[,,c],asp=1,xlab = "",ylab = "",main = paste("Config",c))#original
      points(0,0,pch=3)#center in (0,0)
      points(mre[,,c],pch=20)#reflected
      points((mt[,,c]+mre[,,c]),col="red",pch=20,cex=0.8)#simmetric (residuals + reflected)
      text(mc[,1,c],mc[,2,c],c(seq(1:(dim(mc)[1]))),pos=2,offset=0.5,cex=0.55)#number of landmark
      legend(x=legloc,c("Original","Reflected","Symmetric"),pch=c(1,20,20) ,col = c("black","black","red"), cex=0.7)
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

