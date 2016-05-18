
#### Sequential Regression Imputation algorithm

seqreg<-function(dat,tol,maxit,facs,counts){
  n1<-length(dat[,1]) # Number of obs
  m<-length(dat[1,]) # Number of var
  it<-0 ## Number of iterations
  facind<-unique(c(raoul.factors(dat),facs)) ## Find factors
  countind<-unique(c(raoul.counts(dat),counts))
  ## Expanding factors to dummies
  orgorder<-colnames(dat) ## Perserving the original order of the variables
  orgrows<-rownames(dat) ## Perserving original row names

  facnames<-colnames(dat[,facind])
  levnames<-NULL ## level names
  levnames2<-NULL
  levnames3<-NULL
  facind2<-NULL
  levlen<-NULL ## level lenght
  levlen2<-NULL
  if(!is.null(facind)){ ## If there are more than zero factors
  names(facind)<-colnames(dat)[facind]
  k<-length(facind)
  ## Extract level names
  for(i in 1:k){
    lnm<-levels(as.factor(dat[,facind[i]]))
    levlen<-cbind(levlen,length(lnm)-1)
    length(lnm)<-length(levnames)<-max(length(lnm),length(levnames))
    levnames<-cbind(levnames,lnm)
    colnames(levnames)[i]<-names(facind)[i]
    if(max(is.na(dat[,facind[i]]))==1){
      levlen2<-cbind(levlen2,length(lnm)-1)
      length(levnames3)<-max(length(lnm),length(levnames3))
      levnames3<-cbind(levnames3,lnm)
    }
    
  }
  levnames2<-as.matrix(levnames[-1,]) ## Drop first level (reference)

  
## Expand factors to dummies
  for(i in 1:k){ ###### Consider making this a function
    dum<-model.matrix(~as.factor(dat[,facind[i]])+0)
    dum2<-matrix(,n1,length(dum[1,]))
    rownames(dum2)<-1:n1
    dum3<-dum[match(rownames(dum2),rownames(dum)),]
    colnames(dum3)<-array(names(facind)[i],length(dum3[1,]))
    rownames(dum3)<-orgrows
    if(i==1){
      dat2<-cbind(dat[,seq(1:m)<facind[i]],dum3[,-1, drop=FALSE])
      facind2<-seq(facind[i],length=levlen[i])
    }else{
      dat2<-cbind(dat2,as.matrix(dat[,(seq(1:m)>facind[(i-1)] & seq(1:m)<facind[i])]),dum3[,-1, drop=FALSE])
      facind2<-c(facind2,seq((facind[i]+levlen[i-1]-1),length=levlen[i]))
    }
    
  }
  dat2<-cbind(dat2,dat[,seq(1:m)>facind[k]])
  }else{
  dat2<-dat
}
  

    ## Updating Names for expanded factors
  namelength<-rep(1,m)
  namelength[facind]<-levlen
  orgorder2<-rep(orgorder,namelength)
  if(!is.null(facs)){
  factorlevnames<-paste(rep(facnames,levlen),levnames2[which(!is.na(levnames2))],sep="")
  }
  orgorder2[facind2]<-factorlevnames
  
  
  ## Taking original names and ordering and number of vars; total, numeric, and factors

  dat2<-cbind(1,dat2) # Add constant
  m1<-length(dat2[1,]) # Number of vars
  m3<-length(dat2[1,facind2]) # Number of (expanded) Factors
  
  # Select fully observed / missing variables; total, numeric, and factors
  mispat<-raoul.mispat(dat2) ## Extract missingness pattern
  misvar<-rbind(1:m1,apply(mispat,2,max)) ## Find variables with missing values
  mobs<-m1-sum(misvar[2,]) ## Number of observed variables(+constant)
  facind3<-facind2+1 ## Compensate for addition of constant
  
  ## Finding count variables in dat2
  countind2<-NULL
  for(i in 1:length(countind)){
    coun<-countind[i]+sum((levlen-1)[which(facind<countind[i])])+1
    countind2<-c(countind2,coun)
  }
  
  
  numind2<-setdiff(c(1:m1),c(facind3,countind2)) ## Find numeric variables
  
  
  misnumer<-as.vector(as.matrix(misvar[,numind2])[1,which(misvar[2,numind2]==1)]) ## Numeric variables with missing
  misfac<-as.vector(as.matrix(misvar[,facind3])[1,which(misvar[2,facind3]==1)]) ## Factors with missing
  misfacnames<-colnames(dat2)[misfac] ## Names of factors
  miscount<-as.vector(as.matrix(misvar[,countind2])[1,which(misvar[2,countind2]==1)]) ## Count variables with missing
  
  dat2<-data.matrix(dat2)
  
  ## Creating vecs for fast imputations
  misnumvec<-vec(dat2[,misnumer])
  impnumvec<-misnumvec
  misfacvec<-vec(dat2[,misfac])
  impfacvec<-misfacvec
  miscountvec<-vec(dat2[,miscount])
  impcountvec<-miscountvec
  
  impdat<-dat2
  olsout<-NULL
  logisticout<-NULL
  countout<-NULL
  
  ## First impuation step
  
  # Imputation of numeric
  if(length(misnumer)>0){
  olsout<-raoul.olsfit(mis=as.matrix(dat2[,misnumer]),full=as.matrix(dat2[,which(misvar[2,]==0)]),m=misnumer,it=it)
  impnumvec[which(is.na(misnumvec)==TRUE),]<-vec(olsout[[1]])[which(is.na(misnumvec)==TRUE),]
  imp1<-matrix(impnumvec,n1,length(misnumer))
  impdat[,misnumer]<-imp1
  }
  
  
  
  # Imputation of factors (multinomial probability)
    if(length(misfac)>0){
    misfact<-as.matrix(dat2[,misfac])
    colnames(misfact)<-colnames(dat2)[misfac]
    logisticout<-raoul.logistic1(misfact=misfact,full=as.matrix(dat2[,which(misvar[2,]==0)]),m=facind3,n=n1,it=it,levlen=levlen,misp=mispat[,misfac],ice=NULL)
    impfacvec[which(is.na(misfacvec)==TRUE),]<-vec(logisticout[[1]])[which(is.na(misfacvec)==TRUE),]
    imp2<-matrix(impfacvec,n1,length(misfac))
    colnames(imp2)<-misfacnames
    impdat[,misfac]<-imp2
  }
  
  
  # Imputation of count variables (poisson)
  if(length(miscount)>0){
    countout<-raoul.countfit(mis=as.matrix(dat2[,miscount]),full=as.matrix(dat2[,which(misvar[2,]==0)]),m=miscount,it=it, misp=mispat[,miscount],ice=NULL)
    impcountvec[which(is.na(miscountvec)==TRUE),]<-vec(countout[[1]])[which(is.na(miscountvec)==TRUE),]
    imp3<-matrix(impcountvec,n1,length(miscount))
    impdat[,miscount]<-imp3
  }
  
  ## Setting up evaluation
  lik1<-raoul.lik2(impdat,n=n1) # First likelihood
  lik2<-(1+2*tol)*lik1 #set comparing likelihood to 2
  it<-1
  
  ## Subsequent imputations
  impdat2<-impdat
  while(abs(lik2/lik1-1)>tol & it<maxit){
    
    lik1<-lik2
    
    ### Subsequent imputation of numeric variables
    if(length(misnumer>0)){
    olsout<-raoul.olsfit(imp1,impdat,misnumer,it)
    impnumvec[which(is.na(misnumvec)==TRUE),]<-vec(olsout[[1]])[which(is.na(misnumvec)==TRUE),]
    imp1<-matrix(impnumvec,n1,length(misnumer))
    impdat2[,misnumer]<-imp1
    }
    
    ### Subsequent imputation of factors
    if(length(misfac)>0){
    misfact<-as.matrix(imp2)
    colnames(misfact)<-colnames(dat2)[misfac]
    logisticout<-raoul.logistic1(misfact=misfact,full=impdat,m=misfac,n1=n1,it=it,levlen=levlen2,misp=as.matrix(mispat[,misfac]),ice=logisticout)
    impfacvec[which(is.na(misfacvec)==TRUE),]<-vec(logisticout[[1]])[which(is.na(misfacvec)==TRUE),]
    imp2<-matrix(impfacvec,n1,length(misfac))
    colnames(imp2)<-misfacnames
    impdat2[,misfac]<-imp2
    }
    
    # Subsequent Imputation of count variables (poisson)
    if(length(miscount)>0){
      countout<-raoul.countfit(mis=as.matrix(imp3),full=impdat,m=miscount,it=it,misp=as.matrix(mispat[,miscount]),ice=countout)
      impcountvec[which(is.na(miscountvec)==TRUE),]<-vec(countout[[1]])[which(is.na(miscountvec)==TRUE),]
      imp3<-matrix(impcountvec,n1,length(miscount))
      impdat2[,miscount]<-imp3
    }
    
    ### Merger and calculation of likelihood
    impdat<-impdat2 # Updata full-data matrix
    lik2<-raoul.lik2(impdat,n=n1) # Calculate new likelihood
    it<-it+1
    print(it)
  }
  
  
  
  
  
  ## New name for matrices without intercept
  ml<-impdat[,-1]
  mispat2<-mispat[,-1]

  
  
  ## Obtain SSE of regressions without missing values
  sse<-as.matrix(cbind(matrix(0,1,mobs-1),raoul.sse(fitted=as.matrix(olsout[[1]]),
                        original=as.matrix(dat2[,misnumer])),logisticout[[3]],countout[[3]]))
  
  
  
  ## Obtaining matrix for identifying type of missingness
  idfy<-as.matrix(cbind(matrix(0,1,mobs-1),matrix(1,1,length(olsout[[1]][1,])),
                        matrix(2,1,length(logisticout[[1]][1,])), matrix(3,1,length(countout[[1]][1,]))))
  idfy[1,setdiff(facind2,misfac-1)]<-"5"
  
  
  
  
    ## Obtaining covariance matrix without missing values
  vs<-as.matrix(var(dat2[,-1],na.rm=TRUE))
  
  
  
  ## Obtaining weight-vector
  w<-as.vector(apply(1-mispat2,1,sum)/length(mispat2[1,]))
  w<-w*length(mispat2[,1])/sum(w)
  
  
  
  ## Obtaining weighted X'X
  xtx<-t(impdat*w)%*%impdat
  
  
  
  ## Ordering the MSE and IDFY matrixes
  colnames(sse)<-colnames(idfy)<-c(orgorder2[-c(misnumer-1,misfac-1,miscount-1)],orgorder2[c(misnumer-1)],orgorder2[c(misfac-1)],orgorder2[c(miscount-1)])
  sse<-sse[,orgorder2]
  sse<-rbind(sse,apply((1-mispat2),2,sum))
  idfy<-idfy[,orgorder2]
  idfy<-rbind(idfy,colnames(dat2[,-1]))
  
  
  
  ## Creating matrix of beta values
  betas<-diag(length(impdat[1,]))
  betas[,misnumer]<-olsout[[2]]
  betas[,misfac]<-logisticout[[2]]
  betas[,miscount]<-countout[[2]]
  betas[is.na(betas)]<-0
  betas<-betas[,-1]
  
  
    # Return data to original state and end function
    colnames(ml)<-orgorder2
  out<-list(ml,mispat2,sse,idfy,vs,betas,xtx,levnames3)
  names(out)<-c("ML","Missingness Pattern","SSE","TYPE","Variances","Betas","X'X","facnames")
  print(it)
  return(out)
}









