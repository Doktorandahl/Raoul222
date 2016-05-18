library(Jmisc)
library(fBasics)
library(gtools)
library(ppls)
## Main-Function

raoul<-function(x,m=5,tol=0.0001,maxit=20,facs=NULL,counts=NULL,returncat=FALSE,returncount=FALSE){
  y<-seqreg(x,tol,maxit,facs,counts)
  mi<-raoul.MI(y,m,returncat,returncount,x,facs)
  out<-list(as.data.frame(y[[1]]),mi,y[-1])
  names(out)<-c("Maximum Likelihood","Multiple Imputations","Descriptives")
  class(out)<-"raoul"
  return(out)
}

## Function for extracting mutiple imputed datasets, called by main function

raoul.MI<-function(y,m,returncat,returncount,x,facs){
  sets<-NULL
  for(k in 1:m){
  sets<-c(sets,list(as.data.frame(raoul.imputations(y,returncat,returncount,x,facs))))
  }
  names(sets)<-paste(rep("set",m),seq(1:m),sep="")
  return(sets)
}

## Function for making imputatated datasets, called by raoul.MI
raoul.imputations<-function(y,returncat,returncount,x,facs){
  z<-cbind(1,y[[1]])
  k<-length(z[1,])-1
  n<-length(z[,1])
  s<-y[[3]][1,]/rchisq(k,(y[[3]][2,]-(k+1)))*n/y[[3]][2,]
  bimp<-y[[6]]
  for(i in 1:k){
    ind<-bimp[,i]!=0
    bimp[which(ind==TRUE),i]<-bimp[which(ind==TRUE),i]+sqrt(s[i])*rnorm(sum(ind),1,0)%*%chol(solve(t(z[,which(ind==TRUE)])%*%z[,which(ind==TRUE)]))
    
  }
  set2<-set1<-z%*%bimp ## Matrix of fitted values from random parameters
  for(i in 1:k){
    if(y[[4]][1,i]=="1"|y[[4]][1,i]=="2"){ ## Adding normal error 
      set2[,i]<-set1[,i]+rnorm(n,0,1)*sqrt(s[i])
    }
  }
  set3<-set2
  i<-1
  while(i<=k){
    if(y[[4]][1,i]=="2"){ ## Converting odds ratio to probability
      nl<-sum(y[[4]][2,i]==y[[4]][2,]) ## Identify numer of factor levels
      for(l in i:(i+nl-1)){ ## Converting to probabilities
        set3[,l]<-exp(set2[,l])/as.vector((1+as.matrix(apply(as.matrix(exp(set2[,i:(i+nl-1)])),1,sum))))
        set3[which(is.nan(set3[,l])),l]<-1-1e-10
        set3[which(set3[,l]==0),l]<-1e-10
      }
      i<-i+nl-1

    }
    if(y[[4]][1,i]=="3"){
      set3[,i]<-exp(set2[,i])
    }
    i<-i+1
  }
  if(returncount){
    for(i in 1:k){
      if(y[[4]][1,i]=="3"){
        set3[,i]<-rpois(n,set3[,i])
      }
    }
  }
    if(returncat){
      facs2<-facs[which(apply(is.na(dat[,facs]),2,max)==0)]
      i<-1
      s<-0
      q<-0
      rems=NULL
      if(sum(y[[4]][1,]=="2")!=0){
      repl<-matrix(,n,length(unique(y[[4]][2,which(y[[4]][1,]=="2")])))
      colnames(repl)<-1:length(repl[1,])
      }else{
        repl<-NULL
      }
      if(sum(y[[4]][1,]=="5")!=0){
      repl2<-matrix(,n,length(unique(y[[4]][2,which(y[[4]][1,]=="5")])))
      colnames(repl2)<-1:length(repl2[1,])
      }else{
        repl2<-NULL
      }
      while(i<=k){
        
      if(y[[4]][1,i]=="2"){ ## Converting odds ratio to probability
        s<-s+1
        nl<-sum(y[[4]][2,i]==y[[4]][2,]) ## Identify numer of factor levels
        run<-runif(n,0,1)
        
        for(p in 1:nl){
          if(p==1){
            repl[which(is.na(repl[,s])),s]<-ifelse(set3[which(is.na(repl[,s])),i+p-1]<run[which(is.na(repl[,s]))],y[[8]][2,s],NA)
          }
          if(p==nl){
            repl[which(is.na(repl[,s])),s]<-ifelse(set3[which(is.na(repl[,s])),i+p-1]<run[which(is.na(repl[,s]))],y[[8]][nl+1,s],NA)
            repl[which(is.na(repl[,s])),s]<-ifelse(set3[which(is.na(repl[,s])),i+p-1]>run[which(is.na(repl[,s]))],y[[8]][1,s],NA)
          }
          if(p!=nl & p!=1){
            repl<-ifelse(set3[which(is.na(repl[,s])),i+p-2]>run[which(is.na(repl[,s]))] & set3[,i+p-1]<run[which(is.na(repl[,s]))],y[[8]][s,p+1],NA)
          }
        }
        colnames(repl)[s]<-y[[4]][2,i]
        i<-i+nl
      }else if(y[[4]][1,i]=="5"){
        q<-q+1
        nl<-sum(y[[4]][2,i]==y[[4]][2,]) ## Identify numer of factor levels
        repl2[,q]<-paste(dat[,facs2[q]])
        colnames(repl2)[q]<-colnames(x)[facs2[q]]
        i<-i+nl
        
        }else{
          rems<-c(rems,i)
          i<-i+1
        }
      }
      colnames(set3)<-1:length(set3[1,])
      colnames(set3)[rems]<-colnames(y[[4]][,which(y[[4]][1,]!="2" & y[[4]][1,]!="5")])
      set4<-cbind(set3[,rems],repl,repl2)
      set4<-set4[,colnames(x)]
      set4[!is.na(x)]<-x[!is.na(x)]
      return(set4)
    }else{
      set3[y[[2]]==0]<-y[[1]][y[[2]]==0] ## Replacing fitted values with observed
      colnames(set3)<-colnames(y[[1]])
      return(set3)
    }
    
  
  
  
 
  
}


## Function for making imputatated datasets w. Dirichlet, called by raoul.MI
## NOT USED!
raoul.imputations.dir<-function(y){
  z<-cbind(1,y[[1]])
  k<-length(z[1,])-1
  n<-length(z[,1])
  s<-y[[3]][1,]/rchisq(k,(y[[3]][2,]-(k+1)))*n/y[[3]][2,]
  bimp<-y[[6]]
  for(i in 1:k){
    ind<-bimp[,i]!=0
    bimp[which(ind==TRUE),i]<-bimp[which(ind==TRUE),i]+sqrt(s[i])*rnorm(sum(ind),1,0)%*%chol(solve(t(z[,which(ind==TRUE)])%*%z[,which(ind==TRUE)]))
    
  }
  set2<-set1<-z%*%bimp ## Matrix of fitted values from random parameters
  
  for(i in 1:k){
    if(y[[4]][1,i]=="1"){ ## Adding normal error term to numeric
      set2[,i]<-set1[,i]+rnorm(n,0,1)*s[i]
    }
    if(y[[4]][1,i]=="2"){ ## Converting odds ratio to probability and drawing dirichlet
      nl<-sum(y[[4]][2,i]==y[[4]][2,]) ## Identify numer of factor levels
      for(l in i:(i+nl-1)){ ## Converting to probabilities
        set2[,l]<-exp(set1[,l])/as.vector((1+as.matrix(apply(as.matrix(exp(set1[,i:(i+nl-1)])),1,sum))))
      }
      set2[,i:(i+nl-1)]<-t(as.matrix(apply(cbind(as.matrix(set2[,i:(i+nl-1)]),as.matrix(1-apply(as.matrix(set2[,i:(i+nl-1)]),1,sum))),1,rdirichlet,n=1)))[,-(nl+1)] ## Imputing Dirichlet draws
      i<-i+nl-1
      
    }
  }
  
  set2[y[[2]]==0]<-y[[1]][y[[2]]==0]
  
  return(set2)
  
}



## Function for creating (single) imputed datasets, called by raoul.MI
## NOT USED!
raoul.randomimps<-function(y){
  z<-y[[1]]
  for(i in 1:length(z[1,])){
    if(y[[4]][i]==1){
      z[which(y[[2]][,i]==1),i]<-rnorm(sum(y[[2]][,i]),z[which(y[[2]][,i]==1),i],y[[3]][i])
    } else if(y[[4]][i]==2){
      z[which(y[[2]][,i]==1),i]<-rbeta(sum(y[[2]][,i]),
                                          z[which(y[[2]][,i]==1),i]/y[[3]][i],
                                          (1-z[which(y[[2]][,i]==1),i])/y[[3]][i])
    }
      
  }
  return(z)
  
}  
  
## Fitting function for numeric variables. Called by EM-function
  
raoul.olsfit<-function(mis,full,m,it){
  
    fits<-NULL
   b<-matrix(,length(full[1,]),length(m))
    if(it==0){
      if(length(full[1,])>1){
      b<-apply(mis,2,raoul.b1,obs=full)
      }else{
        b<-t(as.matrix(apply(mis,2,raoul.b1,obs=full)))
        
      }
    for(r in 1:length(mis[1,])){
      fits<-cbind(fits,full%*%b[,r])
    }
    }
    else{
      for(r in 1:length(mis[1,])){
        b[-m[r],r]<-solve(t(full[,-m[r]])%*%full[,-m[r]])%*%t(full[,-m[r]])%*%mis[,r]
        fits<-cbind(fits,full[,-m[r]]%*%b[-m[r],r])
      }
    }
    out<-list(fits,b)
  return(out)
}

## USED?

raoul.imp<-function(mis,full,b){
  x<-full[,which(apply(full,2,identical, y=mis)==FALSE)]
  b2<-b[,which(apply(misall,2,identical, y=mis)==TRUE)]
  mis1<-misorg[,which(apply(misall,2,identical, y=mis)==TRUE)]
  mis3<-mis1
  mis2<-x%*%b2
  mis3[which(is.na(mis1)==TRUE)]<-mis2[which(is.na(mis1)==TRUE)]
  return(mis3)
}

## Function for calculating MSE of the imputation-regressions. Called by EM-function

raoul.sse<-function(fitted,original){
  if(is.null(fitted)){
    return(NULL)
  }
  sse<-NULL
  for(i in 1:length(fitted[1,])){
    sse<-cbind(sse,sum((fitted[which(is.na(original[,i])==FALSE),i]-original[which(is.na(original[,i])==FALSE),i])^2))
  }
  return(sse)
}



## Function for extracting b-values. Called by OLS-fit function

raoul.b1<-function(mis, obs){
  b1<-(solve(t(obs[which(is.na(mis)==FALSE),])%*%obs[which(is.na(mis)==FALSE),])%*%
      t(obs[which(is.na(mis)==FALSE),])%*%mis[which(is.na(mis)==FALSE)])
  return(b1)
}

## USED?

raoul.imp1<-function(obs,b1,mis){
  mis3<-mis
  mis2<-obs%*%b1
  mis3[which(is.na(mis)==TRUE)]<-mis2[which(is.na(mis)==TRUE)]
  return(mis3)
}

## Function for extracting b-values. Called by OLS-fit function
raoul.b2<-function(mis,full){
  b2<-solve(t(full)%*%full)%*%t(full)%*%mis
  return(b2)
}

## USED?
raoul.imp2<-function(mis,full,b2,misall,misorg){
  x<-full[,which(apply(full,2,identical, y=mis)==FALSE)]
  b<-b2[,which(apply(misall,2,identical, y=mis)==TRUE)]
  mis1<-misorg[,which(apply(misall,2,identical, y=mis)==TRUE)]
  mis3<-mis1
  mis2<-x%*%b
  mis3[which(is.na(mis1)==TRUE)]<-mis2[which(is.na(mis1)==TRUE)]
  return(mis3)
}


## Function for evaluating likelihood. Called by EM-function
raoul.lik<-function(mis,full,b2,misall,n1){
  x<-full[,which(apply(full,2,identical, y=mis)==FALSE)]
  b<-as.matrix(b2[,which(apply(misall,2,identical, y=mis)==TRUE)])
  s<-var(mis)
  lik<-(-n1/2)*log(2*pi*s)-1/(2*s)*sum((mis-x%*%b)^2)
  return(lik)
}

## Function for evaluating likelihood. Called by EM-function
raoul.lik2<-function(full,n){
  x<-full[,-1]
  s<-cov(x)
  lik<- n*log(det(s))+sum(diag(solve(s)%*%t(demean(x))%*%demean(x)))
  return(lik)
  }


## Function for identifying factors. Called by EM-function
raoul.factors<-function(dat){
  m<-length(dat[1,])
  facind<-NULL
  for(i in 1:m){
    if(is.numeric(dat[1,i])==FALSE){
      facind<-cbind(facind,i)
    }
  }
  return(facind)
}

## Function for identifying factors. Called by EM-function
raoul.counts<-function(dat){
  m<-length(dat[1,])
  countind<-NULL
  for(i in 1:m){
    if(class(dat[,i])=="integer"){
      countind<-cbind(countind,i)
    }
  }
  return(countind)
}
  
## Function for creating dummies. Called by EM-function

raoul.dummy<-function(data,facts){
  k<-length(facts)
  facid<-matrix(0,1,length(data[1,]))
  for(i in 1:k){
    data<-cbind(data,model.matrix(~facts[,i]+0))
    facid<-cbind(facid,matrix(i,1,length(model.matrix(~facts[,i]+0))))
  }
}



## USED?

raoul.logistic2<-function(misfac,full,n1,it){
  facnames<-unique(colnames(misfac)) ## Unique factors with missing values
  misout<-NULL ## Empty out-matrix
  for(i in 1:length(facnames)){ ## Return all columns which should be imputed
    facimp<-misfac[,which(colnames(misfac)==facnames[i])] ## Factors to be imputed
    r<-length(facimp[1,])
    
    if(it>0){
    full2<-matrix(full[!(colnames(full) %in% facnames[i])],n1,length(full[1,])-r) ## Drop factors to be imputed if not first imputation
    } else {
    full2<-full[which(is.na(facimp[,1])==FALSE),] ## Drop missing values if first imputation
    facimp<-facimp[which(is.na(facimp[,1])==FALSE),] ## Drop missing values if first imputation
    }
    coffs<-matrix(,length(full2[1,]),r)
  for(l in 1:r){ ## Obtain betas
      y<-facimp[,l]
      coffs[,l]<-raoul.logest(full2,y,length(full2[,1]),20,0.0001)
      
      
  }
  for(l in 1:r){ ## Impute fitted values
    if(it>0){
    misout<-cbind(misout,exp(full2%*%coffs[,l])/(1+apply(exp(full2%*%coffs),1,sum)))
    } else{
      misout<-cbind(misout,exp(full%*%coffs[,l])/(1+apply(exp(full%*%coffs),1,sum)))
    }
  }
    
  } 
  return(misout)   
  
  }
  
## Called by EM-main
raoul.logistic1<-function(misfact,full,m,n1,it,levlen,misp,ice){
  
  facnames<-unique(colnames(misfact)) ## Unique factors with missing values
  levlen<-c(1,levlen)
  imps<-NULL ## Empty out-matrix
  sse<-NULL
  if(it>0){
  b<-matrix(,length(full[1,]),length(m))
  for(i in 1:length(facnames)){
    logout<-raoul.logistic3(facimp=as.matrix(misfact[,which(colnames(misfact)==facnames[i])]),full=as.matrix(full),n1=n1,it=it,misp=as.matrix(misp[,which(colnames(misfact)==facnames[i])]),facnames=facnames[i])
    if(!is.list(logout)){
      return(ice)
    }
    imps<-cbind(imps,logout[[1]])
    b[-m[sum(levlen[1:i]):(sum(levlen[1:(i+1)])-1)],sum(levlen[1:i]):(sum(levlen[1:(i+1)])-1)]<-logout[[2]]
    sse<-cbind(sse,logout[[3]])
  }
  
  }else{
    b<-NULL
    for(i in 1:length(facnames)){
      logout<-raoul.logistic3(facimp=as.matrix(misfact[,which(colnames(misfact)==facnames[i])]),full=as.matrix(full),n1=n1,it=it,facnames=facnames[i],misp=misp)
      imps<-cbind(imps,logout[[1]])
      b<-cbind(b,logout[[2]])
      sse<-cbind(sse,logout[[3]])
    }
  }
  out<-list(imps,b,sse)
  return(out)
}
  
## Called by logistic 1
raoul.logistic3<-function(facimp,full,n1,it,facnames,misp){
    misout<-NULL
    r<-length(facimp[1,])
    
    if(it>0){
      full2<-as.matrix(full[,which((colnames(full) %in% facnames)==FALSE)]) ## Drop factors to be imputed if not first imputation
    } else {
      full2<-as.matrix(full[which(is.na(facimp[,1])==FALSE),]) ## Drop missing values if first imputation
      facimp<-as.matrix(facimp[which(is.na(facimp[,1])==FALSE),]) ## Drop missing values if first imputation
    }
    coffs<-matrix(,length(full2[1,]),r)
    sse<-matrix(,1,r)
    for(l in 1:r){ ## Obtain betas
      y<-facimp[,l]
      logout<-raoul.logest(full=full2,y=y,n=length(full2[,1]),maxit=20,toll=0.0001,misp2=misp[,l],it=it)
      if(!is.list(logout)){
        return("ice")
      }
      coffs[,l]<-logout[[1]]
      sse[,l]<-logout[[2]]
      
      
    }
    for(l in 1:r){ ## Impute fitted values
      if(it>0){
        misout<-cbind(misout,exp(full2%*%coffs[,l])/(1+apply(exp(full2%*%coffs),1,sum)))
      } else{
        misout<-cbind(misout,exp(full%*%coffs[,l])/(1+apply(exp(full%*%coffs),1,sum)))
      }
    }
    
 
  out<-list(misout,coffs,sse)
  return(out)   
  
}


## IWLS-estimation function. Called by Logistic3-function
raoul.logest<-function(full,y,n,maxit,toll,misp2,it){ ### IWLS
  b1<-matrix(0,length(full[1,]),1)
  b0<-b1+toll
  its<-1
  while(its<maxit & (abs(sum((b1+toll)/(b0+toll))-length(b1))>toll)){
    b0 <- b1
    eta<-full%*%b1
    mu <- 1/(1 + exp(-eta))
    nu <- (mu*(1 - mu))
    nu[which(nu==0)]<-nu[which(nu==0)]+1e-50
    w <- as.vector(n*nu)
    z <- eta + (y - mu)/nu
    if((abs(det(t(full*w)%*%full))>1e+50 | abs(det(t(full*w)%*%full))<1e-15) & it>1){
      return("singular")
    }
    b1 <- solve(t(full*w)%*%full)%*%t(full*w)%*%z
    its <- its+1
  }
  if(it==0){
    sse<-sum((z-full%*%b1)^2)
  }else{
  sse<-sum((z[which(misp2==0)]-full[which(misp2==0),]%*%b1)^2)
  }
  out<-list(b1,sse)
  return(out)
}

## Missing pattern function. Called by EM-function

raoul.mispat<-function(x){
  z<-as.matrix(apply(x,2,is.na))
  return(z*1)
}


raoul.countfit<-function(mis,full,m,it,misp,ice){
  imps<-matrix(,length(mis[,1]),length(mis[1,]))
  sse<-matrix(,1,length(m))
  b<-matrix(,length(full[1,]),length(m))
  if(it==0){
    for(k in 1:length(m)){
      posout<-raoul.poest(full2=as.matrix(full[which(!is.na(mis[,k])),]),y=as.matrix(mis[which(!is.na(mis[,k])),k]),maxit=1000,toll=0.0001,misp2=misp,it=it)
      b[,k]<-posout[[1]]
    }
    imps<-exp(full%*%b)
  } else{
    for(k in 1:length(m)){
      posout<-raoul.poest(full2=as.matrix(full[,-m[k]]),y=as.matrix(mis[,k]),maxit=1000,toll=0.0001,misp2=misp[,k],it=it)
      if(!is.list(posout)){
        return(ice)
      }
      b[-m[k],k]<-posout[[1]]
      imps[,k]<-exp(full[,-m[k]]%*%b[-m[k],k])
      sse[,k]<-posout[[2]]
    }
    
  }
  out<-list(imps,b,sse)
  return(out)
  }


## IWLS-estimation function for poisson. Called by Raoul.Countfit
raoul.poest<-function(full2,y,maxit,toll,misp2,it){ ### IWLS
  b1<-matrix(0,length(full2[1,]),1)
  b0<-b1+toll
  its<-1
  while(its<maxit & (abs(sum((b1+toll)/(b0+toll))-length(b1))>toll)){
    b0 <- b1
    eta<-full2%*%b1
    mu <- exp(eta)
    w <- normalize.vector(as.vector(mu))
    z <- eta + (y - mu)/mu
    if(abs(det(t(full2*w)%*%full2))>1e+85 & it>1){
      return("singular")
    }
    b1 <- solve(t(full2*w)%*%full2)%*%t(full2*w)%*%z
    its <- its+1
  }
  if(it==0){
    sse<-sum((z-full2%*%b1)^2)
  }else{
    sse<-sum((z[which(misp2==0),]-full2[which(misp2==0),]%*%b1)^2)
  }
  out<-list(b1,sse)
  return(out)
}

raoul.lm<-function(formula, raoul, subset=NULL, w=NULL, na.action=NULL,
                   method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
                   singular.ok = TRUE, contrasts = NULL, offset=NULL){
  if(class(raoul)!="raoul"){
    print("Object not Raoul, running Raoul to impute missing values")
    raoul<-raoul(raoul)
  }
  m<-length(raoul[[2]])
  desks<-NULL
  r2s<-c(0,0)
  for(i in 1:m){
    md<-lm(formula, raoul[[2]][[i]], subset=NULL, w=NULL, na.action=NULL,
           method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
           singular.ok = TRUE, contrasts = NULL, offset=NULL)
    desks<-c(desks,list(summary(md)$coef))
    r2s<-r2s+c(summary(md)[[8]],summary(md)[[9]])
    rnm<-rownames(summary(md)$coef)
    cnm<-colnames(summary(md)$coef)
    
  }
  desks<-array(unlist(desks), dim = c(nrow(desks[[1]]), ncol(desks[[1]]), length(desks)))
  b<-rowMeans(desks[,1,])
  se<-sqrt(1/m*apply(desks[,2,]^2,1,sum)+(1+1/m)*(1/(m-1))*apply((desks[,1,]-b)^2,1,sum))
  t<-b/se
  p<-pt(t,length(raoul[[2]][[1]][,1]))
  avgr2<-r2s/m
  coefs<-cbind(b,se,t,p)
  rownames(coefs)<-rnm
  colnames(coefs)<-cnm
  names(avgr2)<-c("Avg.R-sq.", "Avg.Adj.R-sq.")
  out<-list(coefs,avgr2)
  names(out)<-c("Coefficients","Average R-squared")
  class(out)<-"raoul.lm"
  return(out)
}



raoul.glm<-function(formula, raoul,family="gaussian", subset=NULL, w=NULL, na.action=NULL,
                   method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
                   singular.ok = TRUE, contrasts = NULL, offset=NULL){
  if(class(raoul)!="raoul"){
    print("Object not Raoul, running Raoul to impute missing values")
    raoul<-raoul(raoul)
  }
  m<-length(raoul[[2]])
  desks<-NULL
  r2s<-c(0,0)
  for(i in 1:m){
    md<-lm(formula, raoul[[2]][[i]],family, subset=NULL, w=NULL, na.action=NULL,
           method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
           singular.ok = TRUE, contrasts = NULL, offset=NULL)
    desks<-c(desks,list(summary(md)$coef))
    r2s<-r2s+c(summary(md)[[8]],summary(md)[[9]])
    rnm<-rownames(summary(md)$coef)
    cnm<-colnames(summary(md)$coef)
    
  }
  desks<-array(unlist(desks), dim = c(nrow(desks[[1]]), ncol(desks[[1]]), length(desks)))
  b<-rowMeans(desks[,1,])
  se<-sqrt(1/m*apply(desks[,2,]^2,1,sum)+(1+1/m)*(1/(m-1))*apply((desks[,1,]-b)^2,1,sum))
  t<-b/se
  p<-pt(t,length(raoul[[2]][[1]][,1]))
  avgr2<-r2s/m
  coefs<-cbind(b,se,t,p)
  rownames(coefs)<-rnm
  colnames(coefs)<-cnm
  names(avgr2)<-c("Avg.R-sq.", "Avg.Adj.R-sq.")
  out<-list(coefs,avgr2)
  names(out)<-c("Coefficients","Average R-squared")
  class(out)<-"raoul.lm"
  return(out)
}
