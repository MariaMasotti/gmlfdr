//
//  main.cpp
//  
//
//  Created by Maria on 8/6/18.
//  
//

#include <Rcpp.h>
#include <iostream>
#include <cstdlib>
#include <algorithm> 
#include <array> 
using namespace std;
using namespace Rcpp;

//input null data, alternative data,k=number of components of alternative dist, and initial estimates of pi,mu,sd
//export list 
// [[Rcpp::export]]
List gmmr(NumericVector data,
          NumericVector datac, 
          int k, 
          NumericVector pi0, 
          NumericVector mu0, 
          double sd0) {
  
  //create vectors holding the initial values
  NumericVector pi_new=pi0;
  NumericVector mu_new=mu0;
  double sd_new=sd0;
  int K=k+1;
  
  //create nonzero difference between new and old parameters
  double l=1.0;
  NumericVector v(K);
  NumericVector w(K);
  NumericVector x(K);
  v=rep(l,K);
  w=rep(l,K);
  x=l;
  
  //initiate stuff
  NumericVector pi_old;
  NumericVector mu_old;
  double sd_old;
  NumericMatrix mat(datac.size(),K);
  NumericMatrix t(datac.size(),K);
  NumericMatrix test(datac.size(),K);
  double var_new;
  int j=0;
  
  //em loop
  for(int b=0;b<1000;b++){
    
    //set old parameters equal to new parameters  
    pi_old=clone(pi_new);
    mu_old=clone(mu_new);
    sd_old=sd_new;
    
    
    //calculate t using old parameters 
    for(int a=0;a<K;a++){
      mat(_,a)=pi_old[a]*dnorm(datac,mu_old[a],sd_old);
    }
    for(int a=0;a<K;a++){
      t(_,a)=(pi_old[a]*dnorm(datac,mu_old[a],sd_old))/rowSums(mat);  
    }
    
    //calculate new pi
    pi_new[0]=(data.size()+colSums(t)[0])/(data.size()+sum(colSums(t)));
    for(int a=1;a<K;a++){
      pi_new[a]=(colSums(t)[a])/(data.size()+sum(colSums(t)));
    }
    
    //calculate new mu
    mu_new[0]=(sum(data)+sum(datac*t(_,0)))/(data.size()+colSums(t)[0]);
    for(int a=1;a<K;a++){
      mu_new[a]=(sum(datac*t(_,a)))/(colSums(t)[a]);
    }
    
    //calculate new sd
    for(int a=0;a<K;a++){
      test(_,a)=t(_,a)*pow((datac-mu_old[a]),2);
    }
    var_new=(sum(pow((data-mu_old[0]),2))+ sum(rowSums(test)))/(data.size()+datac.size());
    sd_new=sqrt(var_new);
    
    
    //calc distance between old and new parameters
    NumericVector pidiff=pi_new-pi_old;
    NumericVector mudiff=mu_new-mu_old;
    double sddiff=sd_new-sd_old;
    v=abs(pidiff);
    w=abs(mudiff);
    x=abs(sddiff);
    
    
    //iteration counter
    j=j+1;
    
    //break loop if all parameters converge 
    if(is_true(all(v<.00001))&&is_true(all(w<.00001))&&is_true(all(x<.00001)))
      break;
  }
  
  //if loop doesnt break print error message 
  if(j>999){
    Rprintf("Error: No convergence\\n");
  }
  
  //calc number of parameters
  int p = pi_new.size()-1 + mu_new.size()+1;
  
  //calc marginal likelihood
  NumericMatrix b(datac.size(),K);
  for(int i=0;i<K;i++){
    b(_,i)=pi_new[i]*dnorm(datac,mu_new[i],sd_new);
  }
  double ml;
  ml=sum(log(pi_new[0]*dnorm(data,mu_new[0],sd_new)))+sum(log(rowSums(b)));
  
  //calc aic and bic
  double aic=2*p-2*ml;
  double bic=std::log(data.size()+datac.size())*p-2*ml;  
  
  //output list of all final parameters    
  List all;
  all["pi"]=pi_new;
  all["mu"]=mu_new;
  all["sd"]=sd_new;
  all["k"]=K-1;
  all["iterations"]=j;
  all["aic"]=aic;
  all["bic"]=bic;
  return all;
}