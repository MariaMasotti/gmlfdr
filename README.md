# gmlfdr
 - An R package for fitting gaussian mixture for use in locfdr calculations 
 - Sample R codes
```R
 ## install the package
 devtools::install_github('MariaMasotti/gmlfdr')
 library(gmlfdr)
 ## simulate null and alternative data
null<-rnorm(1000,mean=0,sd=1)
alt<-rnorm(1000,mean=3,sd=1)
#set number of components in alternate dist
k<-1
#set initial estimates
pi0<-c(.5,.5)
mu0<-c(0,0)
sd0<-1
#run em algorithm to obtain parameter extimates 
gmmr(null,alt,k,pi0,mu0,sd0)

```
