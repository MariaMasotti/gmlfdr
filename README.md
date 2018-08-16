# gmlfdr
 - An R package for fitting gaussian mixture for use in locfdr calculations 
 - Sample R codes
```R
 ## install the package
 devtools::install_github('MariaMasotti/gmlfdr')
 library(gmlfdr)
#example on locfdr hivdata
library(locfdr)
data("hivdata")
data<-hivdata
ac<-data[data<quantile(data,.25)|data>quantile(data,.75)]
a<-data[data>=quantile(data,.25)&data<=quantile(data,.75)]
pi0<-c(.5,.25,.25)
mu0<-c(0,-3,3)
sd0<-1
k<-2
gmmr(a,ac,k,pi0,mu0,sd0)


```
