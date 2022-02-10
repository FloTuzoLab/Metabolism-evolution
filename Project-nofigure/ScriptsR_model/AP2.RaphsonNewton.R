require("Rmpfr")
options(digits=20)
##For the moment, only 2 variables
Fxy=function(Fe,x,y) eval(Fe)#evaluating the value of a function
#Implementation of Raphson-Newton process
RaphsonNewton<-function(F_list,init,error){
  x_0=init[1]#should be different from 0
  y_0=init[2]#should be different from 0
  x=x_0
  y=y_0
  x_d=x
  y_d=y
  dF<-list()
  dF[[1]]<-deriv(F_list[[1]],c("x","y"),func=T)
  dF[[2]]<-deriv(F_list[[2]],c("x","y"),func=T)
  while((x_d/x>error | y_d/y>error) | (x_d>error | y_d>error) ){# 
    #F1c<-mpfr(Fxy(F_list[[1]],x,y),80)
    F1c<-Fxy(F_list[[1]],x,y)
    #F2c<-mpfr(Fxy(F_list[[2]],x,y),80)
    F2c<-Fxy(F_list[[2]],x,y)
    #dF1_dx=mpfr(attr(dF[[1]](x,y),"gradient")[1,1],80)
    dF1_dx=attr(dF[[1]](x,y),"gradient")[1,1]
    #dF1_dy=mpfr(attr(dF[[1]](x,y),"gradient")[1,2],80)
    dF1_dy=attr(dF[[1]](x,y),"gradient")[1,2]
    #dF2_dx=mpfr(attr(dF[[2]](x,y),"gradient")[1,1],80)
    dF2_dx=attr(dF[[2]](x,y),"gradient")[1,1]
    #dF2_dy=mpfr(attr(dF[[2]](x,y),"gradient")[1,2],80)
    dF2_dy=attr(dF[[2]](x,y),"gradient")[1,2]
    #y_d=mpfr((-F2c+F1c*dF2_dx/dF1_dx)/(dF2_dy-dF2_dx*dF1_dy/dF1_dx),80)
    y_d=(-F2c+F1c*dF2_dx/dF1_dx)/(dF2_dy-dF2_dx*dF1_dy/dF1_dx)
    #x_d=mpfr((-F1c-y_d*dF1_dy)/dF1_dx,80)
    x_d=(-F1c-y_d*dF1_dy)/dF1_dx
    #x=as.numeric(mpfr(x+x_d,80))
    x=x+x_d
    #y=as.numeric(mpfr(y+y_d,80))
    y=y+y_d
    i=i+1
  }
  res<-list("x"=x,"y"=y)
  return(res)
}

