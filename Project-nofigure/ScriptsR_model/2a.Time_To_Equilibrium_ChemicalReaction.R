source("AP1.Plotting_multiple_images.R")
data<-read.table(file="Time_To_Equilibrium/Results/line24_timeToEq.rj",header=T)
#Defining features of plots
ncol=64
ncontour=5
pal<-colorRampPalette(c("green","red","hot pink"))(ncol)

#Time to reach the equilibrium - Chemical reaction only - using C++ to limit numercial problems


data<-list()
data[[1]]<-read.table(file="Time_To_Equilibrium/Results/line0_timeToEq.rj",header=T)
data[[2]]<-read.table(file="Time_To_Equilibrium/Results/line1_timeToEq.rj",header=T)
##Old: importing values from C++
N_reso=50
DeltaT_eq<-list()
DeltaT_eq[[1]]<-matrix(nrow=N_reso , ncol=N_reso)
DeltaT_eq[[2]]<-matrix(nrow=N_reso , ncol=N_reso)
for (i in 1:N_reso){
  for (j in 1:N_reso){
    DeltaT_eq[[1]][i,j]<-data[[1]]$T_eq[data[[1]]$i==(i-1) & data[[1]]$j==(j-1)]
    DeltaT_eq[[2]][i,j]<-data[[2]]$T_eq[data[[2]]$i==(i-1) & data[[2]]$j==(j-1)]
  }
}

pal<-colorRampPalette(c("green","red","hot pink"))(ncol)
par(mfrow=c(1,1))
image.plot(DeltaT_eq[[1]],nlevel=ncol,col=pal,axes=FALSE)
contour(DeltaT_eq[[1]], add = TRUE,levels=c(0.1,1))
#contour(DeltaT_eq[[2]], add = TRUE,levels=c(0.1,1)) to show the contour is identical
axis(1, at=seq(0,1,0.33),labels=seq(2,8,2))
axis(2, at=seq(0,1,0.33),labels=seq(-2,4,2))

#Useless as the fitness landscape is exactly the same
image.plot(DeltaT_eq[[2]],nlevel=ncol,col=pal,axes=FALSE)
contour(DeltaT_eq[[2]], add = TRUE,levels=c(0.1,1))
axis(1, at=seq(0,1,0.33),labels=seq(2,8,2))
axis(2, at=seq(0,1,0.33),labels=seq(-2,4,2))


###R simulations stuck with slow performance


N_reso=31
#Creation of two sets of parameters
log10kf_var<-c(1,7,1)
log10kcat_var<-c(-1,5,1)
k2f<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
k2cat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
kr=10^3



log10_cSf<-matrix(nrow=N_reso , ncol=N_reso)
cPi<-matrix(nrow=N_reso , ncol=N_reso)
E_tot_conc=10^-3
P1_Phi_eq=10^-6#c(10^-6,10^-3)

#Function determining the time to reach the equilibrium
#dt=0.00001 #timestep in seconds
T_eq<-function(Conc_eq_ratio,t_max){
  DeltaT_eq<-matrix(nrow=N_reso , ncol=N_reso)
  for (i in 1:N_reso){
    for (j in 1:N_reso){
      print(c(i,j))
      dt<-min(1/(2*k2cat[j]),1/(2*kr),1/(5*k2f[i]*P1_Phi_eq),0.01)
      cES=0
      cEf=E_tot_conc
      cSf=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
      ##Initalization of the process
      cP=0
      cSf_c<-cSf
      cEf_c<-cEf
      cES_c<-cES
      t=0
      while((cP/(P1_Phi_eq*dt)<Conc_eq_ratio) && (t<t_max)){
        cSf_c<-cSf
        cEf_c<-cEf
        cES_c<-cES
        FD_dt=P1_Phi_eq#Vm*(Se-cSf_c)/(KT+(Se+cSf_c)+Se*cSf_c/KT) ##Facilitated diffusion
        CR1_dt=k2f[i]*cSf_c*cEf_c-(kr)*cES_c
        CR2_dt=k2cat[j]*cES_c
        cSf=cSf_c+(FD_dt-CR1_dt)*dt 
        cEf=cEf_c+(CR2_dt-CR1_dt)*dt
        cES=cES_c-(CR2_dt-CR1_dt)*dt
        cP=k2cat[j]*cES_c*dt
        DeltaT_eq[i,j]<-t
        t=t+dt
      }
    }
  }
  return(DeltaT_eq)
}

Time_to_eq_CR<-list()
for (p in 1:length(P1_Phi_eq)){
  Time_to_eq_CR[[p]]<-T_eq(0.99,10)
}


x_lim=c(2,6,2)
y_lim=c(0,4,2)

image.plot(Time_to_eq_CR[[1]],nlevel=ncol,col=pal,axes=FALSE)
axis(1, at=seq(0,1,0.5),labels=seq(2,6,2))
axis(2, at=seq(0,1,0.5),labels=seq(0,4,2))
multiplePlot("kr=","","Vm=","M/s",kr,P1_Phi_eq,log10kf_var,log10kcat_var,Time_to_eq_CR[[1]],
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",ncol=ncol,palette=pal)

delta_T_low<-T_eq(0.9,30)
delta_T_high<-T_eq(0.95,15)

xlim<-log10kf_var
ylim<-log10kcat_var
pal<-colorRampPalette(c("green","red","hot pink"))(ncol)
##Need to be contrasted with equilibriums
par(mfrow=c(1,2))

axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]))
axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]))
image.plot(delta_T_low/delta_T_high,nlevel=ncol,col=pal,axes=F)
axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]))
axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]))


image.plot(log10_T_eq,nlevel=ncol,col=pal,axes=F)
axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]))
axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]))
image.plot(log10_T_eq,nlevel=ncol,col=pal,axes=F)
contour(log10_T_eq, add = TRUE,nlevels=ncontour)


k21<-10^2
k22<-10^3
dt=0.0005 #timestep in seconds
for (i in 1:N_reso){
  k_1=0
  for (j in 1:N_reso){
    cES=0
    cEf=E_tot_conc
    cSf=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
    ##Initalization of the process
    cP=0
    cSf_c<-cSf
    cEf_c<-cEf
    cES_c<-cES
    t=0
    while((cP/(P1_Phi_eq*dt)<0.999) && (t<50)){
      t=t+dt
      print(paste("Substrate conc:",cSf_c))
      print(paste("Product flux:",cP))
      cSf_c<-cSf
      cEf_c<-cEf
      cES_c<-cES
      FD_dt=P1_Phi_eq#Vm*(Se-cSf_c)/(KT+(Se+cSf_c)+Se*cSf_c/KT) ##Facilitated diffusion
      CR1_dt=k21*cSf_c*cEf_c-(k_1)*cES_c
      CR2_dt=k22*cES_c
      cSf=cSf_c+(FD_dt-CR1_dt)*dt 
      cEf=cEf_c+(CR2_dt-CR1_dt)*dt
      cES=cES_c-(CR2_dt-CR1_dt)*dt
      cP=k22*cES_c*dt
      DeltaT_eq[i,j]<-t
    }
    #log10_cSf[i,j]<-log(cSf,10)
  }
}