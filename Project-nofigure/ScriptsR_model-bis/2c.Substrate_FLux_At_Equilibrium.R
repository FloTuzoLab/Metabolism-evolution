source("AP1.Plotting_multiple_images.R")
setwd("/Users/florian/Desktop/ScriptsR_modeles")
data<-read.table(file="Time_To_Equilibrium/Results/line0_timeToEq.rj",header=T)

#Determining the flow after 1sec. to quantify the effect of 


N_reso=41
#Creation of two sets of parameters
log10kf_var<-c(0,8,2)
log10kcat_var<-c(-4,4,2)
k2f<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
k2cat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
kr=10^3


E_tot_conc=10^-3
P1_Phi_eq=10^-6#c(10^-6,10^-3)

#Function determining the time to reach the equilibrium
#dt=0.00001 #timestep in seconds
Phi_eq<-function(t_max){
  Phi_ss<-matrix(nrow=N_reso , ncol=N_reso)
  for (i in 1:N_reso){
    for (j in 1:N_reso){
      print(c(i,j))
      dt<-min(1/(2*k2cat[j]),1/(2*kr),1/(5*k2f[i]*P1_Phi_eq),0.00001)
      cES=0
      cEf=E_tot_conc
      cSf=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
      ##Initalization of the process
      cP=0
      cSf_c<-cSf
      cEf_c<-cEf
      cES_c<-cES
      t=0
      while((t<t_max)){
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
        Phi_ss[i,j]<-cP
        t=t+dt
      }
    }
  }
  return(Phi_ss)
}

Phi_eq_CR<-list()
for (p in 1:length(P1_Phi_eq)){
  Phi_eq_CR[[p]]<-Phi_eq(1)
}

Phi_eq_CR[[1]]
image.plot(Phi_eq_CR[[1]],zlim=c(0,10^-11),axes=FALSE)
contour(Phi_eq_CR[[1]], add = TRUE,nlevels = 10)
axis(1, at=seq(0,1,0.25),labels=seq(0,8,2))
axis(2, at=seq(0,1,0.25),labels=seq(-4,4,2))





##Below that, no interest, scripts to determine an equilibrium that is trivial

#Defining features of plots
ncol=128
ncontour=5
pal<-colorRampPalette(c("green","green","red"))(ncol)


N_reso=50
DeltaT_eq<-matrix(nrow=N_reso , ncol=N_reso)
for (i in 1:N_reso){
  for (j in 1:N_reso){
    DeltaT_eq[i,j]<-data$T_eq[data$i==(i-1) & data$j==(j-1)]
  }
}

image.plot(DeltaT_eq,nlevel=ncol,col=pal,axes=F)

#First product concentration at steady-state
E_tot_conc=10^-3
P1_Phi_eq=c(10^-7,10^-5,10^-3)

##Analysis of the influence of reaction chemical parameters on equilibrium flux
N_reso=501
#Creation of two sets of parameters
log10k1_var<-c(0,10,2)
log10k2_var<-c(-4,6,2)
k1<-10^seq(0,10,length.out=N_reso)
k2<-10^seq(-4,6,length.out=N_reso)
lr<-c(0,1/100)

tab_S2_conc_eq <- matrix(nrow=N_reso , ncol=N_reso)
S2_conc_eq<-list()
log10_S2_conc_eq<-list()

for (l in 1:(length(lr))){
  for (s in 1:length(P1_Phi_eq)){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        k_1=lr[l]*k1[i]
        tab_S2_conc_eq[i,j]=P1_Phi_eq[[s]]*(k_1+k2[j])/(k1[i]*(k2[j]*E_tot_conc-P1_Phi_eq[[s]]))
      }
    }
    S2_conc_eq[[s+(l-1)*length(P1_Phi_eq)]]=tab_S2_conc_eq
    log10_S2_conc_eq[[s+(l-1)*length(P1_Phi_eq)]]=log(S2_conc_eq[[s+(l-1)*length(P1_Phi_eq)]],10)
  }
}

multiplePlot("k-1=","(k1)","[Se]=","M",lr,P1_Phi_eq,ncol,log10k1_var,log10k2_var,log10_S2_conc_eq,
             "Fitness landscape for the chemical reaction",
             "log10 k1","log10 k2",scale=c(-12,-2),ncont=ncontour)


##Modifying the way of representing results
tab_S2_conc_eq <- matrix(nrow=N_reso , ncol=N_reso)
S2_conc_eq<-list()
log10_S2_conc_eq<-list()

for (l in 1:(length(lr))){
  for (s in 1:length(P1_Phi_eq)){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        k_1=lr[l]*k1[i]
        Real_conc<-P1_Phi_eq[[s]]*(k_1+k2[j])/(k1[i]*(k2[j]*E_tot_conc-P1_Phi_eq[[s]]))
        if(Real_conc>1 || Real_conc<0){
          tab_S2_conc_eq[i,j]=1
        }
        else if(Real_conc<10^(-10)){
          tab_S2_conc_eq[i,j]=10^(-10)
        }
        else{
          tab_S2_conc_eq[i,j]=Real_conc
        }
      }
    }
    S2_conc_eq[[s+(l-1)*length(P1_Phi_eq)]]=tab_S2_conc_eq
    log10_S2_conc_eq[[s+(l-1)*length(P1_Phi_eq)]]=log(S2_conc_eq[[s+(l-1)*length(P1_Phi_eq)]],10)
  }
}

multiplePlot("k-1=","(k1)","[Se]=","M",lr,P1_Phi_eq,ncol,log10k1_var,log10k2_var,log10_S2_conc_eq,
             "Fitness landscape for the chemical reaction",
             "log10 k1","log10 k2",scale=c(-10,0),ncont=ncontour)
