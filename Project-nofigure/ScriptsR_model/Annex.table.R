setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R")
source("AP2.RaphsonNewton.R") 
require(Rmpfr)
library(RColorBrewer)
setwd(dir="Data")

###1.Simulation of the dynamics to test the analytical resolution accuracy for facilitated diffusion
#Particular conditions to test models
KT=5*10^-5 #Molar, Michaelis-like constant for facilitated diffusion
Vm=10^-6 #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
Se_set=10*KT
N_reso=3
tab_KM<-matrix(nrow=N_reso , ncol=N_reso)
log10k21_var<-c(5,7,2)
log10k22_var<-c(2,4,2)
#log10k1_var<-c(0,10,2)
#log10k2_var<-c(-4,6,2)
k11<-10^7
k12<-10^4
k21<-10^seq(log10k21_var[1],log10k21_var[2],length.out=N_reso)
k22<-10^seq(log10k22_var[1],log10k22_var[2],length.out=N_reso)
#k21<-10^seq(0,8,length.out=N_reso)
#k22<-10^seq(-3,5,length.out=N_reso)
E_tot_conc1=10^-3
E_tot_conc2=10^-3
k_1set<-c(10^3)
eta_d_set<-c(10^2)
#dt=0.00001 #timestep in seconds
#N_steps=2000000##will further need to be replace by a stoping condition of the algorithm, when concentrations nearby equilibrium
## Initializing current concentrations - representing the initial concentrations
cES1=0
cEf1=E_tot_conc1
cSf1=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
cP1=0
cES2=0
cEf2=E_tot_conc2
cSf2=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
##Initalization of the process
tab_P1_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P1_Phi_eq<-list()
tab_P2_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P2_Phi_eq<-list()
log10_P_Phi_eq<-list()
tab_last_it<-matrix(nrow=N_reso , ncol=N_reso)
last_it<-list()

for (l in 1:length(k_1set)){
  for (s in  1:length(Se_set)){
    Se=Se_set[s]
    eta_d<-eta_d_set[1]*Vm
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        cES1=0
        cEf1=E_tot_conc1
        cSf1=0
        cES2=0
        cEf2=E_tot_conc2
        cEP1=0
        cSf2=0
        cES1_c<-1
        cES2_c<-1
        t<-0
        k_1<-k_1set[l]
        #dt<-10^-5
        dt<-min(1/2*min(1/(k_1),1/(10*k21[i]*E_tot_conc2),1/(10*k11*E_tot_conc1),1/(10*k22[j]),1/(10*k12)),0.01)
        print(c(i,j,dt))
        diffdt_c<-1
        while((t*dt)<1){
          print(t*dt)
          diffdt_c<-cES2-cES2_c
          cSf1_c<-cSf1
          cEf1_c<-cEf1
          cES1_c<-cES1
          cSf2_c<-cSf2
          cEf2_c<-cEf2
          cEP1_c<-cEP1
          cES2_c<-cES2
          FD_dt=Vm*(Se-cSf1_c)/(KT+(Se+cSf1_c)+Se*cSf1_c/KT) ##Facilitated diffusion
          CR11_dt=k11*cSf1_c*cEf1_c-(k_1)*cES1_c
          CR12_dt=k12*cES1_c
          CR21_dt=k21[i]*cSf2_c*cEf2_c-(k_1)*cES2_c
          CR22_dt=k22[j]*cES2_c
          cSf1=cSf1_c+(FD_dt-CR11_dt)*dt 
          cEf1=cEf1_c+(CR12_dt-CR11_dt)*dt
          cES1=cES1_c-(CR12_dt-CR11_dt)*dt
          cSf2=cSf2_c+(CR12_dt-CR21_dt)*dt-eta_d*cSf2*dt
          cEf2=cEf2_c+(CR22_dt-CR21_dt)*dt
          cEP1=cEP1_c
          cES2=cES2_c-(CR22_dt-CR21_dt)*dt
          t=t+1
        }
        tab_P1_Phi_eq[i,j]=k12*cES1
        tab_P2_Phi_eq[i,j]=k22[j]*cES2
        tab_last_it[i,j]=(k22[j]*(cES2-cES2_c))/dt
        #print(cSf2)
      }
    }
    P1_Phi_eq[[s+length(Se_set)*(l-1)]]=tab_P1_Phi_eq
    P2_Phi_eq[[s+length(Se_set)*(l-1)]]=tab_P2_Phi_eq
    last_it[[s+length(Se_set)*(l-1)]]=tab_last_it
    #log10_P_Phi_eq[[s]]=log(P_Phi_eq[[s]],10)
  }
}




N_reso=3
#Defining parameters
Vm_set=c(10^-6)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#,10^-3)#
log10kf_var<-c(5,7,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(2,4,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-10^7
kcat1_set<-10^4
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-15


tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_background_eff<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_background_eff<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(eta_d_set))){
  for (s in 1:(length(kcat1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[1]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[1]
    Se=10*KT
    kf1<-kf1_set[s]
    kcat1<-kcat1_set[s]
    kr1<-10^3#kcat1
    for (i in 1:(N_reso)){
      init<-c(10^-20,10^-20)#c(10^-4,10^-10)
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
        tab_Conc_prod_eq[i,j]=SSR$y
        #init<-c(SSR$x,SSR$y)
      }
    }
    P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

for (i in 1:3){
  for (j in 1:3){
    print(paste(i,j))
    print((P2_Phi_eq[[1]][i,j]-P_Phi_eq_sens_2Enz_background_eff[[1]][i,j])/P_Phi_eq_sens_2Enz_background_eff[[1]][i,j])
  }
}

for (i in 1:3){
  for (j in 1:3){
    print(paste(i,j))
    print((P2_Phi_eq[[1]][i,j]-P_Phi_eq_sens_2Enz_background_eff[[1]][i,j]))
  }
}
