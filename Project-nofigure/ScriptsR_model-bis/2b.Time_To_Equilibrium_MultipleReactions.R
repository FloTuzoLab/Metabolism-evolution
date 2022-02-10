source("AP1.Plotting_multiple_images.R")
require("Rmpfr")
#Defining features of plots
ncol=128
ncontour=5
jet.colors <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
palet<-list(jet.colors(ncol))

#Time to reach the equilibrium - Chemical reaction only - using C++ to limit numercial problems

###1.Simulation of the dynamics to test the analytical resolution accuracy for facilitated diffusion+first reaction

##Simulation procedure and results
#Particular conditions to test models
KT=0.1 #Molar, Michaelis-like constant for facilitated diffusion
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
Se=1

N_reso=50
log10kf_var<-c(0,10,2)
log10kcat_var<-c(-4,6,2)
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
kr=10^3
E_tot_conc<-10^3

#Function determining the time to reach the equilibrium
#dt=0.00001 #timestep in seconds
Phi_eq_FD<-function(Vm,t_max){
  Phi_ss<-matrix(nrow=N_reso , ncol=N_reso)
  for (i in 1:N_reso){
    for (j in 1:N_reso){
      print(c(i,j))
      dt<-min(1/(2*kcat[j]),1/(2*kr),1/(E_tot_conc*kf[i]),0.0001)
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
        FD_dt=Vm*(Se-cSf_c)/(KT+(Se+cSf_c)+Se*cSf_c/KT) ##Facilitated diffusion
        CR1_dt=kf[i]*cSf_c*cEf_c-(kr)*cES_c
        CR2_dt=kcat[j]*cES_c
        cSf=cSf_c+(FD_dt-CR1_dt)*dt 
        cEf=cEf_c+(CR2_dt-CR1_dt)*dt
        cES=cES_c-(CR2_dt-CR1_dt)*dt
        cP=kcat[j]*cES_c*dt
        Phi_ss[i,j]<-cP/dt
        t=t+dt
      }
    }
  }
  return(Phi_ss)
}

Phi_eq_T1<-list()
Phi_eq_T01<-list()
tm<-0.1#Prev:1
for (p in 1:length(Vm_set)){
  Vm<-Vm_set[p]
  Phi_eq_T01[[p]]<-Phi_eq_FD(Vm,tm)
}

setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles/Equilibrium_time/")
#save(Phi_eq_T1,file="Phi_eq_T1.Rdata")
#save(Phi_eq_T01,file="Phi_eq_T01.Rdata")
load("/Users/florian/Desktop/ScriptsR_modeles/Equilibrium_time/Phi_eq_T1.Rdata")
load("/Users/florian/Desktop/ScriptsR_modeles/Equilibrium_time/Phi_eq_T01.Rdata")

multiplePlot("kr=","","Vm=","M",kr,Vm,ncol,log10kf_var,log10kcat_var,Phi_eq_T01,
             "Absolute fitness landscape for the system
{passive diffusion + chemical reaction}",
             "log10 kf","log10 kcat",ncont=10,palette=palet,image=TRUE,scale="AUTO") ##The unit and the value of Vm need to be checked and explained


w_kf_TE<-list()
for (l in 1:(length(kr))){
  for (p in 1:length(Vm_set)){
    w_kf_TE[[p+length(kr)*(l-1)]]<-round(Phi_eq_T1[[p+length(kr)*(l-1)]]/max(Phi_eq_T1[[p+length(kr)*(l-1)]]),digits=10)
  }
} 

multiplePlot("kr=","","Vm=","M",kr,Vm,ncol,log10kf_var,log10kcat_var,w_kf_TE,
             "Absolute fitness landscape for the system
{passive diffusion + chemical reaction}",
         "log10 kf","log10 kcat",palette=palet,image=TRUE,lev=c(0.9,0.999999),scale="AUTO") ##The unit and the value of Vm need to be checked and explained

###2.Simulation of the dynamics to test the analytical resolution accuracy for facilitated diffusion+two reactions 
#when both reactions have the same kinetic parameters (most defavourable case)

##Simulation procedure and results
#Particular conditions to test models
KT=0.1 #Molar, Michaelis-like constant for facilitated diffusion
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
Se=1
N_reso=50

log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf_set<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat_set<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
kr_set<-10^3#c(10^-1,10^5) #An intermediate value seems reasonable to show general results
E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]

#Function determining the time to reach the equilibrium
#dt=0.00001 #timestep in seconds
Phi_eq_T2E<-function(Vm,t_max){
  Phi2_ss<-matrix(nrow=N_reso , ncol=N_reso)
  for (i in 1:N_reso){
    for (j in 1:N_reso){
      print(c(i,j))
      kr=kr_set[1]
      #Defining intermediate parameters
      kf<-kf_set[i]
      kcat<-kcat_set[j]
      cES1=0
      cEf1=E_tot_conc
      cSf1=0
      cES2=0
      cEf2=E_tot_conc
      #cEP1=0
      cSf2=0
      #cES1_c<-0
      #cES2_c<-0
      dt<-min(1/(2*kcat),1/(2*kr),1/(E_tot_conc*kf),0.0001)
      t=0
      while((t<t_max)){
        #print(t)
        cSf1_c<-cSf1
        cEf1_c<-cEf1
        cES1_c<-cES1
        cSf2_c<-cSf2
        cEf2_c<-cEf2
        #cEP1_c<-cEP1
        cES2_c<-cES2
        FD_dt=Vm*(Se-cSf1_c)/(KT+(Se+cSf1_c)+Se*cSf1_c/KT) ##Facilitated diffusion
        #print(FD_dt)
        CR11_dt=kf*cSf1_c*cEf1_c-(kr)*cES1_c
        CR12_dt=kcat*cES1_c
        #print(kcat)
        #CR1bis_dt=kout2*cSf2_c*cEf1_c-kin2*cEP1_c
        CR21_dt=kf*cSf2_c*cEf2_c-kr*cES2_c
        
        CR22_dt=kcat*cES2_c
        cSf1=cSf1_c+(FD_dt-CR11_dt)*dt 
        cEf1=cEf1_c+(CR12_dt-CR11_dt)*dt#-CR1bis_dt)*dt
        cES1=cES1_c-(CR12_dt-CR11_dt)*dt
        cSf2=cSf2_c+(CR12_dt-CR21_dt)*dt 
        cEf2=cEf2_c+(CR22_dt-CR21_dt)*dt
        #print(cEf2)
        #cEP1=cEP1_c+CR1bis_dt*dt
        cES2=cES2_c-(CR22_dt-CR21_dt)*dt
        Phi2_ss[i,j]<-CR22_dt
        t=t+dt
      }
    }
  }
  return(Phi2_ss)
}

Phi2_eq_T1<-list()
Phi2_eq_T01<-list()
Phi2_eq_T02<-list()
tm<-0.2#Prev:1
for (p in 1:length(Vm_set)){
  Vm<-Vm_set[p]
  Phi2_eq_T02[[p]]<-Phi_eq_T2E(Vm,tm)
}

setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles/Equilibrium_time/")
save(Phi2_eq_T1,file="Phi2_eq_T1.Rdata")
save(Phi2_eq_T02,file="Phi2_eq_T02.Rdata")
save(Phi2_eq_T01,file="Phi2_eq_T01.Rdata")


multiplePlot("kr=","","Vm=","M",kr,Vm,ncol,log10kf_var,log10kcat_var,Phi2_eq_T01,
             "Absolute fitness landscape for the system
{passive diffusion + chemical reaction}",
             "log10 kf","log10 kcat",ncont=10,palette=palet,image=TRUE,scale="AUTO") ##The unit and the value of Vm need to be checked and explained

multiplePlot("kr=","","Vm=","M",kr,Vm,ncol,log10kf_var,log10kcat_var,Phi2_eq_T02,
             "Absolute fitness landscape for the system
{passive diffusion + chemical reaction}",
             "log10 kf","log10 kcat",ncont=10,palette=palet,image=TRUE,scale="AUTO") ##The unit and the value of Vm need to be checked and explained

multiplePlot("kr=","","Vm=","M",kr,Vm,ncol,log10kf_var,log10kcat_var,Phi2_eq_T1,
             "Absolute fitness landscape for the system
{passive diffusion + chemical reaction}",
             "log10 kf","log10 kcat",ncont=10,palette=palet,image=TRUE,scale="AUTO") ##The unit and the value of Vm need to be checked and explained
