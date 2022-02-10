##Old scripts

###Older simulation results finding equilibrium with euler explicit methods

###1.Simulation of the dynamics to test the analytical resolution accuracy for facilitated diffusion
#Particular conditions to test models
#KT=0.1 #Molar, Michaelis-like constant for facilitated diffusion
#Vm=10^-3 #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
dP1max_dt=c(10^-5,10^-3)
K_s<-10^-2#c(10^-4,10^-2)
N_reso=20
tab_KM<-matrix(nrow=N_reso , ncol=N_reso)
log10k21_var<-c(3,8,1)
log10k22_var<-c(0,5,1)
k21<-10^seq(3,8,length.out=N_reso)
k22<-10^seq(0,5,length.out=N_reso)
k2_1set<-c(10^-2,10^4)
E2_tot_conc=10^-4
#dt=0.00001 #timestep in seconds
#N_steps=2000000##will further need to be replace by a stoping condition of the algorithm, when concentrations nearby equilibrium
## Initializing current concentrations - representing the initial concentrations
cES2=0
cEf2=E2_tot_conc
cPf1=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
##Initalization of the process
tab_P2_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P2_Phi_eq<-list()
tab_P1_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P1_Phi_eq<-list()
tab_P1_conc_eq<- matrix(nrow=N_reso , ncol=N_reso)
P1_conc_eq<-list()
log10_P1_conc_eq<-list()
tab_T_eq<-matrix(nrow=N_reso,ncol=N_reso)
T_eq<-list()

for (p in  1:length(dP1max_dt)){
  for (r in 1:length(k2_1set)){
    #Se=Se_set[s]
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        cES2=0
        cEf2=E2_tot_conc
        cPf1=0
        cES_c<-1
        t<-0
        k2_1<-k2_1set[r]
        #dt<-10^-5
        dt<-min(1/5*min(1/(k2_1),1/(k21[i]*E2_tot_conc),1/k22[j]),0.1)
        print(c(i,j,dt))
        #diffdt_c<-1
        if (k21[i] || k11[j])
          epsilon=0.01*dt
        while(abs((cES2-cES_c)/cES_c)>epsilon && t<100){
          #diffdt_c<-cES-cES_c
          #cSf_c<-cSf
          #cEf_c<-cEf
          if(cES2==0){
            cES_c<-1
          }
          else{
            cES_c<-cES2
          }
          FD_dt=dP1max_dt[p]*K_s/(K_s+cPf1) ##First reaction contribution
          CR21_dt=k21[i]*cPf1*cEf2-(k2_1)*cES2
          CR22_dt=k22[j]*cES2
          cPf1=cPf1+(FD_dt-CR21_dt)*dt 
          cEf2=cEf2+(CR22_dt-CR21_dt)*dt
          cES2=cES2-(CR22_dt-CR21_dt)*dt
          t=t+dt
        }
        tab_P1_Phi_eq[i,j]=FD_dt
        tab_P2_Phi_eq[i,j]=k22[j]*cES2
        tab_P1_conc_eq[i,j]=cPf1
        tab_T_eq[i,j]=t
      }
    }
    P1_Phi_eq[[p+length(dP1max_dt)*(r-1)]]=tab_P1_Phi_eq
    P2_Phi_eq[[p+length(dP1max_dt)*(r-1)]]=tab_P2_Phi_eq
    P1_conc_eq[[p+length(dP1max_dt)*(r-1)]]=tab_P1_conc_eq
    T_eq[[p+length(dP1max_dt)*(r-1)]]=tab_T_eq
    log10_P1_conc_eq[[p+length(dP1max_dt)*(r-1)]]=log(P1_conc_eq[[p+length(dP1max_dt)*(r-1)]],10)
  }
}

wP1_2r_k1<-list()
wP2_2r_k1<-list()
for (r in 1:(length(k2_1set))){
  for (p in 1:(length(dP1max_dt))){
    wP1_2r_k1[[p+length(dP1max_dt)*(r-1)]]<-round(P1_Phi_eq[[p+length(dP1max_dt)*(r-1)]]/max(P1_Phi_eq[[p+length(dP1max_dt)*(r-1)]]),
                                                  digits=8)
    wP2_2r_k1[[p+length(dP1max_dt)*(r-1)]]<-round(P2_Phi_eq[[p+length(dP1max_dt)*(r-1)]]/max(P2_Phi_eq[[p+length(dP1max_dt)*(r-1)]]),
                                                  digits=8)
  }
}

multiplePlot("k-1=","","d[P1]/dt=","M/s",k2_1set,dP1max_dt,ncol,log10k21_var,log10k22_var,P1_Phi_eq,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2",lev=c(0.1,0.95,0.99,0.999,0.999999))

multiplePlot("k-1=","","d[P1]/dt=","M/s",k2_1set,dP1max_dt,ncol,log10k21_var,log10k22_var,P2_Phi_eq,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2",lev=c(0.1,0.95,0.99,0.999))

multiplePlot("k-1=","","d[P1]/dt=","M/s",k2_1set,dP1max_dt,ncol,log10k21_var,log10k22_var,log10_P1_conc_eq,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2",lev=c(0.1,0.95,0.99,0.999))

multiplePlot("k-1=","","d[P1]/dt=","M/s",k2_1set,dP1max_dt,ncol,log10k21_var,log10k22_var,T_eq,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2")

multiplePlot("k-1=","","d[P1]/dt=","M/s",k2_1set,dP1max_dt,ncol,log10k21_var,log10k22_var,wP1_2r_k1,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2",scale=c(0,1),lev=c(0.9,0.99,0.999,0.9999,0.99999))

multiplePlot("k-1=","","d[P1]/dt=","M/s",k2_1set,dP1max_dt,ncol,log10k21_var,log10k22_var,wP2_2r_k1,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2",scale=c(0,1),lev=c(0.1,0.95,0.99,0.999,0.9999,0.99999))

###1.Simulation of the dynamics to test the analytical resolution accuracy for facilitated diffusion
#Particular conditions to test models
KT=0.001 #Molar, Michaelis-like constant for facilitated diffusion
Vm=10^-3 #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
Se_set=10^-3
N_reso=10
tab_KM<-matrix(nrow=N_reso , ncol=N_reso)
log10k21_var<-c(4,8,2)
log10k22_var<-c(1,5,2)
#log10k1_var<-c(0,10,2)
#log10k2_var<-c(-4,6,2)
k11<-10^6
k12<-10^2
k21<-10^seq(log10k21_var[1],log10k21_var[2],length.out=N_reso)
k22<-10^seq(log10k22_var[1],log10k22_var[2],length.out=N_reso)
#k21<-10^seq(0,8,length.out=N_reso)
#k22<-10^seq(-3,5,length.out=N_reso)
E_tot_conc1=10^-4
E_tot_conc2=10^-7
k_1set<-c(10^2)

kin2<-1#A TESTER AVEC 0
kout2<-10^2#A TESTER AVEC 0
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

for (l in 1:length(k_1set)){
  for (s in  1:length(Se_set)){
    Se=Se_set[s]
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
        while(abs((cES2-cES2_c)/diffdt_c)>0.99 || t<10){
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
          CR1bis_dt=kout2*cSf2_c*cEf1_c-kin2*cEP1_c
          CR21_dt=k21[i]*cSf2_c*cEf2_c-(k_1)*cES2_c
          CR22_dt=k22[j]*cES2_c
          cSf1=cSf1_c+(FD_dt-CR11_dt)*dt 
          cEf1=cEf1_c+(CR12_dt-CR11_dt-CR1bis_dt)*dt
          cES1=cES1_c-(CR12_dt-CR11_dt)*dt
          cSf2=cSf2_c+(CR12_dt-CR21_dt)*dt 
          cEf2=cEf2_c+(CR22_dt-CR21_dt)*dt
          cEP1=cEP1_c+CR1bis_dt*dt
          cES2=cES2_c-(CR22_dt-CR21_dt)*dt
          t=t+1
        }
        tab_P1_Phi_eq[i,j]=k12*cES1
        tab_P2_Phi_eq[i,j]=k22[j]*cES2
        #print(cSf2)
      }
    }
    P1_Phi_eq[[s+length(Se_set)*(l-1)]]=tab_P1_Phi_eq
    P2_Phi_eq[[s+length(Se_set)*(l-1)]]=tab_P2_Phi_eq
    #log10_P_Phi_eq[[s]]=log(P_Phi_eq[[s]],10)
  }
}

cP=k2*cES
plot.new()
par(mfrow=c(1,2))
image.plot(P1_Phi_eq[[1]])
image.plot(P1_Phi_eq[[2]])

image.plot(P2_Phi_eq[[1]])
image.plot(P2_Phi_eq[[2]])

wP1_2r_k1<-list()
wP2_2r_k1<-list()
for (l in 1:(length(k_1set))){
  for (s in 1:(length(Se_set))){
    wP1_2r_k1[[s+length(Se_set)*(l-1)]]<-round(P1_Phi_eq[[s+length(Se_set)*(l-1)]]/max(P1_Phi_eq[[s+length(Se_set)*(l-1)]]),digits=3)
    wP2_2r_k1[[s+length(Se_set)*(l-1)]]<-round(P2_Phi_eq[[s+length(Se_set)*(l-1)]]/max(P2_Phi_eq[[s+length(Se_set)*(l-1)]]),digits=3)
  }
}

multiplePlot("k-1=","","[Se]=","M",k_1set,Se_set,ncol,log10k21_var,log10k22_var,wP1_2r_k1,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2",scale=c(0,1),lev=c(0.1,0.95,0.99,0.999))

multiplePlot("k-1=","","[Se]=","M",k_1set,Se_set,ncol,log10k21_var,log10k22_var,wP2_2r_k1,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 k1","log10 k2",scale=c(0,1),lev=c(0.1,0.95,0.99,0.999))



#Testing the equilibrium
init<-c(10^-8,10^-20)
kr1=10^3
kr2=10^-2
Vm=10^-3
kf1<-10^10
kf2<-10^2
kcat1<-10^5
kcat2<-10^-4
A1<-kcat1*kf1*E1_conc
A2<-kcat1+kr1
Ki<-10^5
A3<-kcat2*kf2*E2_conc
A4<-kcat2+kr2
B1<-KT+Se
B2<-1+Se/KT
#Steady-state equations F(x,y)=0 to solve 
F1e<-expression(A1*x/((1+Ki*y)*(A2+kf1*x))-A3*y/(A4+kf2*y))
F2e<-expression(Vm*(Se-x)/(B1+B2*x)-A3*y/(A4+kf2*y))
Flist<-list(F1e,F2e)
SSR<-RaphsonNewton(Flist,init,error)
print(paste(SSR$x,SSR$y))
print(A3*SSR$y/(A4+kf2*SSR$y))
##Detemination de chacune des concentrations pour déceler erreur éventuelle
E1f=E1_conc/((1+Ki*SSR$y)*(1+kf1*SSR$x/(kcat1+kr1)))
E1Sf=kf1*SSR$x*E1_conc/((1+Ki*SSR$y)*((kcat1+kr1)+kf1*SSR$x))
dP1=kcat1*E1Sf
E2f=E2_conc/(1+kf2*SSR$y/(kcat2+kr2))
E2P1f=kf2*SSR$y/(kcat2+kr2)*E2f
dP2=kcat2*E2P1f

E1f
E1Sf
E2f
E2P1f
dP1
dP2

Fxy(F1e,SSR$x,SSR$y)
Fxy(F2e,SSR$x,SSR$y)

library(lattice)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(plot.matrix)

#TESTING inhibition

Vmax1=10^-3*10^5
Vmax2=10^-3*10^1
KM1=(10^5+10^3)/10^9
KM2=(10^0+10^-3)/10^8
KI=10^-7

delta=(Vmax1-Vmax2)^2+4*Vmax1*Vmax2*KM2/KI
P1=(Vmax1-Vmax2+sqrt(delta))/(2*Vmax2/KI)

