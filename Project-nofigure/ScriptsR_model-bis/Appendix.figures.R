setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R")
source("AP2.RaphsonNewton.R") 
require(Rmpfr)
library(RColorBrewer)
setwd(dir="Data")

options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=128
ncontour=5
#palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato"))
palet<-jet.colors(ncol)
pal<-list(palet)

###0.Drawing of fitness landscapes for passive diffusion

r_cell=1*10^-6#unit:m
SA<-4*pi*r_cell^2
V<-4/3*pi*r_cell^3

Se_high=0.3#very high concentration
Se_low=10^-3#very high concentration
Se_set<-c(Se_high,Se_low)
N_reso=100
Pd_set=c(10^-12)#unit:m/s true glucose value:1e-12
E_tot_conc=10^-3#unit:M
log10kf_var<-c(0,10,2)
log10kcat_var<-c(-4,6,2)
kr<-10^3#intermediate level for reverse reaction constant
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
tab_P_Phi_eq_high <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_high<-list()

for(s in 1:length(Se_set)){
  Pd<-Pd_set[1]
  Se<-Se_set[s]
  for (i in 1:N_reso){
    for (j in 1:N_reso){
      a<-mpfr(-Pd*SA/V*kf[i],120)
      b<-mpfr(Pd*SA/V*(Se*kf[i]-(kr+kcat[j]))-kcat[j]*kf[i]*E_tot_conc,120)
      c<-mpfr(Pd*SA/V*Se*(kr+kcat[j]),120)
      delta_pd<-mpfr(b^2-4*a*c,120)
      S_conc_eq=mpfr((-b-delta_pd^(1/2))/(2*a),120)
      ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
      tab_P_Phi_eq_high[i,j]=as.numeric(kcat[j]*ES_conc_eq)
    }
  }
  P_Phi_eq_high[[s]]<-tab_P_Phi_eq_high
}

addtxt<-list(l=0.98,h=0.97,txt=c("A","B","C"),srt = 0,font=2,col="black")
multiplePlot("kr=","","P=","M/s",kr,Se_set[1],ncol,log10kf_var,log10kcat_var,list(P_Phi_eq_high[[1]]),
             abs="log10 kf",ord="log10 kcat",palette=pal,lev=c(0.9*max(P_Phi_eq_high[[1]])),image=TRUE,scale="AUTO",subcex=1,pcex=1,labcex=1,axcex=1.25,meth="flattest",legcex=0.4) ##The unit and the value of Vm need to be checked and explained
text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)


multiplePlot("kr=","","P=","M/s",kr,Se_set[2],ncol,log10kf_var,log10kcat_var,list(P_Phi_eq_high[[2]]),
             abs="log10 kf",ord="log10 kcat",palette=pal,lev=c(0.9*max(P_Phi_eq_high[[2]])),image=TRUE,scale="AUTO",subcex=1,pcex=1,labcex=1,axcex=1.25,meth="flattest",legcex=0.4) ##The unit and the value of Vm need to be checked and explained

text(addtxt$l,addtxt$h,addtxt$txt[2],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)





#I Multiple enzyme pathways

##4.Best value for the first enzyme followed by a variable second enzyme with sensitivy study
##With degradation in the pathway
N_reso=10
#Defining parameters
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1<-10^10
kcat1<-10^6
kr1<-10^3
#Ki<-10^-6
eta_d_set<-c(10^2,10^4)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(KT_set))){
  for (s in 1:(length(eta_d_set))){
    #kr=kr_set[1]
    Vm=Vm_set[s]
    eta_d<-eta_d_set[s]*Vm
    KT=KT_set[l]
    Se=10*KT
    for (i in 1:(N_reso)){
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
      }
    }
    P_Phi_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}


#b.Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Vm_set))){
    w_kf_FD_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]),digits=20)
  }
}

#image.plot(P_Phi_eq[[1]],legend.shrink=1)
options(digits=3)
subtitle<-list("(a) Weak flux, high affinity",
               "(b) Moderate flux, high affinity",
               "(c) High flux, high affinity",
               "(d) Weak flux, moderate affinity",
               "(e) Moderate flux, moderate affinity",
               "(f) High flux, moderate affinity",
               "(g) Weak flux, low affinity",
               "(h) Moderate flux, low affinity",
               "(i) High flux, low affinity")
multiplePlot("KT=","","Vm=","M/s",KT_set,signif(Vm_set,4),ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle)

multiplePlot("KT=","","Vm=","M/s",KT_set,signif(Vm_set,4),ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle)





##1.Same value for both enzymes
##With degradation in the pathway
N_reso=100
#Defining parameters
Vm_set=c(10^-6,10^(-4.5),10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(10^-5,10^-3.5,10^-2)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

kf_set<-kf
kcat_set<-kcat
kr_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Vm_set))){
    #kr=kr_set[1]
    Vm=Vm_set[s]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[l]
    Se=10*KT
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        kf1<-kf_set[i]
        kf2<-kf_set[i]
        kcat1<-kcat_set[j]
        kcat2<-kcat_set[j]
        kr1<-kr_set[1]
        kr2<-kr_set[1]
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
      }
    }
    P_Phi_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}


#b.Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Vm_set))){
    w_kf_FD_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]),digits=20)
  }
}

#image.plot(P_Phi_eq[[1]],legend.shrink=1)
options(digits=3)
subtitle<-list("(a) Weak flux, high affinity",
               "(b) Moderate flux, high affinity",
               "(c) High flux, high affinity",
               "(d) Weak flux, moderate affinity",
               "(e) Moderate flux, moderate affinity",
               "(f) High flux, moderate affinity",
               "(g) Weak flux, low affinity",
               "(h) Moderate flux, low affinity",
               "(i) High flux, low affinity")
multiplePlot("KT=","","Vm=","M/s",signif(KT_set,3),signif(Vm_set,3),ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle)

#Results similar than that with only one enzyme - at least qualitatively

##2.Best value for the first enzyme followed by a variable second enzyme with sensitivy study
##With degradation in the pathway
N_reso=20
#Defining parameters
log10Vm_var<-c(-10,0,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10KT_var<-c(-10,0,2) #catalytic preferred to -1 for same reasons
Vm_set<-10^seq(log10Vm_var[1],log10Vm_var[2],length.out=N_reso)
KT_set<-10^seq(log10KT_var[1],log10KT_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf_set<-c(10^2,10^6,10^10)
kcat_set<-c(10^-2,10^2,10^6)
kf_rat<-c(10^-4,10^-2,1)
kr<-10^3
kr1<-kr
kr2<-kr
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_FDlandscape<-list() #Equilibrium product flux
Conc_prod_eq_FDlandscape<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(kf_set))){
  #for (s in 1:(length(kf_rat))){
  for (s in 1:(length(kcat_set))){
    kf1<-kf_set[l]
    kf2<-kf1
    #kf2<-kf1*kf_rat[s]
    kcat1<-kcat_set[s]
    kcat2<-kcat1
    #Kcat1<-10^4
    #kcat2<-10^1
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        print(paste(i,j))
        Vm=Vm_set[i]
        eta_d<-eta_d_set[1]*Vm
        KT=KT_set[j]
        Se=10*KT
        
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
      }
    }
    P_Phi_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}


#b.Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_FDlandscape<-list()
for (l in 1:(length(kf_set))){
  for (s in 1:(length(kcat_set))){
    w_kf_FD_FDlandscape[[s+length(kcat_set)*(l-1)]]<-round(P_Phi_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]/max(P_Phi_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]),digits=20)
  }
}

#image.plot(P_Phi_eq[[1]],legend.shrink=1)
options(digits=3)
subtitle<-list("(a) Weak flux, high affinity",
               "(b) Moderate flux, high affinity",
               "(c) High flux, high affinity",
               "(d) Weak flux, moderate affinity",
               "(e) Moderate flux, moderate affinity",
               "(f) High flux, moderate affinity",
               "(g) Weak flux, low affinity",
               "(h) Moderate flux, low affinity",
               "(i) High flux, low affinity")
multiplePlot("kf=","/M/s","kcat=","/s",kf_set,kcat_set,ncol,log10Vm_var,log10KT_var,w_kf_FD_FDlandscape,
             abs="log10 (Vm)",ord="log10 (KT)",lev=c(0.01,0.1,0.5,0.9,0.99),palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle)

##2.Best value for the first enzyme followed by a variable second enzyme with sensitivy study
##With degradation in the pathway
N_reso=20
#Defining parameters
log10Vm_var<-c(-10,0,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10KT_var<-c(-10,0,2) #catalytic preferred to -1 for same reasons
Vm_set<-10^seq(log10Vm_var[1],log10Vm_var[2],length.out=N_reso)
KT_set<-10^seq(log10KT_var[1],log10KT_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^2,10^6,10^10)
kcat_set<-c(10^-2,10^2,10^6)
kf2_set<-c(10^2,10^6,10^10)
kr<-10^3
kr1<-kr
kr2<-kr
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_FDlandscape<-list() #Equilibrium product flux
Conc_prod_eq_FDlandscape<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(kf1_set))){
  for (s in 1:(length(kf2_set))){
  #for (s in 1:(length(kcat_set))){
    kf1<-kf1_set[l]
    #kf2<-kf1
    kf2<-kf2_set[s]
    kcat1<-kcat_set[2]
    kcat2<-kcat1
    #Kcat1<-10^4
    #kcat2<-10^1
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        print(paste(i,j))
        Vm=Vm_set[i]
        eta_d<-eta_d_set[1]*Vm
        KT=KT_set[j]
        Se=10*KT
        
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
      }
    }
    P_Phi_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}


#b.Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_FDlandscape<-list()
for (l in 1:(length(kf_set))){
  for (s in 1:(length(kcat_set))){
    w_kf_FD_FDlandscape[[s+length(kcat_set)*(l-1)]]<-round(P_Phi_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]/max(P_Phi_eq_FDlandscape[[s+length(kcat_set)*(l-1)]]),digits=20)
  }
}

#image.plot(P_Phi_eq[[1]],legend.shrink=1)
options(digits=3)
subtitle<-list("(a) Bad 1st enzyme, Bad 2nd enzyme",
               "(b) Bad 1st enzyme, Moderate 2nd enzyme",
               "(c) Bad 1st enzyme, Good 2nd enzyme",
               "(d) Moderate 1st enzyme, Bad 2nd enzyme",
               "(e) Moderate 1st enzyme, Moderate 2nd enzyme",
               "(f) Moderate 1st enzyme, Good 2nd enzyme",
               "(g) Good 1st enzyme, Bad 2nd enzyme",
               "(h) Good 1st enzyme, Moderate 2nd enzyme",
               "(i) Good 1st enzyme, Good 2nd enzyme")
multiplePlot("kf1=","/M/s","kf2=","/s",kf1_set,kf2_set,ncol,log10Vm_var,log10KT_var,w_kf_FD_FDlandscape,
             abs="log10 (Vm)",ord="log10 (KT)",lev=c(0.01,0.1,0.5,0.9,0.99),palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle)

multiplePlot("kf=","/M/s","kf2/kf1=","/s",kf_set,kf_rat,ncol,log10Vm_var,log10KT_var,Conc_prod_eq_FDlandscape,
             abs="log10 (Vm)",ord="log10 (KT)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle)

###6.Effect of background enzyme on the landscape of the second enzyme
N_reso=10
#Defining parameters
Vm_set=c(10^-6)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#,10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^0,10^9)
kcat1_set<-c(10^-4,10^5)
Keq_set<-c(10^-3,10^2)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

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
for (l in 1:(length(Keq_set))){
  for (s in 1:(length(kcat1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[1]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[1]
    Se=10*KT
    kf1<-kf1_set[s]
    kcat1<-kcat1_set[s]
    kr1<-kcat1*sqrt(Keq_set[l])
    kinh1<-kf1*sqrt(Keq_set[l])
    for (i in 1:(N_reso)){
      init<-c(10^-20,10^-20)
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1_pos<-kcat1*E1_conc
        KM1_pos<-(kcat1+kr1)/kf1
        K_I<-kinh1/kf1
        Vm1_neg<-kr1*E1_conc
        KM1_neg<-(kcat1+kr1)/kinh1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1_pos*x/(KM1_pos+x+K_I*y)-Vm1_neg*y/(KM1_neg+y+x/K_I)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y)+eta_d)*y)
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

#Analysing the results
##Maximum of fitness
max(P_Phi_eq_sens_2Enz_background_eff[[1]])/max(P_Phi_eq_sens_2Enz_background_eff[[2]])
multiplePlot("kf=","/M/s","kf2/kf1=","/s",Keq_set,kcat1_set,ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz_background_eff,
             abs="log10 (Vm)",ord="log10 (KT)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle,subcex=1)


##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff<-list()
for (l in 1:(length(Keq_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}

multiplePlot("kf=","/M/s","kf2/kf1=","/s",Keq_set,kcat1_set,ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz_background_eff,
             abs="log10 (Vm)",ord="log10 (KT)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,sub=subtitle,subcex=1)

plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1, 0.75, 1))

contour(w_kf_FD_sens_2Enz_background_eff[[4]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1))
title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 (kf)",ylab="log10 (kcat)",cex.lab=1.25)
  for (s in 1:3){
    contour(w_kf_FD_sens_2Enz_background_eff[[s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=s+1,lty=1,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1.25)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1.25)
    print((l-1)*length(kcat1_set)+s)
}
text(addtxt$l,addtxt$h,addtxt$txt[2],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1.5)
legend("bottomleft",title=c("Equilibrium constant"),legend=c("Very inefficient first enzyme","Inefficient first enzyme","Perfect first enzyme"),lty=1,col=c(2,3,1),ncol=1,cex=0.7)
title(main=bquote(atop("Low degradation rate - "~ eta~" ="~100~"/s,",
                       "reversible first reaction" ~ Keq~" ="~1))
      ,col.main="black",font.main=1,cex.main=1.1)




##5.Best value for the first enzyme followed by a variable second enzyme with reversibility
##And with degradation in the pathway
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)
###With high affinity, testing the effect of reversibility
N_reso=10
#Defining parameters
Vm_set=c(10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^2,10^9)
kcat1_set<-c(10^-2,10^5)
Keq_set<-c(10^-2,1,10^2)
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_rev<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_rev<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(Keq_set))){
  for (s in 1:(length(kf1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[l]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[l]
    Se=10*KT
    kf1<-kf1_set[s]
    kcat1<-kcat1_set[s]
    kr1<-sqrt(Keq_set[l])*kcat1
    kinh1<-sqrt(Keq_set[l])*kf1
    for (i in 1:(N_reso)){
      init<-c(10^-4,10^-10)
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1_pos<-kcat1*E1_conc
        KM1_pos<-(kcat1+kr1)/kf1
        K_I<-kinh1/kf1
        Vm1_neg<-kr1*E1_conc
        KM1_neg<-(kcat1+kr1)/kinh1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1_pos*x/(KM1_pos+x+K_I*y)-Vm1_neg*y/(KM1_neg+y+x/K_I)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y)+eta_d)*y)
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
        tab_Conc_prod_eq[i,j]=SSR$y
      }
    }
    P_Phi_eq_sens_2Enz_rev[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_rev[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_rev<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(kr1_set))){
    w_kf_FD_sens_2Enz_rev[[s+length(kr1_set)*(l-1)]]<-round(P_Phi_eq_sens_2Enz_rev[[s+length(kr1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_rev[[s+length(kr1_set)*(l-1)]]),digits=20)
  }
}

options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=128
ncontour=5
jet.colors <- colorRampPalette(c("blue", "cyan", "yellow", "red")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

#multiplePlot("KT=","","Vm=","M/s",signif(Vm_set,4),kr1_set,ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz_rev,
#abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,sub=subtitle,subcex=2)

plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(1, 1, 1, 1))
contour(w_kf_FD_sens_2Enz_rev[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1.5,col=1,lty=1,method ="edge",xlim=c(0,1))
title(main="",outer=FALSE,line=3,cex.main=1.5,font=2,xlab="log10 (kf)",ylab="log10 (kcat)",cex.lab=1.5)
for (l in 1:(length(Vm_set))){
  for (s in 2:length(kr1_set)){
    contour(w_kf_FD_sens_2Enz_rev[[(l-1)*length(kr1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1.5,col=s,lty=l+0,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1.5)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1.5)
    print((l-1)*length(Vm_set)+s)
  }
}

legend("bottomleft",title=c("Equilibrium constant"),legend=c(expression("Keq="~10^{-4}),expression("Keq="~10^{-2}),expression("Keq="~1),expression("Keq="~10^{2})),lty=1,col=c(3,4,1,2),ncol=2,cex=1,
       main="(c) Influence of ")

#coordonnées des fleches 1
x0<-c(0.2)
y0<-c(0.2)
x1<-c(0.8)
y1<-c(0.8)
arrows(x0,y0,x1,y1,code=2,lwd=2,lty=4)
mtext(text="Effect of increasing
      reversibility",at=c(0.8),line=-5,outer=FALSE,font=2)


##6b.Best value for the first enzyme followed by a variable second enzyme with reversibility: testing the effect of reversibility on the first enzyme ##Useless as results are not very influential
##And with degradation in the pathway

###With high affinity, testing the effect of reversibility
N_reso=10
#Defining parameters
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf2<-10^10
kcat2<-10^6
Keq_set<-c(1,10^-4,10^-2,10^2)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

kf1_set<-kf
kcat1_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_rev<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_rev<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set))){
    #kr=kr_set[1]
    Vm=Vm_set[l]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[l]
    Se=10*KT
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf1<-kf1_set[i]
        #print(paste(i,"kf2",kf2))
        kcat1<-kcat1_set[j]
        kr1<-kcat1*sqrt(Keq_set[s])
        kinh1<-kf1*sqrt(Keq_set[s])
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1_pos<-kcat1*E1_conc
        KM1_pos<-(kcat1+kr1)/kf1
        K_I<-kinh1/kf1
        Vm1_neg<-kr1*E1_conc
        KM1_neg<-(kcat1+kr1)/kinh1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1_pos*x/(KM1_pos+x+K_I*y)-Vm1_neg*y/(KM1_neg+y+x/K_I)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y)+eta_d)*y)
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
        tab_Conc_prod_eq[i,j]=SSR$y
      }
    }
    P_Phi_eq_sens_2Enz_rev[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_rev[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_rev<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set))){
    w_kf_FD_sens_2Enz_rev[[s+length(kr1_set)*(l-1)]]<-round(P_Phi_eq_sens_2Enz_rev[[s+length(kr1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_rev[[s+length(kr1_set)*(l-1)]]),digits=20)
  }
}

options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=128
ncontour=5
palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("steelblue", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

#multiplePlot("KT=","","Vm=","M/s",signif(Vm_set,4),kr1_set,ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz_rev,
#abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,sub=subtitle,subcex=2)

addtxt<-list(l=0.95,h=0.95,txt=c("C","D"),srt = 0,font=2,col="black")

contour(w_kf_FD_sens_2Enz_rev[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1))
title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 (kf)",ylab="log10 (kcat)",cex.lab=1.25)
for (l in 1:(length(Vm_set))){
  for (s in 2:length(Keq_set)){
    contour(w_kf_FD_sens_2Enz_rev[[(l-1)*length(kr1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=s,lty=l+0,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1.25)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1.25)
    print((l-1)*length(kr1_set)+s)
  }
}
text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1.5)
legend("bottomleft",title=c("Equilibrium constant"),legend=c(expression("Keq="~10^{-4}),expression("Keq="~10^{-2}),expression("Keq="~1),expression("Keq="~10^{2})),lty=1,col=c(2,3,1,4),ncol=2,cex=0.7)
title(main=substitute(paste("Low degradation rate - ", eta," =",eta_d,"/s, reversible first reaction"),list(eta_d=eta_d_set[1])),col.main="black",cex.main=1.1,font.main=1)



#coordonnées des fleches 1
x0<-c(0.2)
y0<-c(0.2)
x1<-c(0.65)
y1<-c(0.65)
arrows(x0,y0,x1,y1,code=2,lwd=2,lty=4)
mtext(text="Effect of increasing
      reversibility",at=c(0.7),line=-5.5,outer=FALSE,font=2)


###9.Analysis of Evolutionary Results

kcat_i<-10^-3
kf_i<-10^2
kr_i<-10^3
E_tot_conc=10^-3

KT_set=c(10^-5,10^-1)
Vm_set=c(10^-6,10^-3)

#Flux determination
fit<-function(kf,kr,kcat,E_conc){
  a<- (Vm*kf+kf*E_conc*kcat*(1+Se/KT))
  b<- (Vm*(kr+kcat-Se*kf)+kf*kcat*E_conc*(KT+Se))
  c<- -Vm*Se*(kr+kcat)
  delta<-b^2-4*a*c
  S_conc_eq=(-b+delta^(1/2))/(2*a)
  ES_conc_eq=E_conc*(kf*S_conc_eq)/(kr+kcat+kf*S_conc_eq)
  Flux=kcat*ES_conc_eq
  return(Flux)
}

Max_val<-list()
KT=KT_set[1]
Vm=Vm_set[1]
Se=10*KT
Max_val[["W"]]<-fit(10^10,10^3,10^6,10^-3)

##a.Importing simulation results
setwd(dir="~")
setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("W","H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(0,-0.1,-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,b))
      load(paste("/Users/florian/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,".rda",sep=""))
      if(b==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi00",sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][1])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][2+(long/3-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][3+2*(long/3-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]][3+2*(long/3-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"bias"=rep(mut_bias_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("W","H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
bias_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  bias_set_data[[f]]<-c()
  for (s in 1:12){
    n=param_set$Ne_set[s]
    b=param_set$bias[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    bias_set_data[[f]]<-c(bias_set_data[[f]],rep(b,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))
  }
}

data_fit_eq<-list()
col_headings <- c("Ne","bias","fitness")
for (f in flux_set){
  data_fit_eq[[f]]<-data.frame(cbind(Ne_set_data[[f]],bias_set_data[[f]],fit_evo_eq_summary[[f]]))
  names(data_fit_eq[[f]]) <- col_headings
}

##b.Checking the reach of steady-state
f="W"
n=Ne_set[1]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)

n=Ne_set[2]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)

n=Ne_set[3]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)

n=Ne_set[4]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)


f="H"
n=Ne_set[1]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)

n=Ne_set[2]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set) #Value 35 senseless, numerical bug
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)

n=Ne_set[3]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)

n=Ne_set[4]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)
b=mut_bias_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set)

##c.Plotting fitness values at steady-state

boxplot(log10(1-data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$fitness/Max_val[["W"]])~data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$Ne,col="green")
boxplot(log10(1-data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==-0.1,]$fitness/Max_val[["W"]])~data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==-0.1,]$Ne,add=TRUE,col="blue")
boxplot(log10(1-data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==-0.2,]$fitness/Max_val[["W"]])~data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==-0.2,]$Ne,add=TRUE,col="red")

##d. Plotting kinetic parameters at steady-state #Lack of the isoclines

#Weak flux,no bias
f="W"
b=mut_bias_set[1]
n=Ne_set[1]
plot(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=15,col="blue",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[2]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=16,col="blue",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[3]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=17,col="blue",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[4]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=18,col="blue",xlim=c(0,10),ylim=c(-4,6))

#Weak flux, low bias
f="W"
b=mut_bias_set[2]
n=Ne_set[1]
plot(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=19,col="red",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[2]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=20,col="red")
n=Ne_set[3]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=21,col="red")
n=Ne_set[4]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=22,col="red")

#Weak flux, high bias
f="W"
b=mut_bias_set[3]
n=Ne_set[1]
plot(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=19,col="red",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[2]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=20,col="red")
n=Ne_set[3]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=21,col="red")
n=Ne_set[4]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=22,col="red")

#High flux,no bias
f="H"
b=mut_bias_set[1]
n=Ne_set[1]
plot(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=15,col="blue",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[2]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=16,col="blue",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[3]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=17,col="blue",xlim=c(0,10),ylim=c(-4,6))
n=Ne_set[4]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=18,col="blue",xlim=c(0,10),ylim=c(-4,6))

#High flux, low bias
f="H"
b=mut_bias_set[2]
n=Ne_set[1]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=19,col="red")
n=Ne_set[2]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=20,col="red")
n=Ne_set[3]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=21,col="red")
n=Ne_set[4]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=22,col="red")

#High flux, high bias
f="H"
b=mut_bias_set[3]
n=Ne_set[1]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=19,col="red")
n=Ne_set[2]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=20,col="red")
n=Ne_set[3]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=21,col="red")
n=Ne_set[4]
points(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),pch=22,col="red")


##10.Manifold parameters involved: the case for reversibility
##And with degradation in the pathway

###With high affinity, testing the effect of reversibility
N_reso=10
#Defining parameters
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1<-10^4
kcat1<-10^0
Keq_set_sens<-c(10^-4,10^-2,1,10^2,10^4)
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso, ncol=N_reso)
P_Phi_eq_sens_2Enz_rev<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_rev<-list() #Equilibrium produt concentration

MaxconcP<-0

for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set_sens))){
    #kr=kr_set[1]
    Vm=Vm_set[l]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[l]
    Se=10*KT
    kr1<-sqrt(Keq_set_sens[s])*kcat1
    kinh1<-sqrt(Keq_set_sens[s])*kf1
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        kr2=kr2_set[1]
        kcat2<-kcat2_set[j]
        kf2=kf2_set[i]
        if(kf2/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          Vm1_pos<-kcat1*E1_conc
          KM1_pos<-(kcat1+kr1)/kf1
          K_I<-kinh1/kf1
          Vm1_neg<-kr1*E1_conc
          KM1_neg<-(kcat1+kr1)/kinh1
          Vm2<-kcat2*E2_conc
          KM2<-(kcat2+kr2)/kf2
          T1<-KT+Se
          T2<-1+Se/KT
          #init<-c(10^-15,10^-15)
          F1e<-expression(Vm1_pos*x/(KM1_pos+x+K_I*y)-Vm1_neg*y/(KM1_neg+y+x/K_I)-(Vm2/(KM2+y)+eta_d)*y)
          F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y)+eta_d)*y)
          Flist<-list(F1e,F2e)
          SSR<-RaphsonNewton(Flist,init,error)
          #print(paste(SSR$x,SSR$y))
          if(SSR$y>MaxconcP){
            MaxconcP<-SSR$y
          }
          tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
          tab_Conc_prod_eq[i,j]=SSR$y
        }
      }
    }
    P_Phi_eq_sens_2Enz_rev[[(s)+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_rev[[(s)+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_rev_app<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set_sens))){
    w_kf_FD_sens_2Enz_rev_app[[s+length(Keq_set_sens)*(l-1)]]<-round(P_Phi_eq_sens_2Enz_rev[[s+length(Keq_set_sens)*(l-1)]]/max(P_Phi_eq_sens_2Enz_rev[[s+length(Keq_set_sens)*(l-1)]]),digits=20)
  }
}

#Defining features of plots
plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1, 0.75, 1.5))


ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)
addtxt<-list(l=0.98,h=0.98,txt=c("A","B","C"),srt = 0,font=2,col="black")
addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)

contour(w_kf_FD_sens_2Enz_rev_app[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1))
title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 (kf)",ylab="log10 (kcat)",cex.lab=1)
for (s in 1:2){
  contour(w_kf_FD_sens_2Enz_rev_app[[s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=4-s,method = "edge",add=TRUE,xlim=c(0,1))
}
for (s in 4:5){
  contour(w_kf_FD_sens_2Enz_rev_app[[s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=s-2,method = "edge",add=TRUE,xlim=c(0,1))
}

text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
legend(xpd=TRUE,x=1.1,y=0.6,title=c("Equilibrium constant"),legend=c(expression("Keq="~10^{-4}),expression("Keq="~10^{-2}),expression("Keq="~1),expression("Keq="~10^{2}),expression("Keq="~10^{4})),lty=c(3,2,1,2,3),col=c(2,2,1,3,3),ncol=1,cex=0.7)

#coordonnées des fleches 1
x0<-c(0.3)
y0<-c(0.2)
x1<-c(0.8)
y1<-c(0.7)
arrows(x0,y0,x1,y1,code=2,lwd=2,lty=4)
mtext(text="Effect of increasing
      reversibility",at=c(0.85),line=-5.4,outer=FALSE,font=2,cex=0.75)
axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
     labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1.25)
axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
     labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1.25)

