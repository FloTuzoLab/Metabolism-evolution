setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R") #To plot multiple gradient and color plots
source("AP2.RaphsonNewton.R") #To find equilibrium for 2 variables problems
require("Rmpfr") #to find equilibrium when very different values (dozens orders of magnitude) matter

#Defining features of plots
ncol=128
ncontour=5
jet.colors <- colorRampPalette(c("blue", "cyan", "yellow", "red")) #palette for fitness
fit.pal<-jet.colors(ncol)
conc.pal<-colorRampPalette(c("green","red"))(ncol)
###2 pathways enzymes: finding the equilibrium using Raphson-Newton method

##1.With both enzymes set to the same value, maximizing the infuence of different parameters

#a.Equilibrium calculus
N_reso=100
##Defining parameters
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=10^-2#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-10^3#c(10^-1,10^5) #An intermediate value seems reasonable to show general results
E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Inhibition through non-competitive inhibitors
ko<-10^-2
ki<-1
eta_d<-5*10^-2#degradation rate

##Defining the whole variable set of parameters for both inside reactions
kf1_set<-kf
kcat1_set<-kcat
kr1_set<-kr_set
kf2_set<-kf
kcat2_set<-kcat
kr2_set<-kr_set
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

##Error ceiling for the equilibrium finding
error=10^-10

##Finding the steady-state with numerical Raphson Newton method
tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_id<-list() #Equilibrium product flux
Conc_prod_eq<-list() #Equilibrium produt concentration

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    for (i in 1:N_reso){
      #k_1=lr[l]*k2[j]
      kr1=kr1_set[l]
      kr2=kr2_set[l]
      Vm=Vm_set[s]
      for (j in 1:N_reso){
        print(paste(i,j))
        init<-c(10^-14,10^-14)
        #Defining intermediate parameters
        kf1<-kf1_set[i]
        kf2<-kf2_set[i]
        kcat1<-kcat1_set[j]
        kcat2<-kcat2_set[j]
        A1<-kcat1*kf1*E1_conc
        A2<-kcat1+kr1
        Ki<-ko/ki #inhibition constant
        A3<-kcat2*kf2*E2_conc
        A4<-kcat2+kr2
        B1<-KT+Se
        B2<-1+Se/KT
        #Steady-state equations F(x,y)=0 to solve
        F1e<-expression(A1*x/((1+Ki*y)*(A2+kf1*x))-A3*y/(A4+kf2*y))
        F2e<-expression(Vm*(Se-x)/(B1+B2*x)-A3*y/(A4+kf2*y))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        while (SSR$y>10^-2){
          #print("YES")
          Ki<-Ki*10 #increasing inhibition if product concentration overcomes 10mM (around )
          SSR<-RaphsonNewton(Flist,init,error)
        }#No inhibition needed
        #print(SSR$y)
        tab_P_Phi_eq[i,j]=A3*SSR$y/(A4+kf2*SSR$y)
        tab_Conc_prod_eq[i,j]=SSR$y
      }
    }
    P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    #Conc_prod_eq[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
    #log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10) abandoned variable that masks flux differences
  }
}

Res_enz_id<-list(tab_P_Phi_eq,tab_Conc_prod_eq)
#b.Analysing the results
##Ploting the equilibrium flux with regards to concentrations
Vm_rep<-list(Vm_set,Vm_set)
pal<-list(fit.pal,conc.pal)
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_rep,ncol,log10kf_var,log10kcat_var,Res_enz_id,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",scale="AUTO",palette=pal,ncont=10,image=TRUE,sub=c("Flux (M/s)","Product concentration (M)")) ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    w_kf[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]),digits=10)
  }
}

pal=list(fit.pal)
##Plotting the relative fitness
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set,ncol,log10kf_var,log10kcat_var,w_kf,
             "Relative fitness landscape for the system
{facilitated diffusion+ 2 reactions}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,scale=c(0,1),lev=c(0.5,
                                                        0.99,
                                                        0.9999,#below any Ne values
                                                        0.999999#,#lower values for Ne
             )
             ,image=TRUE)#rounding to 8 digits, higher Ne values



##2.Best value for the first enzyme followed by a variable second enzyme with sensitivy study
##With degradation in the pathway
N_reso=50
#Defining parameters
Se=10^-3
KT=0.0001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3 Vm=c(1e-06,1e-03) à montrer
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
eta_d_set<-c(10^1,10^4)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_Vm6<-list() #Equilibrium product flux
P_Phi_eq_Vm3<-list()
Conc_prod_eq_Vm6<-list() #Equilibrium produt concentration
Conc_prod_eq_Vm3<-list()
MaxconcP<-0

#for (l in 1:(length(kr2_set))){
for (d in 1:length(eta_d_set)){
  for (s in 1:(length(Vm_set))){
    for (i in 1:(N_reso)){
      Vm=Vm_set[s]
      eta_d<-eta_d_set[d]*Vm
      for (j in 1:N_reso){
        print(paste(i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[l]
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
    if(Vm==Vm_set[1]){
      P_Phi_eq_Vm6[[d+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
      Conc_prod_eq_Vm6[[d+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
    }
    if(Vm==Vm_set[2]){
      P_Phi_eq_Vm3[[d+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
      Conc_prod_eq_Vm3[[d+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
    }
  }
}

MaxconcP
Res_enz_mut_Vm6<-list(P_Phi_eq_Vm6[[1]],Conc_prod_eq_Vm6[[1]],P_Phi_eq_Vm6[[2]],Conc_prod_eq_Vm6[[2]])
pal<-list(fit.pal,conc.pal,fit.pal,conc.pal)
Vm_rep<-list(Vm_set[1],Vm_set[1])
subtitle<-c("Flux (M/s)","Product concentration (M)","Flux (M/s)","Product concentration (M)")

#b.Analysing the results
##Ploting the equilibrium flux with regards to concentrations #Vm=1e-06
multiplePlot("Degradation rate=","x Vm","Vm=","M/s",eta_d_set,Vm_rep,ncol,log10kf_var,log10kcat_var,Res_enz_mut_Vm6,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,ncont=10,image=TRUE,scale="AUTO",sub=subtitle) ##The unit and the value of Vm need to be checked and explained

Res_enz_mut_Vm3<-list(P_Phi_eq_Vm3[[1]],Conc_prod_eq_Vm3[[1]],P_Phi_eq_Vm3[[2]],Conc_prod_eq_Vm3[[2]])
Vm_rep<-list(Vm_set[2],Vm_set[2])

##Ploting the equilibrium flux with regards to concentrations #Vm=1e-03

multiplePlot("Degradation rate=","x Vm","Vm=","M/s",eta_d_set,Vm_rep,ncol,log10kf_var,log10kcat_var,Res_enz_mut_Vm3,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,ncont=10,image=TRUE,scale="AUTO",sub=subtitle) ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux for a low degradation rate and low flux
w_kf_Vm6_eta1<-list()
#for (l in 1:(length(kr2_set))){
  #for (s in 1:(length(eta_d))){
    w_kf_Vm6_eta1[[1]]<-round(P_Phi_eq_Vm6[[1]]/max(P_Phi_eq_Vm6[[1]]),digits=10)
  #}
#}
pal<-list(fit.pal,conc.pal)
Res_enz_mut_Vm6_eta1<-list(w_kf_Vm6_eta1[[1]],Conc_prod_eq_Vm6[[1]])
subtitle<-c("Flux (M/s)","Product concentration (M)")
##Plotting the relative fitness
multiplePlot("kr=","","Vm=","M/s",kr2_set,Vm_set,ncol,log10kf_var,log10kcat_var,Res_enz_mut_Vm6_eta1,
             "Relative fitness landscape for the system
{facilitated diffusion+ 2 reactions}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,scale="AUTO",lev=c(0.5,
                                                        0.99,
                                                        0.9999,#below any Ne values
                                                        0.999999#,#lower values for Ne
             ),image=TRUE,sub=subtitle)#rounding to 8 digits, higher Ne values


##2.Best value for the first enzyme followed by a variable second enzyme with sensitivy study
##With degradation in the pathway
N_reso=50
#Defining parameters
Vm_set=c(10^-6,10^(-4.5),10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(10^-5,10^-3,10^-1)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1<-10^10
kcat1<-10^6
kr1<-10^6
kinh1<-10^10
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


##3.Best value for the first enzyme followed by a variable second enzyme with sensitivy study
##With degradation in the pathway and testing the influence of backward reactions

N_reso=100
#Defining parameters
Se=10^-2
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3 Vm=c(1e-06,1e-03) à montrer
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1<-10^9
kcat1<-10^5
kr1<-c(10^4,10^5,10^6)
kinh1<-c(10^8,10^9,10^10)
K_I<-kf1*kcat1/(kr1*kinh1)
eta_d<-10^2
#eta_d_set<-c(10^1,10^4)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_Vm6<-list() #Equilibrium product flux
P_Phi_eq_Vm3<-list()
Conc_prod_eq_Vm6<-list() #Equilibrium produt concentration
Conc_prod_eq_Vm3<-list()
MaxconcP<-0

for (d in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
  for (s in 1:(length(Vm_set))){
    for (i in 1:(N_reso)){
      Vm=Vm_set[s]
      eta_d<-eta_d_set[d]*Vm
      for (j in 1:N_reso){
        print(paste(i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[l]
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
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y))*y)
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
    if(Vm==Vm_set[1]){
      P_Phi_eq_Vm6[[d+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
      Conc_prod_eq_Vm6[[d+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
    }
    if(Vm==Vm_set[2]){
      P_Phi_eq_Vm3[[d+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
      Conc_prod_eq_Vm3[[d+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
    }
  }
}

MaxconcP
Res_enz_mut_Vm6<-list(P_Phi_eq_Vm6[[1]],Conc_prod_eq_Vm6[[1]],P_Phi_eq_Vm6[[2]],Conc_prod_eq_Vm6[[2]])
pal<-list(fit.pal,conc.pal,fit.pal,conc.pal)
Vm_rep<-list(Vm_set[1],Vm_set[1])
subtitle<-c("Flux (M/s)","Product concentration (M)","Flux (M/s)","Product concentration (M)")

#b.Analysing the results
##Ploting the equilibrium flux with regards to concentrations #Vm=1e-06
multiplePlot("Degradation rate=","x Vm","Vm=","M/s",eta_d_set,Vm_rep,ncol,log10kf_var,log10kcat_var,Res_enz_mut_Vm6,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,ncont=10,image=TRUE,scale="AUTO",sub=subtitle) ##The unit and the value of Vm need to be checked and explained

Res_enz_mut_Vm3<-list(P_Phi_eq_Vm3[[1]],Conc_prod_eq_Vm3[[1]],P_Phi_eq_Vm3[[2]],Conc_prod_eq_Vm3[[2]])
Vm_rep<-list(Vm_set[2],Vm_set[2])

##Ploting the equilibrium flux with regards to concentrations #Vm=1e-03

multiplePlot("Degradation rate=","x Vm","Vm=","M/s",eta_d_set,Vm_rep,ncol,log10kf_var,log10kcat_var,Res_enz_mut_Vm3,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,ncont=10,image=TRUE,scale="AUTO",sub=subtitle) ##The unit and the value of Vm need to be checked and explained



##Calculating the relative fitness with regards to the maximum achievable flux for a low degradation rate and low flux
w_kf_Vm6_eta1<-list()
#for (l in 1:(length(kr2_set))){
#for (s in 1:(length(eta_d))){
w_kf_Vm6_eta1[[1]]<-round(P_Phi_eq_Vm6[[1]]/max(P_Phi_eq_Vm6[[1]]),digits=10)
#}
#}
pal<-list(fit.pal,conc.pal)
Res_enz_mut_Vm6_eta1<-list(w_kf_Vm6_eta1[[1]],Conc_prod_eq_Vm6[[1]])
subtitle<-c("Flux (M/s)","Product concentration (M)")
##Plotting the relative fitness
multiplePlot("kr=","","Vm=","M/s",10^3,Vm_set,ncol,log10kf_var,log10kcat_var,Res_enz_mut_Vm6_eta1,
             "Relative fitness landscape for the system
{facilitated diffusion+ 2 reactions}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,scale="AUTO",lev=c(0.5,
                                                                    0.99,
                                                                    0.9999,#below any Ne values
                                                                    0.999999#,#lower values for Ne
             ),image=TRUE,sub=subtitle)#rounding to 8 digits, higher Ne values




##4.Best value for the first enzyme followed by a variable second enzyme, inhibition by products
##With degradation in the pathway
#a.Equilibrium calculus
N_reso=10
##Defining parameters
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=10^-2#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-10^3#c(10^-1,10^5) #An intermediate value seems reasonable to show general results
E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Inhibition through non-competitive inhibitors
ko<-10^3
ki<-1
eta_d<-5*10^-2#degradation rate

##Defining the whole variable set of parameters for both inside reactions
kf1_set<-kf
kcat1_set<-kcat
kr1_set<-kr_set
kf2_set<-kf
kcat2_set<-kcat
kr2_set<-kr_set
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

##Error ceiling for the equilibrium finding
error=10^-8

##Finding the steady-state with numerical Raphson Newton method
tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_id<-list() #Equilibrium product flux
Conc_prod_eq<-list() #Equilibrium produt concentration

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    for (i in 1:N_reso){
      #k_1=lr[l]*k2[j]
      kr1=kr_set
      kr2=0
      Vm=Vm_set[s]
      for (j in 1:N_reso){
        print(paste(i,j))
        init<-c(10^-5,10^-10)
        #Defining intermediate parameters
        kf1<-kf1_set[8]
        kf2<-kf2_set[i]
        kcat1<-kcat1_set[8]
        kcat2<-kcat2_set[j]
        A1<-kcat1*kf1*E1_conc
        A2<-kcat1+kr1
        Ki<-ko/ki #inhibition constant
        A3<-kcat2*kf2*E2_conc
        A4<-kcat2+kr2
        B1<-KT+Se
        B2<-1+Se/KT
        #Steady-state equations F(x,y)=0 to solve
        F1e<-expression(A1*x/((1+Ki*10^-1)*(A2+kf1*x))-A3*y/((1+Ki*10^-1)*(A4+kf2*y)))
        F2e<-expression(Vm*(Se-x)/((B1+B2*x)*(1+Ki*y))-A3*y/((1+Ki*10^-1)*(A4+kf2*y)))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #while (SSR$y>10^-2){
          #print("YES")
          #Ki<-Ki*10 #increasing inhibition if product concentration overcomes 10mM (around )
          #SSR<-RaphsonNewton(Flist,init,error)
        #}#No inhibition needed
        #print(SSR$y)
        tab_P_Phi_eq[i,j]=A3*SSR$y/(A4+kf2*SSR$y)
        tab_Conc_prod_eq[i,j]=SSR$y
      }
    }
    P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    #Conc_prod_eq[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
    #log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10) abandoned variable that masks flux differences
  }
}

Res_enz_id<-list(tab_P_Phi_eq,tab_Conc_prod_eq)
#b.Analysing the results
##Ploting the equilibrium flux with regards to concentrations
Vm_rep<-list(Vm_set,Vm_set)
pal<-list(fit.pal,conc.pal)
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_rep,ncol,log10kf_var,log10kcat_var,Res_enz_id,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",scale="AUTO",palette=pal,ncont=10,image=TRUE,sub=c("Flux (M/s)","Product concentration (M)")) ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    w_kf[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]),digits=10)
  }
}

pal=list(fit.pal)
##Plotting the relative fitness
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set,ncol,log10kf_var,log10kcat_var,w_kf,
             "Relative fitness landscape for the system
{facilitated diffusion+ 2 reactions}, kr=1e-03, Vm=1e-06",
             "log10 kf","log10 kcat",palette=pal,scale=c(0,1),lev=c(0.5,
                                                                    0.99,
                                                                    0.9999,#below any Ne values
                                                                    0.999999#,#lower values for Ne
             )
             ,image=TRUE)#rounding to 8 digits, higher Ne values

