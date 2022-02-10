#3.Best value for the first enzyme followed by a variable second enzyme
##Using inhibition by reverse reaction and by limiting the flux through the membrane
N_reso=51
#Defining parameters
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(2,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-3,5,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
Vm_set=c(10^-6,10^-3)
kf1<-10^10
kcat1<-10^5
kr1<-10^3
Ki<-10^-6
kinh1<-10^8
#Max_conc<-c(10^-1,10^-4)

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^-2,10^5)
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
P_Phi_eq<-list() #Equilibrium product flux

for (l in 1:(length(kr2_set))){
  for (s in 1:(length(Vm_set))){
    for (i in 1:(N_reso)){
      Vm=Vm_set[s]
      for (j in 1:N_reso){
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[l]
        #print(paste(j,"kcat2",kcat2))
        Vm1p<-kcat1*E1_conc
        KM1p<-(kcat1+kr1)/kf1
        Vm1n<-kr1*E1_conc
        KM1n<-(kcat1+kr1)/kinh1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression((Vm1p*x/(KM1p+x+kinh1/kf1*y)-Vm1n*y/(KM1n+y+kf1/kinh1*x))-Vm2*y/(KM2+y))
        F2e<-expression(Vm/(1+y/Ki)*(Se-x)/(T1+x*T2)-Vm2*y/(KM2+y))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        print(paste(SSR$x,SSR$y))
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
      }
    }
    P_Phi_eq[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
  }
}

#b.Analysing the results
##Ploting the equilibrium flux with regards to concentrations
multiplePlot("kr=","","Vm=","M/s",kr2_set,Vm_set,ncol,log10kf_var,log10kcat_var,P_Phi_eq,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",ncont=10) ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    w_kf[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq[[s+length(Vm_set)*(l-1)]]),digits=10)
  }
}

##Plotting the relative fitness
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set,ncol,log10kf_var,log10kcat_var,w_kf,
             "Relative fitness landscape for the system
             {facilitated diffusion+ 2 reactions}",
             "log10 kf","log10 kcat",scale=c(0,1),lev=c(0.5,
                                                        0.99,
                                                        0.9999,#below any Ne values
                                                        0.999999#,#lower values for Ne
             ))#rounding to 8 digits, higher Ne values


#4.Best value for the first enzyme followed by a variable second enzyme
##Using inhibition by reverse reaction
N_reso=51
#Defining parameters
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
Vm_set=c(10^-6,10^-3)
kf1<-10^10
kcat1<-10^5
kr1<-10^3
kinh1<-10^15
Max_conc<-c(10^-1,10^-4)

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^-2,10^5)
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-8,10^-8)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
P_Phi_eq<-list() #Equilibrium product flux
P_Phi_eq_g<-list()

for (c in 1:length(Max_conc)){
  for (l in 1:(length(kr2_set))){
    for (s in 1:(length(Vm_set))){
      for (i in 1:(N_reso)){
        Vm=Vm_set[s]
        for (j in 1:N_reso){
          kinh1<-10^15
          kf2<-kf2_set[i]
          print(paste(i,"kf2",kf2))
          kcat2<-kcat2_set[j]
          kr2=kr2_set[l]
          print(paste(j,"kcat2",kcat2))
          A11<-kcat1*kf1*E1_conc
          A12<-kr1*kinh1*E1_conc
          A2<-kcat1+kr1
          A3<-kcat2*kf2*E2_conc
          A4<-kcat2+kr2
          B1<-KT+Se
          B2<-1+Se/KT
          init<-c(10^-15,10^-15)
          F1e<-expression((A11*x-A12*y)/(A2+kf1*x+kinh1*y)-A3*y/(A4+kf2*y))
          F2e<-expression(Vm*(Se-x)/(B1+B2*x)-A3*y/(A4+kf2*y))
          Flist<-list(F1e,F2e)
          SSR<-RaphsonNewton(Flist,init,error)
          print(paste(SSR$x,SSR$y))
          while ((SSR$x<10^-2 & SSR$y<(Max_conc[c])) & kinh1>10^3){
            print(paste(kinh1,SSR$x,SSR$y))
            kinh1=kinh1/10
            A12<-kr1*kinh1*E1_conc
            SSR<-RaphsonNewton(Flist,init,error)
          }
          tab_P_Phi_eq[i,j]=A3*SSR$y/(A4+kf2*SSR$y)
        }
      }
      P_Phi_eq[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    }
  }
  P_Phi_eq_g[[c]]<-P_Phi_eq
}

##Ploting the equilibrium flux with regards to concentrations
multiplePlot("kr=","","Vm=","M/s",kr2_set,Vm_set,ncol,log10kf_var,log10kcat_var,P_Phi_eq_g[[1]],
             "Absolute fitness landscape for the system
{facilitated diffusion + 2 chemical reactions}",
             "log10 kf","log10 kcat",ncont=5) ##The unit and the value of Vm need to be checked and explained

multiplePlot("kr=","","Vm=","M/s",kr2_set,Vm_set,ncol,log10kf_var,log10kcat_var,P_Phi_eq_g[[2]],
             "Absolute fitness landscape for the system
{facilitated diffusion + 2 chemical reactions}",
             "log10 kf","log10 kcat",ncont=5) ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    w_kf[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq[[s+length(Vm_set)*(l-1)]]),digits=10)
  }
}

##Plotting the relative fitness
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set,ncol,log10kf_var,log10kcat_var,w_kf,
             "Relative fitness landscape for the system
{facilitated diffusion + 2 chemical reactions}",
             "log10 kf","log10 kcat",scale=c(0,1),lev=c(0.5,
                                                        0.99,
                                                        0.9999,#below any Ne values
                                                        0.999999#,#lower values for Ne
             ))#rounding to 8 digits, higher Ne values
