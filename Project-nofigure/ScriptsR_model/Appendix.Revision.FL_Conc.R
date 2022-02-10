library(stringr)
#install.packages("deming")
library(deming)
library(plot3D)

###Influence de la concentration sur la fitness en fonction des propriétés cinétiques de l'enzyme
fit<-function(kf,kr,kcat,E_conc_m,cost_Enz){
  ##If deterministic relationship
    E_conc=E_conc_m
    Etot_back<-2.5*10^-3#broadly corresponds to E.coli or S.cerevisae
    #Net_Flux=0
    kf_act<-kf*10^(-(E_conc+Etot_back)/(5*10^-3))
    a<- (Vm*kf_act+kf_act*E_conc*kcat*(1+Se/KT))
    b<- (Vm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_conc*(KT+Se))
    c<- -Vm*Se*(kr+kcat)
    delta<-b^2-4*a*c
    S_conc_eq=(-b+delta^(1/2))/(2*a)
    ES_conc_eq=E_conc*(kf_act*S_conc_eq)/(kr+kcat+kf_act*S_conc_eq)
    Flux=kcat*ES_conc_eq
    #print(Flux)
    #Net_Flux=Flux*K/(E_conc+2*K)
    if(Flux-E_conc*cost_Enz>0){
      Net_Flux=Flux-E_conc*cost_Enz
    }
    else{
      Net_Flux=10^-20
    }
  ##If noisiness introduced
    return(Net_Flux)
  }
  
  
kr<-10^-3
KT=10^-1
Vm=10^-3
N_resoSM=1000
cost_set<-c(10^-5,10^-4,10^-3,10^-2)
log10E_tot<-seq(-10,0,length=N_resoSM)
E_tot_set<-10^log10E_tot
Se=10*KT
kf_set<-c(10^4,10^7,10^10)
kcat_set<-c(10^0,10^3,10^6)

w<-c()


t<-c()
len1<-3
len2<-3
for(i in 1:len1){
  if(i<2){
    t<-c(t,(i-1)*(len2)+seq(1,len2),len1*len2+1)
  }
  else{
    t<-c(t,(i-1)*len2+seq(1,len2),len1*len2+1)
  }
}
layout(matrix(t, len1, len2+1, byrow = TRUE), 
       widths=c(rep(1,len2),0.5), heights=c(rep(1,len1)))
col_set<-c("green","blue","grey","red")
addtxt<-list(x=-10,y=1,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")

par(xpd=FALSE)
for (f in 1:length(kf_set)){
  for (cat in 1:length(kcat_set)){
    for(c in 1:length(cost_set)){
      cost<-cost_set[c]
      kf<-kf_set[f]
      kcat<-kcat_set[cat]
      for (i in 1:N_resoSM){
        E_tot<-E_tot_set[i]
        w[i]<-fit(kf,kr,kcat,E_tot,cost)
      }
      if(c==1){
        plot(w/max(w)~log10E_tot,col=col_set[c],type="l",lty=1,xlab="log10 [E]",ylab="Relative fitness")
        abline(v=-0.45+log(f,10)*0.8,col=" light grey",lwd=40)
        if(f==1){
          title(main=substitute(paste("kcat=",10^{log10kcat},"/s"),list(log10kcat=log(kcat,10))))
        }
        if(cat==3){
          text(x=2,y=0.5,labels=substitute(paste("kf=",10^{log10kf},"M/s"),list(log10kf=log(kf,10))),srt=-90,xpd=NA,font=2,cex=1)
        }
        text(x=addtxt$x,y=addtxt$y,labels=addtxt$txt[(f-1)*length(kcat_set)+cat],font=2)
      }
      else{
        points(w/max(w)~log10E_tot,col=col_set[c],,type="l")
      }
      text(x=-0.5,y=0.5,labels="Crowding limit",font=2,srt=90)
    }
  }
}

par(xpd=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(0, 20), ylim=c(0, 20),xaxt='n', ann=FALSE,bty="n",yaxt='n')
legend("center",legend=c(expression(paste("c="~10^{-5})),
                           expression(paste("c="~10^{-4})),
                           expression(paste("c="~10^{-3})),
                           expression(paste("c="~10^{-2}))),cex=1.25,title="Protein cost",col=col_set,lty=1)

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "1DFit_Lanscape_Concentration.jpeg", width = 1000*3,height=900*3,res=100*3,type="cairo")

par(xpd=FALSE)
addtxt$y=9
for (f in 1:length(kf_set)){
  for (cat in 1:length(kcat_set)){
    for(c in 1:length(cost_set)){
      cost<-cost_set[c]
      kf<-kf_set[f]
      kcat<-kcat_set[cat]
      for (i in 1:N_resoSM){
        E_tot<-E_tot_set[i]
        w[i]<-fit(kf,kr,kcat,E_tot,cost)
      }
      if(c==1){
        plot(-log10(1-w/max(w+10^-12))~log10E_tot,col=col_set[c],type="l",lty=1,lwd=1,
             xlab="log10 [E]",ylab="-log 10 spread to max fitness",cex.lab=1.25)
        abline(v=-0.45+log(f,10)*0.8,col=" light grey",lwd=40)
        text(x=addtxt$x,y=9,labels=addtxt$txt[(f-1)*length(kcat_set)+cat],font=2)
        if(f==1){
          title(main=substitute(paste("kcat=",10^{log10kcat},"/s"),list(log10kcat=log(kcat,10))))
        }
        if(cat==3){
          text(x=2,y=5,labels=substitute(paste("kf=",10^{log10kf},"M/s"),list(log10kf=log(kf,10))),srt=-90,xpd=NA,font=2,cex=1)
        }
      }
      else{
        points(-log10(1-w/max(w+10^-12))~log10E_tot,col=col_set[c],,type="l")
      }
      text(x=-0.5,y=5,labels="Crowding limit",font=2,srt=90)
    }
  }
}
par(xpd=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(0, 20), ylim=c(0, 20),xaxt='n', ann=FALSE,bty="n",yaxt='n')
legend("center",legend=c(expression(paste("c="~10^{-5})),
                         expression(paste("c="~10^{-4})),
                         expression(paste("c="~10^{-3})),
                         expression(paste("c="~10^{-2}))),cex=1.25,title="Protein cost",col=col_set,lty=1)

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "1DFit_Lanscape_Scaling_Concentration.jpeg", width = 1000*3,height=900*3,res=100*3,type="cairo")

###5.Effect of background enzyme on the landscape of the second enzyme (assuming no reversibility)
N_reso=50
#Defining parameters
Vm_set=c(10^-6)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#,10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-5##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^2,10^4,10^10)
kcat1_set<-c(10^-2,10^0,10^6)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^-6,10^-4,10^-2)#degradation rate #Sensitivity study to be done relatively to Vm

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
    eta_d<-eta_d_set[l]
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

#Analysing the results
##Maximum of fitness
max(P_Phi_eq_sens_2Enz_background_eff[[1]])/max(P_Phi_eq_sens_2Enz_background_eff[[2]])

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff<-list()
for (l in 1:(length(eta_d_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}
eff_enz<-c("Low","Moderate","High")
multiplePlot("Degradation rate=","/s","First enzyme efficiency=","",signif(eta_d_set,4),eff_enz,ncol,log10kf_var,log10kcat_var,P_Phi_eq_sens_2Enz_background_eff,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=20,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey="COMMON",globcex=0.5,contcex=1,contourlab=TRUE,axcex=1,legcex=1)

multiplePlot("Degradation rate=","/s","First enzyme efficiency=","",signif(eta_d_set,4),eff_enz,ncol,log10kf_var,log10kcat_var, w_kf_FD_sens_2Enz_background_eff,
             abs="log10 (kf)",ord="log10 (kcat)",lev=c(0.9,0.99,0.999),palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey="COMMON",globcex=0.5,contcex=1,contourlab=TRUE,axcex=1)

multiplePlot("Degradation rate=","/s","First enzyme efficiency=","",signif(eta_d_set,4),eff_enz,ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz_background_eff,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey=TRUE,globcex=0.5,contcex=1,contourlab=TRUE,axcex=1)
sublegend<-c(expression(paste(eta,"=",10^-6,"/s")),
             expression(paste(eta,"=",10^-4,"/s")),
             expression(paste(eta,"=",10^-2,"/s")))


plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1.125, 0.75, 1.125))

contour(w_kf_FD_sens_2Enz_background_eff[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
title(main="",outer=FALSE,line=2.5,cex.main=1.25,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
for (l in 2:(length(eta_d_set))){
  for(s in 1:length(kcat1_set)){
    contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=s,lty=l+1,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
    print((l-1)*length(kcat1_set)+s)
  }
}


addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","  Perfect"),bty="n",lty=c(rep(1,3)),col=c(1,2),ncol=1,cex=0.65)
L2=legend(0.71,0.84,title="Degradation rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)
rect(L2$rect$left-0.02,L1$rect$top-L2$rect$h-L1$rect$h+0.02, L2$rect$left+L1$rect$w+0.03, L1$rect$top+0.02,col = "white",border="black", lty = 1, lwd = 1)
L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","  Perfect"),bty="n",lty=c(rep(1,2)),col=c(1,2),ncol=1,cex=0.65)
L2=legend(0.71,0.84,title="  Degradation rate",legend=sublegend,bty="n",lty=c(1,3,4),col=c(1),ncol=1,cex=0.65)

dev.print(device = jpeg, file = "2DFit_Landscape_2Enz_First_Enz_Influence_.jpeg", width = 575*3,height=475*3,res=100*3,type="cairo")

###5ter.Plotting fitness with a sigmoid toxicity
##High flux
N_reso=50
#Defining parameters
Vm_set=c(10^-3)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-3)#,10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^2,10^5,10^10,10^7)
kcat1_set<-c(10^-2,10^1,10^6,10^3)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^-6)#degradation rate #Sensitivity study to be done relatively to Vm
T_half_set<-c(10^-4,10^-2)#Toxicity penality

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-10^-6
E2_conc<-E_tot_conc

error=10^-15


tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_background_eff<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_background_eff<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[1]
    eta_d<-eta_d_set[1]
    T_half<-T_half_set[l]
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
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)#+eta_d
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)*(1-SSR$y/(SSR$y+T_half))
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

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff_tox<-list()
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff_tox[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C","D","E","F"),srt = 0,font=2,col="black",cex=1)
eff_enz<-c("Low","Moderate","High")
multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,P_Phi_eq_sens_2Enz_background_eff,TEXT_to_Add = addtxt2,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9,legcex=1,cextext=1.5)


multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz_background_eff_tox,TEXT_to_Add = addtxt2,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9,legcex=1,cextext=1.5)

plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1.125, 0.75, 1.125))

contour(w_kf_FD_sens_2Enz_background_eff_tox[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[4]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=4,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
title(main="",outer=FALSE,line=2.5,cex.main=1.25,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
for (l in 2:(length(eta_d_set))){
  for(s in 1:length(kcat1_set)){
    contour(w_kf_FD_sens_2Enz_background_eff_tox[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=s,lty=l+1,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
    print((l-1)*length(kcat1_set)+s)
  }
}

sublegend<-c(expression(10^{-4}~M),expression(10^{-2}~M))
addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
L1=legend(0.05,0.4,title="First enzyme",legend=c("Inefficient","Moderate","Good","Perfect"),bty="n",lty=c(rep(1,4)),col=c(1,2,4,3),ncol=1,cex=0.65)
L2=legend(0.05,0.175,title="    Toxicity rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)
rect(L2$rect$left-0.02,L1$rect$top-L2$rect$h-L1$rect$h+0.02, L2$rect$left+L1$rect$w+0.03, L1$rect$top+0.02,col = "white",border="black", lty = 1, lwd = 1)
L1=legend(0.05,0.4,title="First enzyme",legend=c("Inefficient","Moderate","Good","Perfect"),bty="n",lty=c(rep(1,4)),col=c(1,2,4,3),ncol=1,cex=0.65)
L2=legend(0.06,0.175,title="     Toxicity rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)

dev.print(device = jpeg, file = "2DFit_Landscape_2Enz_First_Enz_Influence&Tox&DegLow&Conc_Low.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")

##With a higher degradation rate

##High flux
N_reso=50
#Defining parameters
Vm_set=c(10^-3)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-3)#,10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^2,10^5,10^10,10^7)
kcat1_set<-c(10^-2,10^1,10^6,10^3)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^-3)#degradation rate #Sensitivity study to be done relatively to Vm
T_half_set<-c(10^-4,10^-2)#Toxicity penality

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-10^-6
E2_conc<-E_tot_conc

error=10^-15


tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_background_eff<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_background_eff<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[1]
    eta_d<-eta_d_set[1]
    T_half<-T_half_set[l]
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
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)#+eta_d
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)*(1-SSR$y/(SSR$y+T_half))
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

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff_tox<-list()
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff_tox[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C","D","E","F"),srt = 0,font=2,col="black",cex=1)
eff_enz<-c("Low","Moderate","High")
multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,P_Phi_eq_sens_2Enz_background_eff,TEXT_to_Add = addtxt2,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9,legcex=1,cextext=1.5)


multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz_background_eff_tox,TEXT_to_Add = addtxt2,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9,legcex=1,cextext=1.5)

plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1.125, 0.75, 1.125))

contour(w_kf_FD_sens_2Enz_background_eff_tox[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[4]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=4,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
title(main="",outer=FALSE,line=2.5,cex.main=1.25,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
for (l in 2:(length(eta_d_set))){
  for(s in 1:length(kcat1_set)){
    contour(w_kf_FD_sens_2Enz_background_eff_tox[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=s,lty=l+1,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
    print((l-1)*length(kcat1_set)+s)
  }
}

sublegend<-c(expression(10^{-4}~M),expression(10^{-2}~M))
addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
L1=legend(0.75,1,title="First enzyme",legend=c("Inefficient","Moderate","Good","Perfect"),bty="n",lty=c(rep(1,4)),col=c(1,2,4,3),ncol=1,cex=0.65)
L2=legend(0.77,0.75,title="    Toxicity rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)
rect(L2$rect$left-0.02,L1$rect$top-L2$rect$h-L1$rect$h+0.02, L2$rect$left+L1$rect$w+0.03, L1$rect$top+0.02,col = "white",border="black", lty = 1, lwd = 1)
L1=legend(0.75,1,title="First enzyme",legend=c("Inefficient","Moderate","Good","Perfect"),bty="n",lty=c(rep(1,4)),col=c(1,2,4,3),ncol=1,cex=0.65)
L2=legend(0.77,0.75,title="     Toxicity rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)

dev.print(device = jpeg, file = "2DFit_Landscape_2Enz_First_Enz_Influence&Tox&DegMod&Conc.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")




plot.new()
log10Phi<-seq(-15,-6,length=100)
Phi<-10^log10Phi
kcat=10^-0
E_tot<-10^-3

kf=10^5
kr=10^3

eta_d<-10^-2

P1_eq<-c()
P1_test<-c()
f<-c()
f_test<-c()
for (i in 1:100){
  Vm=mpfr(kcat*E_tot,120)
  KM<-mpfr((kcat+kr)/kf,120)
  a=mpfr(eta_d,120)
  b=mpfr((Vm+eta_d*KM-Phi[i]),120)
  c=mpfr(-Phi[i]*KM,120)
  Delta=mpfr(b^2-4*a*c,120)
  P1_eq[i]<-as.numeric(mpfr((-b+sqrt(Delta))/(2*a),120))
  P1_test[i]=as.numeric(-c/b)
  f[i]=as.numeric(mpfr(Vm*(P1_eq[i])/(P1_eq[i]+KM),120))
  f_test[i]=as.numeric(1-mpfr(eta_d*KM/(Vm+eta_d*KM),120))#-Phi[i]
}
f/Phi
P1_eq
P1_test
(P1_eq-P1_test)/P1_test
plot(log(eta_d*P1_eq/Phi,10)~log10Phi)
plot(log(f/Phi,10)~log10Phi)
KM
eta_d*P1_eq[100]





##4.Best value for the first enzyme followed by a variable second enzyme with sensitivy study : drawing isoclines of concentration
##With degradation in the pathway
N_reso=100
#Defining parameters
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5,5*10^-3)#
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
eta_d_set<-c(10^-6,10^-4,10^-2)#degradation rate #Sensitivity study to be done relatively to Vm

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
    Vm=Vm_set[l]
    eta_d<-eta_d_set[s]
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
    P_Phi_eq_sens_2Enz[[s+length(eta_d_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz[[s+length(eta_d_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}


#b.Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(eta_d_set))){
    w_kf_FD_sens_2Enz[[s+length(eta_d_set)*(l-1)]]<-round(P_Phi_eq_sens_2Enz[[s+length(eta_d_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz[[s+length(eta_d_set)*(l-1)]]),digits=20)
  }
}

#image.plot(P_Phi_eq[[1]],legend.shrink=1)
options(digits=10)
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B","C"),srt = 0,font=2,col="black",cex=1)

options(digits=3)
deg=round(eta_d_set*Vm_set[1],5)
round(deg,5)
options(digits=10)
subtitle<-list(substitute(paste("Low degradation rate - ",, eta," =",eta_d,"/s, irreversible first reaction"),list(eta_d=eta_d_set[1])),
               substitute(paste("Moderate degradation rate - ", eta," =",eta_d,"/s, irreversible first reaction"),list(eta_d=eta_d_set[2])),
               substitute(paste("High degradation rate - ", eta," =",eta_d,"/s, irreversible first reaction"),list(eta_d=eta_d_set[3])))

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
FitLandscape_etahigh_noreverse<- multiplePlot("KT=","","Vm=","M/s",KT_set[1],signif(Vm_set[1],4),ncol,log10kf_var,log10kcat_var,list(w_kf_FD_sens_2Enz[[3]]),
                                              abs="log10 kf",ord="log10 kcat",lev=c(0.9),palette=pal,image=TRUE,scale="AUTO",labcex=1,globcex=0.5,axcex=1,,contcex=1,
                                              TEXT_to_Add=addtxt2,meth="flattest",colorkey="TRUE",contourlab=FALSE)+
  contour(w_kf_FD_generic[[1]],lev=c(0.9),add=TRUE,lty=4,labcex=1,method='edge',col="white",labels=c("1st enz")) 
dev.print(device = jpeg, file = "2DFitLandscape_etahigh_noreverse.jpeg", width = 575*3,height=460*3,res=300,type="cairo")

Conc_prod_eq_sens_2EnzWF<-list(Conc_prod_eq_sens_2Enz[[1]],Conc_prod_eq_sens_2Enz[[2]],Conc_prod_eq_sens_2Enz[[3]])
Conc_prod_eq_sens_2EnzHF<-list(Conc_prod_eq_sens_2Enz[[4]],Conc_prod_eq_sens_2Enz[[5]],Conc_prod_eq_sens_2Enz[[6]])

multiplePlot("","","","",c(""),c("","",""),ncol,log10kf_var,log10kcat_var, Conc_prod_eq_sens_2EnzHF,
             abs="log10 kf",ord="log10 kcat",palette=pal,image=TRUE,scale="AUTO",labcex=1.5,globcex=0.5,axcex=1,contcex=1,
             TEXT_to_Add=addtxt2,meth="flattest",colorkey="TRUE",contourlab=FALSE,cextext=1.5)

#contour(Conc_prod_eq_sens_2Enz[[1]],xaxt="n",yaxt="n",levels=c(10^-4),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)

dev.print(device = jpeg, file = "2DConcentrationLandscape_degratesHF.jpeg", width = 1400*3,height=390*3,res=300,type="cairo")


