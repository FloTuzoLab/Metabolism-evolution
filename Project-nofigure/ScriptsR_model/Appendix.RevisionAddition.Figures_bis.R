setwd(dir="~")
setwd(dir="enzyme-evolution/ScriptsR_model")
setwd(dir="Desktop/enzyme-evolution/ScriptsR_model")
source("AP11.Plotting_multiple_images.R")
ncol=128
ncontour=5
palet<-colorRampPalette(c("white"))(ncol)
jet.colors <- colorRampPalette(c("steelblue", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

Max_val<-list()
KT=KT_set[2]
Vm=Vm_set[2]
Se=10*KT
Max_val[["H"]]<-fit(10^10,10^3,10^6,10^-3)

##a.Importing simulation results
setwd(dir="~")
#setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("H")#,"H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(-0.1,-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,c))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,"EnzVar.rda",sep=""))
      if(c==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      tempo_Econc<-c()
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]])
        tempo_Econc<-c(tempo_Econc,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][1])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][2+(long/4-1)])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][3+2*(long/4-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][4+3*(long/4-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]][4+3*(long/4-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("Econc_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_Econc)
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
flux_set=c("H")#,"H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
bias_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  bias_set_data[[f]]<-c()
  for (s in 1:8){
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

par(mfrow=c(5,2),mai=c(0.55,0.75,0.45,0.45))
##b.Checking the reach of steady-state
options(digits=20)
f="H"
n=Ne_set[1]
#logn<-log(n,10)
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main="Moderate mutational bias (b = - 0.1)")
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main="High mutational bias (b = - 0.2)")
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)))/2
text(x=12,y=-3.6,labels=expression(paste("Ne=",10^2)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[2]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^3)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[3]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^4)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[4]
b=mut_bias_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
b=mut_bias_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^5)),srt=-90,xpd=NA,font=2,cex=1)

##Fit
#boxplot(log10(1-data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$fitness/max(data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$fitness))~data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$Ne,col=c("green"),ylim=c(-10,-1),xlab="Ne",ylab="Drift strength")
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$bias==-0.1,]$fitness/Max_val[["H"]])~data_fit_eq[["H"]][data_fit_eq[["H"]]$bias==-0.1,]$Ne,col=c("blue"),ylim=c(-7,-1),xlab="Ne",ylab="Drift strength")
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$bias==-0.2,]$fitness/Max_val[["H"]])~data_fit_eq[["H"]][data_fit_eq[["H"]]$bias==-0.2,]$Ne,col=c("red"),ylim=c(-7,-1),xlab="Ne",ylab="Drift strength")
text(x=5,y=-4,labels=paste("Steady-state flux"),srt=-90,xpd=NA,font=1,cex=1)

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "Evo_SteadyState_ConcCrow_HighF.jpeg", width = 700*3,height=1400*3,res=100*3,type="cairo")

N_reso=250#250
KT_set=c(10^-1) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^(2.8)
E_tot_conc=10^-3
##Defining variables and set of parameters
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        KT=KT_set[s]
        Se=10*KT
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),80)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),80)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),80)
          delta<-mpfr(b^2-4*a*c,80)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),80)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),80)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_kcatKM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kcatKM[[s+length(Vm)*(l-1)]]<-w_kcatKM_tab
  }
}


##Plot with low mutational bias
f="H"
b=mut_bias_set[1]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.75,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),
             cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=2.75,globcex=0.5)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=8)
text(0.325,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.425,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
text(0.525,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
text(0.625,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
axis(side=1,labels=FALSE,tck=FALSE)
axis(side=4,labels=FALSE,tck=FALSE)
legend(x=1.08,y=0.9,title="Mutational bias",legend=c("b = - 0.1","b = - 0.2"),pch=c(19,17),col=1,xpd=NA,cex=1.75)


#High mutational bias
b=mut_bias_set[2]
i=1
for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}

xyData=Evo_data
xlim=c(0,10)
ylim=c(-4,6)
for (i in 1:length(xyData)){
  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
         col=i,pch=17,cex=0.8)
}

text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
title(main="Mutation-selection-drift balance",cex.main=1.5)

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "Evo_Results_Crowding.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")

addtxt2<-list(l=0.6,h=-2.05,txt=c("A","B","C"),srt = 0,font=2,col="black",cex=2)

b=mut_bias_set[1]
i=1
Evo_Conc_data<-c()
Evo_Conc<-list()
for (n in Ne_set){
  Evo_Conc[[i]]<-as.data.frame(cbind(n,log10(get(paste("Econc_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  Evo_Conc_data<-rbind(Evo_Conc_data,Evo_Conc[[i]])
  i=i+1
}
par(mfrow=c(2,1),mai = c(1, 0.8, 0.75, 1))
boxplot(Evo_Conc_data$V2~log(Evo_Conc_data$n,10),xlab="log10 (Ne)",ylab="Enzyme concentration (log10 M)",pch=19,col=c("blue"),cex.lab=1,cex.axis=1)
title(main="Moderate mutational bias (b = - 0.1)",cex.main=0.75)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
b=mut_bias_set[2]
i=1
Evo_Conc_data<-c()
Evo_Conc<-list()
for (n in Ne_set){
  Evo_Conc[[i]]<-as.data.frame(cbind(n,log10(get(paste("Econc_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  Evo_Conc_data<-rbind(Evo_Conc_data,Evo_Conc[[i]])
  i=i+1
}
boxplot(Evo_Conc_data$V2~log(Evo_Conc_data$n,10),xlab="log10 (Ne)",ylab="Enzyme concentration (log10 M)",pch=19,col=c("red"),cex.lab=1,cex.axis=1)
title(main="High mutational bias (b = - 0.2)",cex.main=0.75)
text(addtxt2$l,addtxt2$h,addtxt2$txt[3],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "Evo_Conc_Results_Crowding.jpeg", width = 530*3,height=900*3,res=100*3,type="cairo")

##Simulations with cost of concentration
###Rev1.a ConcCost
kcat_i<-10^-3
kf_i<-10^2
kr_i<-10^3
#E_tot_conc_i=10^-3

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
Max_val[["H"]]<-fit(10^10,10^3,10^6,10^-3)
cost_set<-c(10^2,10^3,10^4,10^5)

##a.Importing simulation results
setwd(dir="~")
#setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("H")#,"H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(c in cost_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,c))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_co1_",c,"EnzVarCostONLY.rda",sep=""))
      if(c==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""),evo_results)
      }
      tempo_Econc<-c()
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]])
        tempo_Econc<-c(tempo_Econc,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][1])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][2+(long/4-1)])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][3+2*(long/4-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][4+3*(long/4-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[g]][[e]][4+3*(long/4-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("Econc_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_Econc)
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"cost"=rep(cost_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("H")#,"H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
cost_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  cost_set_data[[f]]<-c()
  for (s in 1:16){
    n=param_set$Ne_set[s]
    c=param_set$cost[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    cost_set_data[[f]]<-c(cost_set_data[[f]],rep(c,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))
  }
}

data_fit_eq<-list()
col_headings <- c("Ne","cost","fitness")
for (f in flux_set){
  data_fit_eq[[f]]<-data.frame(cbind(Ne_set_data[[f]],cost_set_data[[f]],fit_evo_eq_summary[[f]]))
  names(data_fit_eq[[f]]) <- col_headings
}

par(mfrow=c(5,4),mai=c(0.55,0.55,0.35,0.3))
##b.Checking the reach of steady-state
options(digits=20)
f="H"
n=Ne_set[1]
#logn<-log(n,10)
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-2})))
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-3})))
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-4})))
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-5})))
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^2)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[2]
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^3)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[3]
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^4)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[4]
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^5)),srt=-90,xpd=NA,font=2,cex=1)

##Fit
#boxplot(log10(1-data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$fitness/max(data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$fitness))~data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$Ne,col=c("green"),ylim=c(-10,-1),xlab="Ne",ylab="Drift strength")
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^2,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^2,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^2,]$Ne,col=c("red"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^3,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^3,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^3,]$Ne,col=c("light grey"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^4,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^4,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^4,]$Ne,col=c("blue"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^5,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^5,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^5,]$Ne,col=c("green"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
text(x=5,y=-3.5,labels=paste("Steady-state flux"),srt=-90,xpd=NA,font=1,cex=1)

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "Evo_SteadyState_ConcCost_HighF.jpeg", width = 1100*3,height=1400*3,res=100*3,type="cairo")

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B","C"),srt = 0,font=2,col="black",cex=1)

Evo_data<-list()

##Plot with low mutational bias
par(mfrow=c(4,2),mai=c(0.55,0.75,0.45,0.45))
f="H"
b=mut_bias_set[1]
j=1
Evo_data_full_high_b<-list()
for (c in cost_set){
  i=1
  for (n in Ne_set){
    Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))))
    i=i+1
  }
  Evo_data_full_high_b[[j]]<-Evo_data
  j=j+1
}


addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B","C","D","E","F","G","H"),srt = 0,font=2,col="black",cex=2)
legtext1=c(2,3,4,5)
ctext_set<-c("H","M","L","VL")
cost_text<-c("High","Moderate","Low","Very low")
log10c_set<-c(2,3,4,5)

for (cost in 1:length(ctext_set)){
  print(cost)
  multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
               abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
               xyData=Evo_data_full_high_b[[cost]],TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=1,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),
               cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=1.5,globcex=0.5)
  levels=c()
  for(cont in 1:4){
    levels=c(levels,1-10^(-cont-1))
  }
  contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
  abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=8)
  text(0.325,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
  text(0.425,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
  text(0.525,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
  text(0.625,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
  title(main=substitute(paste(Cost_adj," protein cost (c =",10^{-c},")"),list(Cost_adj=cost_text[cost],c=log10c_set[cost])),cex.main=1.75)
  text(addtxt2$l,addtxt2$h,addtxt2$txt[-1+2*cost],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
  axis(side=1,labels=FALSE,tck=FALSE)
  axis(side=4,labels=FALSE,tck=FALSE)
  
  text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
  setwd(dir="~")
  setwd(dir="enzyme-evolution/Paper/Figures")
  dev.print(device = jpeg, file = gsub(" ", "",paste("Evo_Results_Cost_",ctext_set[cost],".jpeg")), width = 575*3,height=460*3,res=100*3,type="cairo")
  dev.off()
}

col_set<-c("grey","red","green","blue")
par(mfrow=c(4,1),mai = c(0.55, 0.75, 0.75, 0.75))
for (cost in 1:length(ctext_set)){
  c=cost_set[cost]
  i=1
  for (n in Ne_set){
    Evo_Conc[[i]]<-as.data.frame(cbind(n,log10(get(paste("Econc_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))))
    Evo_Conc_data<-rbind(Evo_Conc_data,Evo_Conc[[i]])
    i=i+1
  }
  
  b=mut_bias_set[1]
  i=1
  Evo_Conc_data<-c()
  Evo_Conc<-list()
  for (n in Ne_set){
    Evo_Conc[[i]]<-as.data.frame(cbind(n,log10(get(paste("Econc_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))))
    Evo_Conc_data<-rbind(Evo_Conc_data,Evo_Conc[[i]])
    i=i+1
  }
  boxplot(Evo_Conc_data$V2~log(Evo_Conc_data$n,10),xlab="log10 (Ne)",ylab="Enzyme concentration (log10 M)",pch=19,col=col_set,cex.lab=1.25,cex.axis=1.25)
  title(main=substitute(paste(Cost_adj," protein cost (c =",10^{-c},")"),list(Cost_adj=cost_text[cost],c=log10c_set[cost])),cex.main=1.25)
  pos_text<-max(Evo_Conc_data$V2)
  #text(x=6,y=pos_text,labels=expression(paste("Ne=",10^5)),srt=-90,xpd=NA,font=2,cex=1)
  text(x=4.4,y=pos_text-0.05,addtxt2$txt[2*cost],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1.5)
}
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = gsub(" ", "",paste("Evo_Conc_Results_Cost.jpeg")), width = 500*3,height=1860*3,res=100*3,type="cairo")



###Cost&Crowding
##a.Importing simulation results
setwd(dir="~")
#setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("H")#,"H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(c in cost_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,c))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_co1_",c,"EnzVarCostCrow.rda",sep=""))
      if(c==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""),evo_results)
      }
      tempo_Econc<-c()
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]])
        tempo_Econc<-c(tempo_Econc,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][1])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][2+(long/4-1)])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][3+2*(long/4-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[10]][[e]][4+3*(long/4-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_c",c,sep=""))[[g]][[e]][4+3*(long/4-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("Econc_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_Econc)
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"cost"=rep(cost_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("H")#,"H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
cost_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  cost_set_data[[f]]<-c()
  for (s in 1:16){
    n=param_set$Ne_set[s]
    c=param_set$cost[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    cost_set_data[[f]]<-c(cost_set_data[[f]],rep(c,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))
  }
}

data_fit_eq<-list()
col_headings <- c("Ne","cost","fitness")
for (f in flux_set){
  data_fit_eq[[f]]<-data.frame(cbind(Ne_set_data[[f]],cost_set_data[[f]],fit_evo_eq_summary[[f]]))
  names(data_fit_eq[[f]]) <- col_headings
}

par(mfrow=c(5,4),mai=c(0.55,0.55,0.35,0.3))
##b.Checking the reach of steady-state
options(digits=20)
f="H"
n=Ne_set[1]
#logn<-log(n,10)
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-2})))
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-3})))
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-4})))
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
title(main=expression(paste("c = ",10^{-5})))
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^2)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[2]
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^3)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[3]
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^4)),srt=-90,xpd=NA,font=2,cex=1)

n=Ne_set[4]
c=cost_set[1]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[2]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[3]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
c=cost_set[4]
boxplot(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)~get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$g_set,xlab="Time",ylab="Log10 Flux(M)",cex.axis=0.9)
pos_text<-(min(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit))+max(log10(get(paste("Fit_evo_time",f,"_Ne",n,"_c",c,sep=""))$tempo_time_Fit)))/2
text(x=12,y=pos_text,labels=expression(paste("Ne=",10^5)),srt=-90,xpd=NA,font=2,cex=1)

##Fit
#boxplot(log10(1-data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$fitness/max(data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$fitness))~data_fit_eq[["W"]][data_fit_eq[["W"]]$bias==0,]$Ne,col=c("green"),ylim=c(-10,-1),xlab="Ne",ylab="Drift strength")
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^2,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^2,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^2,]$Ne,col=c("red"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^3,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^3,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^3,]$Ne,col=c("light grey"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^4,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^4,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^4,]$Ne,col=c("blue"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
boxplot(log10(1-data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^5,]$fitness/max(data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^5,]$fitness))~
          data_fit_eq[["H"]][data_fit_eq[["H"]]$cost==10^5,]$Ne,col=c("green"),ylim=c(-6,-1),xlab="Ne",ylab="Drift strength",cex.axis=0.9)
text(x=5,y=-3.5,labels=paste("Steady-state flux"),srt=-90,xpd=NA,font=1,cex=1)

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "Evo_SteadyState_ConcCostCrow_HighF.jpeg", width = 1100*3,height=1400*3,res=100*3,type="cairo")

f="H"
b=mut_bias_set[1]
j=1
Evo_data<-list()
Evo_data_full_high_b<-list()
for (c in cost_set){
  i=1
  for (n in Ne_set){
    Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_c",c,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))))
    i=i+1
  }
  Evo_data_full_high_b[[j]]<-Evo_data
  j=j+1
}

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B","C","D","E","F","G","H"),srt = 0,font=2,col="black",cex=2)
legtext1=c(2,3,4,5)
ctext_set<-c("H","M","L","VL")
cost_text<-c("High","Moderate","Low","Very low")
log10c_set<-c(2,3,4,5)
plot.new()
for (cost in 1:length(ctext_set)){
  print(cost)
  multiplePlot("","","","","","",ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
               abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
               xyData=Evo_data_full_high_b[[cost]],TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=1,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),
               cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=1.5,globcex=0.5)
  levels=c()
  for(cont in 1:4){
    levels=c(levels,1-10^(-cont-1))
  }
  contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
  abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=8)
  text(0.325,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
  text(0.425,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
  text(0.525,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
  text(0.625,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
  title(main=substitute(paste(Cost_adj," protein cost (c =",10^{-c},")"),list(Cost_adj=cost_text[cost],c=log10c_set[cost])),cex.main=1.75)
  text(addtxt2$l,addtxt2$h,addtxt2$txt[-1+2*cost],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
  axis(side=1,labels=FALSE,tck=FALSE)
  axis(side=4,labels=FALSE,tck=FALSE)
  
  text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
  setwd(dir="~")
  setwd(dir="enzyme-evolution/Paper/Figures")
  dev.print(device = jpeg, file = gsub(" ", "",paste("Evo_Results_CostCrow_",ctext_set[cost],".jpeg")), width = 575*3,height=460*3,res=100*3,type="cairo")
  dev.off()
}

N_reso=250#250
KT_set=c(10^-1) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^(2.8)
E_tot_conc=10^-5
##Defining variables and set of parameters
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        KT=KT_set[s]
        Se=10*KT
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),80)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),80)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),80)
          delta<-mpfr(b^2-4*a*c,80)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),80)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),80)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_kcatKM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kcatKM[[s+length(Vm)*(l-1)]]<-w_kcatKM_tab
  }
}


col_set<-c("grey","red","green","blue")
par(mfrow=c(4,1),mai = c(0.55, 0.75, 0.75, 0.75))
for (cost in 1:length(ctext_set)){
  c=cost_set[cost]
  i=1
  for (n in Ne_set){
    Evo_Conc[[i]]<-as.data.frame(cbind(n,log10(get(paste("Econc_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))))
    Evo_Conc_data<-rbind(Evo_Conc_data,Evo_Conc[[i]])
    i=i+1
  }
  
  b=mut_bias_set[1]
  i=1
  Evo_Conc_data<-c()
  Evo_Conc<-list()
  for (n in Ne_set){
    Evo_Conc[[i]]<-as.data.frame(cbind(n,log10(get(paste("Econc_ev_eq_fl",f,"_Ne",n,"_c",c,sep="")))))
    Evo_Conc_data<-rbind(Evo_Conc_data,Evo_Conc[[i]])
    i=i+1
  }
  boxplot(Evo_Conc_data$V2~log(Evo_Conc_data$n,10),xlab="log10 (Ne)",ylab="Enzyme concentration (log10 M)",pch=19,col=col_set,cex.lab=1.25,cex.axis=1.25)
  title(main=substitute(paste(Cost_adj," protein cost (c =",10^{-c},")"),list(Cost_adj=cost_text[cost],c=log10c_set[cost])),cex.main=1.25)
  pos_text<-max(Evo_Conc_data$V2)
  #text(x=6,y=pos_text,labels=expression(paste("Ne=",10^5)),srt=-90,xpd=NA,font=2,cex=1)
  text(x=4.4,y=pos_text-0.05,addtxt2$txt[2*cost],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1.5)
}
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = gsub(" ", "",paste("Evo_Conc_Results_CostCrow.jpeg")), width = 500*3,height=1860*3,res=100*3,type="cairo")

xlim=c(0,10)
ylim=c(-4,6)
###For the main paper
ctext_set<-c("M")
cost_text<-c("Moderate")
log10c_set<-c(3)
for (cost in 1:length(ctext_set)){
  print(cost)
  multiplePlot("","","","","","",ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
               abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
               xyData=Evo_data_full_high_b[[2]],TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.8,subcex=1,labcex=1.8,axcex=0.9,
               cextext=2,legtitle="log10(Ne)",legtext = legtext1,legpos=c(1.125,0.5),col_data = c(1,2,3,4),pch_dat=c(rep(5,4)),colorkey=FALSE,legcex=2,globcex=0.5)
  
  levels=c()
  for(cont in 1:4){
    levels=c(levels,1-10^(-cont-1))
  }
  #contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels[4],labcex=1,col=c(4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
  #contour(w_kcatKM_p[[1]],xaxt="n",yaxt="n",levels=levels[3],labcex=1,col=c(3),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
  #contour(w_kcatKM_p[[2]],xaxt="n",yaxt="n",levels=levels[2],labcex=1,col=c(2),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
  #contour(w_kcatKM_p[[3]],xaxt="n",yaxt="n",levels=levels[1],labcex=1,col=c(1),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
  abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=8)
  #text(0.325,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
  #text(0.425,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
  #text(0.525,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
  #text(0.625,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
  #title(main=substitute(paste(Cost_adj," protein cost (c =",10^{-c},")"),list(Cost_adj=cost_text[cost],c=log10c_set[cost])),cex.main=1.75)
  #text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
  axis(side=1,labels=FALSE,tck=FALSE)
  axis(side=4,labels=FALSE,tck=FALSE)
  legend(x=1.08,y=0.9,title="Production cost",legend=c("c = high","c = low"),pch=c(5,8),col=1,xpd=NA,cex=1.5)
  text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
  xyData=Evo_data_full_high_b[[3]]
  for (i in 1:length(xyData)){
    points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
           (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
           col=i,pch=8,cex=0.8)
  }
  setwd(dir="~")
  setwd(dir="enzyme-evolution/Paper/Figures")
  dev.print(device = jpeg, file = gsub(" ", "",paste("Evo_Results_CostCrow_PaperMod",ctext_set[cost],".jpeg")), width = 575*3,height=460*3,res=100*3,type="cairo")
  dev.off()
}

  

N_reso=50#250
KT_set=c(10^-1) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=c(10^(2.8),10^(2.8),10^(2.8))
E_tot_conc_set=c(10^-4,10^-3.5,10^-3)
##Defining variables and set of parameters
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        E_tot_conc<-E_tot_conc_set[l]
        KT=KT_set[s]
        Se=10*KT
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),80)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),80)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),80)
          delta<-mpfr(b^2-4*a*c,80)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),80)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),80)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_kcatKM_tab_c<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM_p<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_kcatKM_tab_c[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_kcatKM_tab_c[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kcatKM_p[[s+length(Vm)*(l-1)]]<-w_kcatKM_tab_c
  }
}
