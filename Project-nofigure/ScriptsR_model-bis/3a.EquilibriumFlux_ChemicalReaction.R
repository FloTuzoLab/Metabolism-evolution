setwd("/Users/florian/Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R")
#Defining features of plots
ncol=128
ncontour=5
pal<-colorRampPalette(c("hot pink","red","green"))(ncol)

##Analysis of the influence of reaction chemical parameters on equilibrium flux
N_reso=501
#Creation of two sets of parameters
log10kf_var<-c(0,10,2)
log10kcat_var<-c(-4,6,2)
kf<-10^seq(0,10,length.out=N_reso)
kcat<-10^seq(-4,6,length.out=N_reso)
kr_set=c(10^-1,10^5)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()
log10_P_Phi_eq<-list()
log10_ratio_Phi_Schem<-list()

####Reaction without accounting for facilitated diffusion

####a)In relation to k1 and k2
S_chem<-c(10^-6,10^-3) #Constant concentration (very inaccurate definition of the process)
E_tot_conc=10^-3

for (r in 1:(length(kr_set))){
  for (s in 1:length(S_chem)){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        kr<-kr_set[r]
        ES_conc_chem=E_tot_conc*(kf[i]*S_chem[[s]])/(kr+kcat[j]+kf[i]*S_chem[[s]])
        tab_P_Phi_eq[i,j]=kcat[j]*ES_conc_chem
      }
    }
    P_Phi_eq[[s+(r-1)*length(S_chem)]]=tab_P_Phi_eq
    log10_P_Phi_eq[[s+(r-1)*length(S_chem)]]=log(P_Phi_eq[[s+(r-1)*length(S_chem)]],10)
    #log10_ratio_Phi_Schem[[s+(r-1)*length(S_chem)]]=log10_P_Phi_eq[[s+(r-1)*length(S_chem)]]-log(S_chem[[s]],10)
  }
}

##Conclusion: with no saturation effect at the enter of the cell, there should be no saturation for the effect of enzymes
              #the saturation effect comes from the impossibility of increasing too much kcat and kf

##Fitness landscapes: different point of views
#multiplePlot("kr=","(kf)","[Se]=","M",kr_set,S_chem,ncol,log10kf_var,log10kcat_var,log10_ratio_Phi_eq,
             #"Fitness landscape for the chemical reaction",
             #"log10 kf","log10 kcat",ncont=ncontour,palette=pal)

multiplePlot("kr=","(kf)","[Se]=","M",kr_set,S_chem,ncol,log10kf_var,log10kcat_var,log10_P_Phi_eq,
             "Fitness landscape for the chemical reaction",
             "log10 kf","log10 kcat",ncont=ncontour,palette=pal)


###b)In relation to KM and k2  ##not reanalyzed as it does not represent the true phenomenon
log10KM_var<-c(-9,1,2)
log10k2_var<-c(-4,6,2)
k2<-10^seq(log10k2_var[1],log10k2_var[2],length.out=N_reso)
KM<-10^seq(log10KM_var[1],log10KM_var[2],length.out=N_reso)
lr<-c(0,1/10)
tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()
log10_P_Phi_eq<-list()

S_chem<-c(10^-8,10^-6,10^-4)

for (l in 1:(length(lr))){
  for (s in 1:length(S_chem)){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        k1<-(k_1[lr]+k2[j])/KM[i] #No interest as k1 compensates
        #the difference to obtain the same KM which is the only variable
        ES_conc_chem=E_tot_conc*(S_chem[[s]])/(KM[i]+S_chem[[s]])
        tab_P_Phi_eq[i,j]=k2[j]*ES_conc_chem
      }
    }
    P_Phi_eq[[s+(l-1)*3]]=tab_P_Phi_eq
    log10_P_Phi_eq[[s+(l-1)*3]]=log(P_Phi_eq[[s+(l-1)*3]],10)
    log10_ratio_Phi_Schem[[s+(l-1)*3]]=log10_P_Phi_eq[[s+(l-1)*3]]-log(S_chem[[s]],10)
  }
}

##Fitness landscapes: different point of views
multiplePlot("k-1=","(k1)","[Se]=","M",lr,S_chem,ncol,log10KM_var,log10k2_var,log10_ratio_Phi_Schem,
             "Fitness landscape for the chemical reaction relatively to KM",
             "log10 KM","log10 k2",ncont=ncontour)

multiplePlot("k-1=","(k1)","[Se]=","M",lr,S_chem,ncol,log10k1_var,log10k2_var,log10_ratio_Phi_Schem,
             "Fitness landscape for the chemical reaction relatively to KM",
             "log10 k1","log10 k2",c(-9,5))



##APPENDIX: Drawing of the fitness landscapes in 4D old version

plot1<-list()
grobplot<-list()
plist<-list()

for (l in 1:length(lr)){
  for (s in 1:length(S_chem)){
    plot1[[s+(l-1)*3]]<-{levelplot(log10_ratio_Phi_Schem[[s+(l-1)*3]],at=seq(min(-9), max(5),length=20),col.regions=pal,
                                   scales=list(x=list(at=c(1,1+(N_reso-1)/5,1+2*(N_reso-1)/5,1+3*(N_reso-1)/5,1+4*(N_reso-1)/5,1+5*(N_reso-1)/5),labels=c(0,2,4,6,8,10)),
                                               y=list(at=c(1,1+(N_reso-1)/5,1+2*(N_reso-1)/5,1+3*(N_reso-1)/5,1+4*(N_reso-1)/5,1+5*(N_reso-1)/5),labels=c(-4,-2,0,2,4,6))),
                                   xlab ="log k1 (base 10)",
                                   ylab ="log k2 (base 10)",
                                   contour=T,
                                   line = list(lwd=1,lty=2),
                                   auto.key=list(space="top", columns=3,title=expression(m^3/m^3)))
    }
    if(l==1){
      if(s==1){
        grobplot[[s+(l-1)*3]]<-arrangeGrob(plot1[[s+(l-1)*3]],
                                           top=paste("Substrate concentration=",S_chem[s]),
                                           left=paste("k-1=",lr[l],"(k1)"),padding = unit(0.5, "line"))
      }
      else{
        grobplot[[s+(l-1)*3]]<-arrangeGrob(plot1[[s+(l-1)*3]],
                                           top=paste("Substrate concentration=",S_chem[s]),padding = unit(0.5, "line"))
      }
    }
    else{
      if(s==1){
        grobplot[[s+(l-1)*3]]<-arrangeGrob(plot1[[s+(l-1)*3]],
                                           left=paste("k-1=",lr[l],"(k1)"),padding = unit(0.5, "line"))
      }
      else{
        grobplot[[s+(l-1)*3]]<-arrangeGrob(plot1[[s+(l-1)*3]],padding = unit(0.5, "line"))
      }
    }
  }
}

plots1<-grid.arrange(grobs=grobplot,ncol=3,
                     top=textGrob("Fitness landscapes for chemical reactions",
                                  gp=gpar(fontface = "bold", cex = 1.5)))
plots1

plot1<-list()

for (l in 1:length(lr)){
  for (s in 1:length(S_chem)){
    plot1[[s+(l-1)*3]]<-{levelplot(cuts=10,log10_P_Phi_eq[[s+(l-1)*3]],at=seq(min(-15), max(2),length=100),main ="Fitness landscape of
  cells", col.regions=colorRampPalette(c("purple","red" ,"dark green","green","blue")),
                                   scales=list(x=list(at=c(1,1+(N_reso-1)/5,1+2*(N_reso-1)/5,1+3*(N_reso-1)/5,1+4*(N_reso-1)/5,1+5*(N_reso-1)/5),labels=c(-9,-7,-5,-3,-1,1)),
                                               y=list(at=c(1,1+(N_reso-1)/5,1+2*(N_reso-1)/5,1+3*(N_reso-1)/5,1+4*(N_reso-1)/5,1+5*(N_reso-1)/5),labels=c(-3,-1,1,3,5,7))),
                                   xlab ="log KM (base 10)",
                                   ylab ="log k2 (base 10)",
                                   auto.key=list(space="top", columns=3,title=expression(m^3/m^3)))
    }
  }
}
plots1<-grid.arrange(plot1[[1]],plot1[[2]],plot1[[3]],plot1[[4]],plot1[[5]],plot1[[6]],ncol=3)
plots1