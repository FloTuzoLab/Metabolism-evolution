##Function to plot multiple images
library(lattice)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(plot.matrix)
library(fields)

##Function to plot multiple images
multiplePlot<-function(nvar1,unitV1,nvar2,unitV2,v1,v2,ncol,xlim,ylim,xyZ,tit,
                       abs,ord,palette,scale,ncont,lev,xyData,rego,TEXT_to_Add,col_data,pch_data,image,ncat,pcex,labcex,legtext,legpos,sub,subcex){
  if(length(palette)==1){
    plot.new()
    pchlegend<-c()
    collegend<-c()
    len1<-length(v1)
    len2<-length(v2)
    par(mfrow=c(len1,len2))
    for (l in 1:len1){
      for (s in 1:len2){
        if(l==1 && s==1){
          if(missing(tit)){
            par(mar=c(6,6,2,1.5),mai = c(0.75, 1, 0.75, 0.75),cex=0.5)
          }
          else{
            par(mar=c(6,6,2,1.5),mai = c(0.75, 1, 1.5, 0.75),cex=0.5)
          }
        }
        else if(l==1){
          if(missing(tit)){
            par(mar=c(6,5,2,1.5),mai = c(0.75, 1, 0.75, 0.75),cex=0.5)
          }
          else{
            par(mar=c(6,5,2,1.5),mai = c(0.75, 1, 1.5, 0.75),cex=0.5)
          }
        }
        else if(s==1){
          if(missing(tit)){
            par(mar=c(6,6,2,1.5),mai = c(0.75, 1, 0.75, 0.75),cex=0.5)
          }
          else{
            par(mar=c(6,6,2,1.5),mai = c(0.75, 1, 1.5, 0.75),cex=0.5)
          }
        }
        else{
          if(missing(tit)){
            par(mar=c(6,5,2,1.5),mai = c(0.75, 1, 0.75, 0.75),cex=0.5)
          }
          else{
            par(mar=c(6,5,2,1.5),mai = c(0.75, 1, 1.5, 0.75),cex=0.5)
          }
        }
        if(image==TRUE){
          if(missing(scale)){
            Maximal<-max(xyZ[[s+(l-1)*len2]])
            Minimal<-min(xyZ[[s+(l-1)*len2]])
            if(Minimal<0){
              newz.na <- 0-(Maximal)/ncol
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<0)] <- newz.na
            #ncol=as.integer((Maximal-Minimal)/Minimal)
            #palette<-colorRampPalette(palette)(ncol-1)
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n")
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex)
                }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black")
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                }
              }
            }
          }
          else if(length(scale)==1){
            if(scale=="AUTO"){
              Maximal<-max(xyZ[[s+(l-1)*len2]])
              Minimal<-min(xyZ[[s+(l-1)*len2]])
              if(Minimal<0){
                newz.na <- 0-(Maximal)/ncol
                xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<0)] <- newz.na
              #ncol=as.integer((Maximal-Minimal)/Minimal)
              #palette<-colorRampPalette(palette)(ncol-1)
                if(missing(sub)){
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n")
                  }
                  else{
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex)
                  }
                }
                else{
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black")
                  }
                  else{
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                  }
                }
              }
              else{
                if(missing(sub)){
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n")
                  }
                  else{
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1)
                  }
                }
                else{
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black")
                  }
                  else{
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol+1,col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                  }
                }
              }
            }
          }
          else{
            if(length(which(xyZ[[s+(l-1)*len2]][]>scale[2]))>0){
              newz.na <- scale[2]+(scale[2]-scale[1])/ncol # new z for NA
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]>scale[2])] <- newz.na
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n")
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],
                             zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                             #breaks=seq(min(scale[1]),max(scale[2])),
                             nlevel=ncol,
                             col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex)
                }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black")
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],
                             zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                             #breaks=seq(min(scale[1]),max(scale[2])),
                             nlevel=ncol,
                             col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                }
              }
            ##mtext()
            }
            else if(length(which(xyZ[[s+(l-1)*len2]][]<scale[1]))>0){
              newz.na <- scale[1]-(scale[2]-scale[1])/ncol # new z for NA
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<scale[1])] <- newz.na
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n")
                }
                  else{
                    image.plot(xyZ[[s+(l-1)*len2]],
                               zlim=c(newz.na,scale[2]),
                               #breaks=seq(min(scale[1]),max(scale[2])),
                               nlevel=ncol,
                               col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=0.5,cex.main=subcex)
                  }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black")
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],
                             zlim=c(newz.na,scale[2]),
                             #breaks=seq(min(scale[1]),max(scale[2])),
                             nlevel=ncol,
                             col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                }
              }
            }
            else{
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n")
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex)
                }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),ncol=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black")
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                }
              }
            }
          }
          if(missing(ncont)){
          }
          else{
            contour(xyZ[[s+(l-1)*len2]],nlevels=ncont, add = TRUE,labcex=1,method = "edge")
          }
          if(missing(lev)){
          }
          else{
            contour(xyZ[[s+(l-1)*len2]],levels=lev, add = TRUE,labcex=1,method = "edge")
          }
        }
        else{
          par(mfrow=c(1,1))
          if(missing(ncont)){
          }
          else{
            if(s==1 & l==1){
              contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,xaxt="n",yaxt="n",labcex=1,col=1,lty=4,method = "flattest")
              collegend[1]=2
            }
            else{
              contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,add=TRUE,xaxt="n",yaxt="n",labcex=1,col=s+(l-1)*len2,lty=4,method = "flattest")
              collegend[s+(l-1)*len2]=1+s+(l-1)*len2
            }
          }
          if(missing(lev)){
          }
          else{
            if(s==1 & l==1){
              contour(xyZ[[s+(l-1)*len2]],levels=lev,xaxt="n",yaxt="n",labcex=1,col=2,lty=4,method = "flattest")
              collegend[1]=2
            }
            else{
              contour(xyZ[[s+(l-1)*len2]],levels=lev,add=TRUE,xaxt="n",yaxt="n",labcex=1,col=1+s+(l-1)*len2,lty=4,method = "flattest")
              collegend[s+(l-1)*len2]=1+s+(l-1)*len2
            }
          }
        }
        
        if(missing(TEXT_to_Add)){
        }
        else{
          text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=1)
        }
        
        if(missing(rego)){
        }
        else{
          abline(coef=c(((rego$coefficients[1]-rego$coefficients[2]
                          *(0-xlim[1]))-ylim[1])/(ylim[2]-ylim[1]),
                        (rego$coefficients[2])),
                lwd=1,col="white")
        }
    
        if(missing(xyData)){
        }
        else{
          if(missing(col_data)){
            if(missing(ncat)){
              for (i in 1:length(xyData)){
                if(missing(pcex)){
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                       (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                       col=i+1,pch=20)
                  pchlegend[1]=20
                  collegend[i]=1
                }
                else{
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                       (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                       col=i+1,pch=20,cex=pcex)
                  pchlegend[1]=20
                  collegend[i]=1
                }
              }
            }
            else{
              print(pch_data)
              #for (i in 1:(length(xyData)/ncat)){
                for (j in 1:length(xyData)){
                  #for(i in 1:length(xyData[[j]][[1]])){
                    if(missing(pcex)){
                      if(missing(pch_data)){
                        points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                               (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                               col=col_data[j],pch=j)
                        pchlegend[j]=j
                        collegend[j]=col_data[j]
                      }
                      else{
                        points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                               (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                               col=col_data[j],pch=pch_data[j])
                        pchlegend[j]=pch_data[j]
                        collegend[j]=col_data[j]
                      }
                     }
                    else{
                      if(missing(pch_data)){
                      points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                          (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                          col=j,pch=j,cex=pcex)
                      pchlegend[j]=j
                      collegend[j]=col_data[j]
                      }
                      else{
                        points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                               (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                               col=col_data[j],pch=pch_data[j],cex=pcex)
                        pchlegend[j]=pch_data[j]
                        collegend[j]=col_data[j]
                      }
                    }
                  #}
              }
            }
          }
          else{
            if(missing(ncat)){
              for (i in 1:length(xyData)){
                if(missing(pcex)){
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                     (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                     col=col_data[i],pch=20)
                  pchlegend[1]=20
                  collegend[i]=col_data[i]
                }
              
                else{
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                       (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                       col=col_data[i],pch=20,cex=pcex)
                  pchlegend[1]=20
                  collegend[i]=col_data[i]
                }
              }
            }
            else{
              #for (i in 1:(length(xyData)/ncat)){
                for (j in 1:length(xyData)){
                    if(missing(pcex)){
                     points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                         (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                         col=col_data[j],pch=j)
                      pchlegend[j]=j
                      collegend[j]=col_data[j]
                    }
                    else{
                      points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                          (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                          col=col_data[j],pch=j,cex=pcex)
                      pchlegend[j]=j
                      collegend[j]=col_data[j]
                  }
                #}
              }
            }
          }
        #print(pchlegend[1])
          if(missing(TEXT_to_Add)){
          
          }
          else{
            text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=1)
          }
        }
        if(missing(legtext)){
        }
        else{
          if(image==TRUE){
            legend(legpos,legend=legtext,col=collegend,pch=pchlegend, cex=pcex)
          }
          else{
            legend(legpos,legend=legtext,col=collegend,lty=6, cex=pcex)
          }
        }
        if(image==TRUE){
          axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]),cex.axis=1.25)
          axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]),cex.axis=1.25)
        }
        else{
          if(s==1 & l==1){
            axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]),cex.axis=1.25)
            axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]),cex.axis=1.25)
          }
        }
        if(missing(labcex)){
          if(image==TRUE){
            if(missing(tit)){
              title(xlab=abs,ylab=ord,line=2.5)
            }
            else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
              title(main=tit,outer=T,line=-3,cex.main=1.25,font=2)
            }
            else{
              title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord)
            }
          }
          else{
            if(s==1 & l==1){
              if(missing(tit)){
                title(xlab=abs,ylab=ord,line=2.5)
              }
              else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
                title(main=tit,outer=T,line=-3,cex.main=1.25,font=2,xlab=abs,ylab=ord)
              }
              else{
                title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord)
              }
            }
          }
        }
        else{
          if(image==TRUE){
            if(missing(tit)){
              title(xlab=abs,ylab=ord,line=2.5,cex.lab=labcex)
            }
            else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
              title(main=tit,outer=T,line=-3,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
            }
            else{
              title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
            }
          }
          else{
            if(s==1 & l==1){
              if(missing(tit)){
                title(xlab=abs,ylab=ord,line=2.5,cex.lab=labcex)
              }
              else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
                title(main=tit,outer=T,line=-3,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
              }
              else{
                title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
              }
            }
          }
        }
      
        if (len1==1 & len2==1){
        }
        else{
          if(l==len1){
            if(image==TRUE){
              mtext(paste(nvar2,v2[s],unitV2),side=1,line=4,cex=0.75,font=1)
            }
            else{
            }
          }
          if(s==1){
            if(image==TRUE){
              mtext(paste(nvar1,v1[l],unitV1),side=2,line=4.5,cex=0.75,font=1)
            }
            else{
            }
          }
        }
      }
    }
  }
  else{
    plot.new()
    pchlegend<-c()
    collegend<-c()
    len1<-length(v1)
    len2<-length(v2)
    par(mfrow=c(len1,len2))
    for (l in 1:len1){
      for (s in 1:len2){
        if(l==1 && s==1){
          if(missing(tit)){
            par(mar=c(6,6,2,1.5),mai = c(0.75, 1.25, 0.75, 0.5),cex=0.5)
          }
          else{
            par(mar=c(6,6,2,1.5),mai = c(0.75, 1.25, 1.5, 0.5),cex=0.5)
          }
        }
        else if(l==1){
          par(mar=c(6,5,5,2.5),mai = c(0.75, 1, 1.5, 0.75),cex=0.5)
        }
        else if(s==1){
          par(mar=c(6,6,2,1.5), mai = c(1.25, 1.25, 1, 0.5),cex=0.5)
        }
        else{
          par(mar=c(6,5,2,2.5),mai = c(1.25, 1, 1, 0.75),cex=0.5)
        }
        if(image==TRUE){
          if(missing(scale)){
            Maximal<-max(xyZ[[s+(l-1)*len2]])
            Minimal<-min(xyZ[[s+(l-1)*len2]])
            if(Minimal<0){
              newz.na <- 0-(Maximal)/ncol
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<0)] <- newz.na
              #ncol=as.integer((Maximal-Minimal)/Minimal)
              #palette<-colorRampPalette(palette)(ncol-1)
              if(missing(sub)){
                image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[s+(l-1)*len2]]),xaxt="n",yaxt="n",legend.cex=1.5,cex.main=subcex)
              }
              else{
                image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[s+(l-1)*len2]],legend.cex=1.5),
                           xaxt="n",yaxt="n")
                title(main=sub[s+(l-1)*len2],line=1,cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
              }
            }
          }
          else if(length(scale)==1){
            if(scale=="AUTO"){
              Maximal<-max(xyZ[[s+(l-1)*len2]])
              Minimal<-min(xyZ[[s+(l-1)*len2]])
              if(Minimal<0){
                newz.na <- 0-(Maximal)/ncol
                xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<0)] <- newz.na
                #ncol=as.integer((Maximal-Minimal)/Minimal)
                #palette<-colorRampPalette(palette)(ncol-1)
                if(missing(sub)){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[s+(l-1)*len2]]),xaxt="n",yaxt="n",legend.cex=1.5,cex.main=subcex)
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[s+(l-1)*len2]]),legend.cex=1.5,
                             xaxt="n",yaxt="n")
                  title(main=sub[s+(l-1)*len2],line=1,cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                }
              }
              else{
                if(missing(sub)){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol,col=c(palette[[s+(l-1)*len2]]),xaxt="n",yaxt="n",legend.cex=1.5,cex.main=subcex)
                }
                else{
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol,col=c(palette[[s+(l-1)*len2]]),legend.cex=1.5,
                             xaxt="n",yaxt="n")
                  title(main=sub[s+(l-1)*len2],line=1,cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
                }
              }
            }
          }
          else{
            if(length(which(xyZ[[s+(l-1)*len2]][]>scale[2]))>0){
              newz.na <- scale[2]+(scale[2]-scale[1])/ncol # new z for NA
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]>scale[2])] <- newz.na
              if(missing(sub)){
                image.plot(xyZ[[s+(l-1)*len2]],
                           zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                           #breaks=seq(min(scale[1]),max(scale[2])),
                           nlevel=ncol,
                           col=c(palette[[s+(l-1)*len2]],"light gray"),xaxt="n",yaxt="n",legend.cex=1.5,cex.main=subcex)
              }
              else{
                image.plot(xyZ[[s+(l-1)*len2]],
                           zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                           #breaks=seq(min(scale[1]),max(scale[2])),
                           nlevel=ncol,
                           col=c(palette[[s+(l-1)*len2]],"light gray"),
                           xaxt="n",yaxt="n",legend.cex=1.5)
                title(main=sub[s+(l-1)*len2],line=1,cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
              }
              ##mtext()
            }
            else if(length(which(xyZ[[s+(l-1)*len2]][]<scale[1]))>0){
              newz.na <- scale[1]-(scale[2]-scale[1])/ncol # new z for NA
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<scale[1])] <- newz.na
              if(missing(sub)){
                image.plot(xyZ[[s+(l-1)*len2]],
                           zlim=c(newz.na,scale[2]),
                           #breaks=seq(min(scale[1]),max(scale[2])),
                           nlevel=ncol,
                           col=c("light gray",palette[[s+(l-1)*len2]]),xaxt="n",yaxt="n",legend.cex=1.5,cex.main=subcex)
              }
              else{
                image.plot(xyZ[[s+(l-1)*len2]],
                           zlim=c(newz.na,scale[2]),
                           #breaks=seq(min(scale[1]),max(scale[2])),
                           nlevel=ncol,
                           col=c("light gray",palette[[s+(l-1)*len2]]),
                           xaxt="n",yaxt="n",legend.cex=1.5)
                title(main=sub[s+(l-1)*len2],line=1,cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
              }
            }
            else{
              if(missing(sub)){
                image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[s+(l-1)*len2]]),xaxt="n",yaxt="n",legend.cex=1.5,cex.main=subcex)
              }
              else{
                image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[s+(l-1)*len2]]),
                           xaxt="n",yaxt="n",legend.cex=1.5)
                title(main=sub[s+(l-1)*len2],line=1,cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1)
              }
            }
          }
          if(missing(ncont)){
          }
          else{
            contour(xyZ[[s+(l-1)*len2]],nlevels=ncont, add = TRUE,labcex=1)
          }
          if(missing(lev)){
          }
          else{
            contour(xyZ[[s+(l-1)*len2]],levels=lev, add = TRUE,labcex=1)
          }
        }
        else{
          par(mfrow=c(1,1))
          if(missing(ncont)){
          }
          else{
            if(s==1 & l==1){
              contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,xaxt="n",yaxt="n",labcex=1,col=2,lty=6,method = "edge")
              collegend[1]=2
            }
            else{
              contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,add=TRUE,xaxt="n",yaxt="n",labcex=1,col=1+s+(l-1)*len2,lty=6,method = "edge")
              collegend[s+(l-1)*len2]=1+s+(l-1)*len2
            }
          }
          if(missing(lev)){
          }
          else{
            if(s==1 & l==1){
              contour(xyZ[[s+(l-1)*len2]],levels=lev,xaxt="n",yaxt="n",labcex=1,col=2,lty=6,method = "edge")
              collegend[1]=2
            }
            else{
              contour(xyZ[[s+(l-1)*len2]],levels=lev,add=TRUE,xaxt="n",yaxt="n",labcex=1,col=1+s+(l-1)*len2,lty=6,method = "edge")
              collegend[s+(l-1)*len2]=1+s+(l-1)*len2
            }
          }
        }
        
        if(missing(TEXT_to_Add)){
        }
        else{
          text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=1)
        }
        
        if(missing(rego)){
        }
        else{
          abline(coef=c(((rego$coefficients[1]-rego$coefficients[2]
                          *(0-xlim[1]))-ylim[1])/(ylim[2]-ylim[1]),
                        (rego$coefficients[2])),
                 lwd=1,col="white")
        }
        
        if(missing(xyData)){
        }
        else{
          if(missing(col_data)){
            if(missing(ncat)){
              for (i in 1:length(xyData)){
                if(missing(pcex)){
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                         col=i+1,pch=20)
                  pchlegend[1]=20
                  collegend[i]=i+1
                }
                else{
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                         col=i+1,pch=20,cex=pcex)
                  pchlegend[1]=20
                  collegend[i]=i+1
                }
              }
            }
            else{
              for (i in 1:(length(xyData)/ncat)){
                for (j in 1:ncat){
                  if(missing(pcex)){
                    points((xyData[[j+(i-1)*(length(xyData)/ncat)]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                           (xyData[[j+(i-1)*(length(xyData)/ncat)]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                           col=i+1,pch=15+j)
                    pchlegend[j]=15+j
                    collegend[i]=i+1
                  }
                  else{
                    points((xyData[[j+(i-1)*(length(xyData)/ncat)]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                           (xyData[[j+(i-1)*(length(xyData)/ncat)]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                           col=i+1,pch=15+j,cex=pcex)
                    pchlegend[j]=15+j
                    collegend[i]=i+1
                  }
                }
              }
            }
          }
          else{
            if(missing(ncat)){
              for (i in 1:length(xyData)){
                if(missing(pcex)){
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                         col=col_data[i],pch=20)
                  pchlegend[1]=20
                  collegend[i]=col_data[i]
                }
                
                else{
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                         col=col_data[i],pch=20,cex=pcex)
                  pchlegend[1]=20
                  collegend[i]=col_data[i]
                }
              }
            }
            else{
              for (i in 1:(length(xyData)/ncat)){
                for (j in 1:ncat){
                  if(missing(pcex)){
                    points((xyData[[j+(i-1)*(length(xyData)/ncat)]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                           (xyData[[j+(i-1)*(length(xyData)/ncat)]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                           col=col_data[i],pch=15+j)
                    pchlegend[j]=15+j
                    collegend[i]=col_data[i]
                  }
                  else{
                    points((xyData[[j+(i-1)*(length(xyData)/ncat)]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                           (xyData[[j+(i-1)*(length(xyData)/ncat)]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                           col=col_data[i],pch=15+j,cex=pcex)
                    pchlegend[j]=15+j
                    collegend[i]=col_data[i]
                  }
                }
              }
            }
          }
          #print(pchlegend[1])
          if(missing(TEXT_to_Add)){
            
          }
          else{
            text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=1)
          }
        }
        if(missing(legtext)){
        }
        else{
          if(image==TRUE){
            legend(legpos,legend=legtext,col=collegend,pch=pchlegend, cex=pcex)
          }
          else{
            legend(legpos,legend=legtext,col=collegend,lty=6, cex=pcex)
          }
        }
        if(image==TRUE){
          axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]),cex.axis=1.25)
          axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]),cex.axis=1.25)
        }
        else{
          if(s==1 & l==1){
            axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]),cex.axis=1.25)
            axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]),cex.axis=1.25)
          }
        }
        if(missing(labcex)){
          if(image==TRUE){
            if(missing(tit)){
              title(xlab=abs,ylab=ord,line=2.5)
            }
            else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
              title(main=tit,outer=T,line=-3,cex.main=1.25,font=1,xlab=abs,ylab=ord)
            }
            else{
              title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=1,xlab=abs,ylab=ord)
            }
          }
          else{
            if(s==1 & l==1){
              if(missing(tit)){
                title(xlab=abs,ylab=ord,line=2.5)
              }
              else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
                title(main=tit,outer=T,line=-3,cex.main=1.25,font=1,xlab=abs,ylab=ord)
              }
              else{
                title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=1,xlab=abs,ylab=ord)
              }
            }
          }
        }
        else{
          if(image==TRUE){
            if(missing(tit)){
              title(xlab=abs,ylab=ord,line=2.5,cex.lab=labcex)
            }
            else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
              title(main=tit,outer=T,line=-3,cex.main=1.25,font=1,xlab=abs,ylab=ord,cex.lab=labcex)
            }
            else{
              title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=1,xlab=abs,ylab=ord,cex.lab=labcex)
            }
          }
          else{
            if(s==1 & l==1){
              if(missing(tit)){
                title(xlab=abs,ylab=ord,line=2.5,cex.lab=labcex)
              }
              else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
                title(main=tit,outer=T,line=-3,cex.main=1.25,font=1,xlab=abs,ylab=ord,cex.lab=labcex)
              }
              else{
                title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=1,xlab=abs,ylab=ord,cex.lab=labcex)
              }
            }
          }
        }
        
        if (len1==1 & len2==1){
        }
        else{
          if(l==len1){
            if(image==TRUE){
              mtext(paste(nvar2,v2[s],unitV2),side=1,line=4,cex=0.75,font=1)
            }
            else{
            }
          }
          if(s==1){
            if(image==TRUE){
              mtext(paste(nvar1,v1[l],unitV1),side=2,line=4.5,cex=0.75,font=1)
            }
            else{
            }
          }
        }
      }
    }
  }
}

#Defining colors
#ncol=64
#pal<-colorRampPalette(c("hot pink","red","green"))(ncol)

