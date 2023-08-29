library(fChange)
library(reshape2)
library(fda)
library(sde)
library(MASS)

data_temp=list()
data_temp[[1]]=Australian_Temp$Sydney
data_temp[[2]]=Australian_Temp$Melbourne
data_temp[[3]]=Australian_Temp$Boulia
data_temp[[4]]=Australian_Temp$Cape_Otway
data_temp[[5]]=Australian_Temp$Gayndah
data_temp[[6]]=Australian_Temp$Gunnedah
data_temp[[7]]=Australian_Temp$Hobart
data_temp[[8]]=Australian_Temp$Robe


loc.years=list()
loc.years[[1]]=2:153
loc.years[[2]]=2:157
loc.years[[3]]=c(8:39,41:77,79:117)
loc.years[[4]]=2:148
loc.years[[5]]=13:116
loc.years[[6]]=2:96
loc.years[[7]]=15:120
loc.years[[8]]=2:128



data_FD=list()

for(j in 1:8){
  fun_data_S = data_temp[[j]]
  D = 20
  #basis = create.fourier.basis(rangeval = c(0, 1), nbasis = D)
  basis<-create.bspline.basis(rangeval=c(0,1),nbasis=D,norder=4)

  nas = which(is.na(fun_data_S$Days.of.accumulation.of.minimum.temperature))
  if(length(nas)>0)
    fun_data_S = fun_data_S[-nas, ]
  yy = unique(fun_data_S$Year)
  yy=yy[loc.years[[j]]]
  print(yy[1])
  print(yy[length(yy)])


  mat.S = matrix(0, D, length(yy))
  bads=1:length(yy)
  for (i in 1:length(yy)){
    aa = subset(fun_data_S, Year==yy[i])
    cc = aa$Minimum.temperature..Degree.C.
    bads[i]=length(cc)
    #print(length(cc))
    a = which(is.na(cc))
    if (length(a)>0){
      cc = cc[-which(is.na(cc))]
    }else{
      cc = cc
    }
    f_Obs = Data2fd(argvals=seq(0, 1, length = length(cc)) , cc, basisobj = basis)
    mat.S[, i] = f_Obs$coefs
  }
  fdata = fd(mat.S, basis)
  # note that the last year, has data only up to 6 months
  # therefore we remove it
  #fdata = fdata[-length(yy)]
  data_FD[[j]]=fdata
  plac=eval.fd(seq(0,1,by=1/360), fdata)

  print(which(bads<300))
  print(length(bads))
  data_temp[[j]]=plac
}

pdf("gyn.pdf")
  num=153
  plot(data_FD[[5]],col=rainbow(153),ylab="Temperature C",xlim=c(0.01,.99),xlab="",xaxt="n",lty=1)
  axis(side=1, at=c(0.01,3/12,6/12,9/12,11/12), labels=c("Jan","Mar","Jun","Sep","Dec"))
dev.off()


pdf("gyn.pdf")
plot_fd(eval.fd(seq(0,1,0.05),data_FD[[5]]),interactive = FALSE,val_axis_title = '',res_axis_title = '',FD_axis_title = '',showticklabels = F)
dev.off()

tmp <- linear_imputatation(temperature$Sydney)
pdf("sydney.pdf")
plot_fd(tmp[5:10],interactive = FALSE,
        val_axis_title = '',res_axis_title = '',FD_axis_title = '',showticklabels = F)
dev.off()
