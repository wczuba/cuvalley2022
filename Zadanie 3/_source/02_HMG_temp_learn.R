#Lib ----

library("lubridate")
library("hexbin")
library("latticeExtra")
library("dplyr")
library("pheatmap")
library("caTools")
library("gridExtra")
library(randomForestSRC)

#Load adn subset learning dataset  ----

setwd("E:/R_hack2")
smc<-readRDS("learning_dset.rds")
smc<-na.omit(smc)

smc<-smc[,c(1:4,6,7,8,12,24,25,26,27,28,30,31,32,36,48,49,50,51,52,54,55,56,60,72,73:78)]
colnames(smc)[31]<-"temp_zuz"
colnames(smc)[32]<-"dT"
colnames(smc)[29]<-"godzina"
colnames(smc)[30]<-"dzien_tyg"
colnames(smc)[33]<-"odst_pom"
smc$last_temp<-lag(smc$temp_zuz,1)
smc$godzina<-as.character(floor(smc$godzina))
smc$dzien_tyg<-as.character(floor(smc$dzien_tyg))
smc[,29]<-as.factor(smc[,29])
smc[,30]<-as.factor(smc[,30])
smc<-na.omit(smc)
smca<-smc

for (i in c(1:27))
{
smca[,i]<-c(0,diff(smc[,i]))
}

smc<-smc[smc$'odst_pom'<4000,]
smc<-smc[smc$'dT'>=-25&smc$'dT'<=25,]
smc<-smc[smc$'temp_zuz'>=1270,]

smca<-smca[smca$'odst_pom'<4000,]
smca<-smca[smca$'dT'>=-25&smca$'dT'<=25,]
smca<-smca[smca$'temp_zuz'>=1270,]




#SMC Long term Temp Parallel plots ----

smc$temp_kat<-cut(smc$'temp_zuz',breaks=c(quantile(smc$'temp_zuz', probs = seq(0, 1, by = 0.2),na.rm=T)))

smcg<-na.omit(smc) %>% group_by(temp_kat,`godzina`,`dzien_tyg`) %>% summarise(across(everything(), list(mean)))
parallel(~smcg[,4:30]|temp_kat, data = smcg, layout=c(5,1),panel=function(x,y,...)
{
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})

smcg<-na.omit(smc) %>% group_by(temp_kat,`dzien_tyg`) %>% summarise(across(everything(), list(mean)))
parallel(~smcg[,3:29]|temp_kat, data = smcg, layout=c(5,1),panel=function(x,y,...)
{
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})

smcg<-na.omit(smc) %>% group_by(temp_kat,`godzina`) %>% summarise(across(everything(), list(mean)))
parallel(~smcg[,3:29]|temp_kat, data = smcg, layout=c(5,1),panel=function(x,y,...)
{
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})

smcg<-na.omit(smc) %>% group_by(temp_kat) %>% summarise(across(everything(), list(mean)))
parallel(~smcg[,2:28], data = smcg,lwd=1,panel=function(x,y,...)
{
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})

#SMCA Short term dT Parallel plots ----

smca$dT_kat<-cut(smca$'dT',breaks=c(quantile(smca$'dT', probs = seq(0, 1, by = 0.2),na.rm=T)))

smcag<-na.omit(smca) %>% group_by(dT_kat,`godzina`,`dzien_tyg`) %>% summarise(across(everything(), list(mean)))
parallel(~smcag[,4:30]|dT_kat, data = smcag, layout=c(5,1),panel=function(x,y,...)
{
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})

smcag<-na.omit(smca) %>% group_by(dT_kat,`dzien_tyg`) %>% summarise(across(everything(), list(mean)))
parallel(~smcag[,3:29]|dT_kat, data = smcag, layout=c(5,1),panel=function(x,y,...)
{
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})

smcag<-na.omit(smca) %>% group_by(dT_kat,`godzina`) %>% summarise(across(everything(), list(mean)))
parallel(~smcag[,3:29]|dT_kat, data = smcag, layout=c(5,1),auto.key=T,panel=function(x,y,...)
{
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})

smcag<-na.omit(smca) %>% group_by(dT_kat) %>% summarise(across(everything(), list(mean)))
parallel(~smcag[,2:28], data = smcag,lwd=2,auto.key=T,panel=function(x,y,...)
  {
  panel.parallel(x,y,...)
  panel.abline(h=c(9.5,18.5),lty=2,lwd=2)
})



#Train model sampling learn/test data ----


sample_list<-sample(c(1:nrow(smc)))
smc<-smc[sample_list,]
smca<-smca[sample_list,]

train_factor=0.8
train_count=floor(train_factor*nrow(smc))


#SMC Long term temp model ----

model.rf1 <- rfsrc(temp_zuz ~ ., data=smc[c(1:train_count),c(1:27,29,30,31)],ntree = 5000, nodedepth=100,proximity =F,distance = T, forest.wt = F,importance = c("none", "permute", "random", "anti")[3])

test_result1<-predict(model.rf1,newdata=smc[c(train_count:nrow(smc)),c(1:27,29,30)])
#model.rf1.cor=sd(smc[c(1:train_count),c(31)])/sd(test_result1$predicted)
#test_result1$predicted=test_result1$predicted*model.rf1.cor

hexbinplot(test_result1$predicted~smc[c(train_count:nrow(smc)),31],xbin=32,aspect=1,trans = sqrt, inv = function(x) x^2,
           ylab="Predicted",xlab="Real",panel=function(x,y,...){
             panel.hexbinplot(x,y,...)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.25, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.50, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.75, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.95, segments = 50,)
           })

#SMC Short term dT model ----

model.rf2 <- rfsrc(dT ~ ., data=smca[c(1:train_count),c(1:27,29,30,32,31)],ntree = 5000, nodedepth=100,proximity =F,distance = T, forest.wt = F,importance = c("none", "permute", "random", "anti")[3])

smca_t_l<-smca[c(train_count:nrow(smca)),c(1:27,29,30,31)]
smca_t_l[,30]<-test_result1$predicted
  
test_result2<-predict(model.rf2,newdata=smca_t_l)
model.rf2.cor=sd(smca[c(1:train_count),c(32)])/sd(test_result2$predicted)
test_result2$predicted=test_result2$predicted*model.rf2.cor

hexbinplot(test_result2$predicted~smca[c(train_count:nrow(smca)),32],xbin=32,aspect=1,trans = sqrt, inv = function(x) x^2,
           ylab="Predicted",xlab="Real",panel=function(x,y,...){
             panel.hexbinplot(x,y,...)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.25, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.50, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.75, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.95, segments = 50,)
           })


hexbinplot((test_result2$predicted+smca[c(train_count:nrow(smca)),34])~smca[c(train_count:nrow(smca)),31],xbin=32,aspect=1,trans = sqrt, inv = function(x) x^2,
           ylab="Predicted",xlab="Real",panel=function(x,y,...){
             panel.hexbinplot(x,y,...)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.25, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.50, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.75, segments = 50,)
             panel.ellipse(x,y,...,col="blue",lwd=2,alpha=0.8,level = 0.95, segments = 50,)
           })





res<-cbind((test_result2$predicted+smca[c(train_count:nrow(smca)),34]),test_result1$predicted,smca[c(train_count:nrow(smca)),c(31,28)],smca[c(train_count:nrow(smca)),34])
res<-res[order(res$'1_czas'),]
xyplot(res[,1]+res[,2]+res[,3]~res[,4],type="spline")

res<-res[order(res[,3]),]

xyplot(res[,1]~c(1:nrow(res)),ylim=c(1280,1325),ylab="Temp. żuzla [st.C]",xlab="n. of obs.",panel=function(x,y,...)
{
  panel.xyplot(x,y,...,type="p",col="gray")
  panel.xyplot(x,runquantile(y,200,probs=0.25),...,type="l",col="black",lwd=2)
  panel.xyplot(x,runquantile(y,200,probs=0.75),...,type="l",col="black",lwd=2)
  #panel.xyplot(x,runquantile(y,200,probs=0.5),...,type="l",col="black",lwd=2)
})+as.layer(xyplot(res[,3]~c(1:nrow(res)),type="l",lwd=3,col="red"),x.same = T)




xyplot(res[,2]~c(1:nrow(res)),ylim=c(1280,1325),panel=function(x,y,...)
{
  panel.xyplot(x,y,...,type="p",col="gray")
})+as.layer(xyplot(res[,3]~c(1:nrow(res)),type="l",lwd=2,col="red"),x.same = T)

xyplot(res[,5]~res[,3])


#SMC base test

#model.rf <- rfsrc(temp_zuz ~ ., data=smc[c(1:train_count),c(1:27,29,30,31)],ntree = 10, block.size =10,proximity = F, distance = F, forest.wt = F)

#test_result<-predict(model.rf,newdata=smc[c(1:train_count),c(1:27,29,30,31)])
#hexbinplot(test_result$predicted~smc[c(1:train_count),c(31)],xbin=60,aspect=1,
#           panel=function(x,y,...){
#             panel.hexbinplot(x,y,...)
#             panel.ellipse(x,y,...,col="red",lwd=2)
#           })
#
#mean(test_result$predicted/smc[c(1:train_count),c(31)])


