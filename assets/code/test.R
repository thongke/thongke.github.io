rm(list = ls())
library(ElemStatLearn)
library(foreign)
library(randomForest)
library(xtable)
library(e1071) # naiveBayes(.), svm(.), tune.svm(.)
library(rpart) # tree
# library(neuralnet) # neural net
library(nnet) # neural net
# library(tree)
# library("gbm")
library("glmnet") # for regulated regression
# library(class) # for k-NN knn(.)
# install.packages("RGtk2")
library("rattle")
rattle()
library(arules)

diabete = read.dta("diabete.dta",convert.factors=FALSE)
names(diabete) = c("id","sex","married","job","edu","income","alcohol","smoke","height","weight","systolic","diastolic","glucose","age","riskscore","diabetic","waist","hip")
str(diabete,give.attr = FALSE)
exclude = c(which(colnames(diabete)==c("id")),
			which(colnames(diabete)==c("glucose")),
			which(colnames(diabete)==c("diabetic")))
exclude

diabete$sex = as.factor(diabete$sex)
diabete$married = as.factor(diabete$married)
diabete$job = as.factor(diabete$job)
diabete$edu = as.factor(diabete$edu)
diabete$alcohol = as.factor(diabete$alcohol)
diabete$smoke = as.factor(diabete$smoke)
diabete$diabetic = as.factor(diabete$diabetic)

diabete_bin = diabete[,-c(1,13)]
diabete_reg = diabete[,-c(1,16)]
diabete_un = diabete[,-c(13,16)]

summary(diabete$glucose)
# diabete$glucose[diabete$glucose==72] == 7.2
table(diabete$diabetic)
plot(glucose~diabetic,data=diabete)

# Missing
sapply(diabete, function(x) sum(is.na(x)))

N = nrow(diabete)
n.cols = ncol(diabete)

set.seed(2014)
idx = sample(N,size=2*N/3,replace=F)
trainset = diabete[idx,]
testset = diabete[-idx,]

# Remove missing data
trainset_narm = na.omit(trainset)
testset_narm = na.omit(testset)
nrow(trainset_narm);nrow(testset_narm)
# Impute using mean/median
trainset_rf = na.roughfix(trainset)
testset_rf = na.roughfix(testset)
nrow(trainset_rf);nrow(testset_rf)

# model_classb = diabetic~sex+married+job+edu+income+alcohol+smoke+height+weight+systolic+diastolic+age+riskscore+waist+hip
model_class = diabetic~sex+married+job+edu+income+alcohol+smoke+height+weight+systolic+diastolic+age+riskscore+waist+hip
model_reg = glucose~sex+married+job+edu+income+alcohol+smoke+height+weight+systolic+diastolic+age+riskscore+waist+hip

#=============================================
#=            Naive Bayes
#=============================================
# Can handle both categorical and numeric input, 
# but output must be categorical
modelnaive <- naiveBayes(model_class, data=trainset)
prediction <- predict(modelnaive, testset[,-exclude])
table(prediction, testset$diabetic)

#=====================================
#=            SVM
#=====================================
tune <- tune.svm(model_class, data=trainset, 
                   gamma=10^(-6:-1), 
                   cost=10^(1:4))
summary(tune)
modelsvm <- svm(model_class, 
               data=trainset, 
               method="C-classification", 
               kernel="sigmoid", 
               probability=T, 
               gamma=0.000001, 
               cost=10)
summary(modelsvm)
prediction <- predict(modelsvm, testset_narm[,-exclude])
table(prediction,testset_narm$diabetic)

#================================================
#=            Neural network
#================================================
require(nnet, quietly=TRUE)
set.seed(199)
modelnn <- nnet(glucose ~ .,
    data=diabete_reg,
    size=10, linout=TRUE, skip=TRUE, MaxNWts=10000, trace=FALSE, maxit=100)
# Print the results of the modelling.

cat(sprintf("A %s network with %d weights.\n",
    paste(crs$nnet$n, collapse="-"),
    length(crs$nnet$wts)))
cat(sprintf("Inputs: %s.\n",
    paste(crs$nnet$coefnames, collapse=", ")))
cat(sprintf("Output: %s.\n",
    names(attr(crs$nnet$terms, "dataClasses"))[1]))
cat(sprintf("Sum of Squares Residuals: %.4f.\n",
    sum(residuals(crs$nnet) ^ 2)))
cat("\n")
print(summary(modelnn))
cat('\n')

#=============================================
#=            MARS model
# ============================================
library(mda)
library(earth)

modelmars_rm <- earth(model_reg, trainset_narm,trace=1,nk=20,degree=1,penalty=2,thresh=0.001,minspan=1,fast.k=0,fast.beta=0,pmethod="backward")
modelmars_rminter <- earth(model_reg, trainset_narm,trace=1,nk=20,degree=2,penalty=3,thresh=0.001,minspan=1,fast.k=0,fast.beta=0,pmethod="backward")
modelmars_rf <- earth(model_reg, trainset_rf,trace=1,nk=20,degree=1,penalty=2,thresh=0.001,minspan=1,fast.k=0,fast.beta=0,pmethod="backward")
modelmars_rfinter <- earth(model_reg, trainset_rf,trace=1,nk=20,degree=2,penalty=3,thresh=0.001,minspan=1,fast.k=0,fast.beta=0,pmethod="backward")

summary(modelmars_rm,digits=4)
summary(modelmars_rminter,digits=4)
summary(modelmars_rf,digits=4)
summary(modelmars_rfinter,digits=4)
# make predictions for the test data
# and compute mean squared error
prediction <- predict(modelmars_rf,newdata=testset_rf[,-exclude])
a = mean((testset_rf$glucose - prediction)^2)
prediction <- predict(modelmars_rfinter,newdata=testset_rf[,-exclude])
a1 = mean((testset_rf$glucose - prediction)^2)
prediction <- predict(modelmars_rm,newdata=testset_narm[,-exclude])
a2 = mean((testset_narm$glucose - prediction)^2)
prediction <- predict(modelmars_rminter,newdata=testset_narm[,-exclude])
a3 = mean((testset_narm$glucose - prediction)^2)
print(c(a,a1,a2,a3),digits=4)

print(xtable(modelmars_rf$coefficients,label="tab:modelmars_rf$coefficients",caption="tab:modelmars_rf$coefficients", align = c("@{\\extracolsep{\\fill}}","lc")),booktabs=TRUE, caption.placement = "top", table.placement = "htb",tabular.environment = "tabular*",width = '0.75\\textwidth',comment=FALSE,latex.environments="")
# plot diagnostics of the model plot() plotd()
# plot the line of best fit for the training data
pdf("modelmars_rm.pdf")
plotmo(modelmars_rm,cex=1)
dev.off()
# importance of input variables evimp(modelmars_rm)
#---------------  End of MARS  ----------------#