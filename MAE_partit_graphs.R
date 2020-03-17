library(dplyr)
library(lme4)
library(ggplot2)
library(caret)

cook_comp_data<-read.csv("/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Moisture_Uptake_Data.csv")
cook_comp_data<-na.omit(cook_comp_data)
cook_comp_data$Moisture_Avg<-cook_comp_data$Moisture_Avg*100
head(cook_comp_data)

MAE_test<-vector(length = 99)
MAESDA_test<-vector(length = 99)

for(i in 1:88){ #above 88 you get a warning saying the new data has 0 rows.
mup_index<-createDataPartition(cook_comp_data$Genotype, p=(i/100), list = FALSE)
mup_test<-cook_comp_data[-mup_index,]
mup_train<-cook_comp_data[mup_index,]

control<-trainControl(method = "cv", number = 10)
metric<-"MAE"

set.seed(7)
fit.rf<-train(Moisture_Avg~
                Protein_As_is+ 
                Starch_As_is+ 
                Fiber_As_is+ 
                Fat_As_is+
                Ash_As_is+
                Starch_As_is:Fat_As_is+
                Starch_As_is:Ash_As_is+ 
                Protein_As_is:Fiber_As_is, 
              data = mup_train, method = "rf", metric = metric, trControl = control)


predictions<-predict(fit.rf,mup_test)

MAE_test[i]<-mean(abs(mup_test$Moisture_Avg-predictions)) #MAE of test

MAESDA_test[i]<-(MAE_test-fit.rf$results$MAE[3])/fit.rf$results$MAESD[3]
print(i)
}
lm(MAE_test[1:88]~c(1:88))
lm(MAESDA_test[1:88]~c(1:88))

plot(MAE_test)
abline(1.89435,-0.00801)
abline(-1.17122,0.04655)
plot(x = NA, type = "n", ylim = c(0,6), xlim = c(0,88))
abline(1.89435,-0.00801)
abline(-1.17122,0.04655)
locator()
