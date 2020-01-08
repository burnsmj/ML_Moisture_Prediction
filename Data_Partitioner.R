kernel_comp_data<-read.csv("/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Matlab_Tests/Matlab_all_nonred_traits_full.csv", check.names = F)
head(kernel_comp_data)

mup_index<-caret::createDataPartition(kernel_comp_data$Moisture_Avg, p=0.90, list = FALSE)
mup_test<-kernel_comp_data[-mup_index,]
mup_train<-kernel_comp_data[mup_index,]

write.csv(mup_train, "/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Matlab_Tests/Matlab_all_nonred_traits_train.csv", row.names = F)
write.csv(mup_test, "/Users/michael/Desktop/Grad_School/Research/Datasets/Machine Learning/Matlab_Tests/Matlab_all_nonred_traits_test.csv", row.names = F)

#Added this note as a test