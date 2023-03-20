############################################################
##             K-fold cross-validation，分层分组          ##
############################################################
cv_group <- function(k, label, seed)  # label：标签向量，且赋值为0或-1和1
{
  set.seed(seed)
  temp <- data.frame(obs=c(1:length(label)), label=label)
  temp_classone <- temp[temp$label==0 | temp$label==-1,]
  temp_classtwo <- temp[temp$label==1,]
  n_classone <- rep(1:k, ceiling(nrow(temp_classone)/k))[1:nrow(temp_classone)]
  n_classtwo <- rep(1:k, ceiling(nrow(temp_classtwo)/k))[1:nrow(temp_classtwo)]
  temp_classone$k <- sample(n_classone, nrow(temp_classone), replace=F)
  temp_classtwo$k <- sample(n_classtwo, nrow(temp_classtwo), replace=F)
  x <- 1:k
  cvlist <- list()
  cvlist <- lapply(x, function(x) 
  {
    obs <- c(temp_classone$obs[temp_classone$k==x], temp_classtwo$obs[temp_classtwo$k==x])
    obs <- obs[order(obs)]
  })
}
############################################################
##                     Hold-out，分层分组                 ##
############################################################
ho_group <- function(train_num, label, seed)  # # label：标签向量，且赋值为0或-1和1
{
  set.seed(seed)
  temp <- data.frame(obs=c(1:length(label)), label=label)
  sample_rate <- train_num / nrow(temp)  # 抽样比
  temp_class_one <- temp[temp$label==0 | temp$label==-1,]
  temp_class_two <- temp[temp$label==1,]
  sample_num_one <- ceiling(nrow(temp_class_one) * sample_rate)
  sample_num_two <- train_num - sample_num_one
  train_index_one <- sample(temp_class_one$obs, sample_num_one, replace=F)
  train_index_two <- sample(temp_class_two$obs, sample_num_two, replace=F)
  train_index <- c(train_index_one, train_index_two)
  train_index <- train_index[order(train_index)]
  test_index <- setdiff(temp$obs, train_index)
  ho_list <- list()
  ho_list[[1]] <- train_index
  ho_list[[2]] <- test_index
  return(ho_list)
}
############################################################
##                discrimination measure                  ##
############################################################
measure_disc <- function(true_label, pre_prob)  # 只计算AUC
{
  library(ROCR)
  temp <- data.frame(true_label, pre_prob)
  pred <- prediction(temp$pre_prob, temp$true_label)
  AUC <- performance(pred,"auc")@y.values[[1]]
  return(AUC)
}
############################################################
##                  calibration measure                   ##
############################################################
measure_cali <- function(true_label, pre_prob, num_bin=10, plot_title="calibration plot")
{
  library(PredictABEL)
  Data_temp <- data.frame(true_label, pre_prob)
  calibration <- plotCalibration(data=Data_temp, cOutcome=1, predRisk=Data_temp$pre_prob, groups=num_bin, plottitle=plot_title)  # outcome=true label
  chi_squre <- calibration$Chi_square
  p_hosleme <- calibration$p_value
  MSE <- sum((Data_temp$true_label - Data_temp$pre_prob)^2)/nrow(Data_temp)
  result <- calibration$Table_HLtest
  ECE <- sum(abs(result[, 2] - result[ ,3]) * (result[, 1] / nrow(Data_temp)))
  MCE <- max(abs(result[, 2] - result[, 3]))
  measure_result <- c(chi_squre=chi_squre, p_hosleme=p_hosleme, MSE=MSE, ECE=ECE, MCE=MCE)
  return(measure_result)
}
############################################################
##                    platt calibration                   ##
############################################################
platt_sacling <- function(pre_prob, true_label)  # y与true_label相同
{
  y = true_label
  temp_platt <- data.frame(pre_prob, y, true_label)
  num_pos <- length(temp_platt$true_label[temp_platt$true_label==1])
  num_neg <- length(temp_platt$true_label[temp_platt$true_label==0])
  temp_platt$y[temp_platt$y==1] <- (num_pos + 1) / (num_pos + 2)
  temp_platt$y[temp_platt$y==0] <- 1 / (num_neg + 2)
  temp_platt$y <- as.factor(temp_platt$y)
  platt <- glm(y~pre_prob, data=temp_platt, family=binomial)
  A <- platt$coefficients[2]
  B <- platt$coefficients[1]
  para <- c(A, B)
  names(para) <- c("para_A", "para_B")
  return(para)
}
############################################################
##                   isotonic calibration                 ##
############################################################
iso_reg <- function(pre_prob, true_label)
{
  temp_iso <- data.frame(pre_prob, true_label)
  temp_iso <- temp_iso[order(temp_iso$pre_prob),]
  iso <- isoreg(x=temp_iso$pre_prob, y=temp_iso$true_label)
  temp_cutpoint <- data.frame(p_pre=iso$x, p_cali=iso$yf)
  cut <- temp_cutpoint$p_cali[!duplicated(temp_cutpoint$p_cali)]
  cut_mat <- matrix(0, nrow=length(cut), ncol=3)
  r <- 1
  for(i in cut)
  {
    temp <- temp_cutpoint[temp_cutpoint$p_cali==i,]
    cut_mat[r,] <- c(min(temp$p_pre), max(temp$p_pre), i)
    r <- r+1
  }
  cut_mat <- as.data.frame(cut_mat)
  colnames(cut_mat) <- c("Lower", "Upper", "p_cali")
  return(cut_mat)
}
