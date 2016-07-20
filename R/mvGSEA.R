#uses the edited_speed_glm_L1 and edited_L1_speedglm.wfit to perform a Wald Test and to predict the outcome
#will return pval and odds ratio as a data frame

mvGSEA <- function(y,X,v_col,conf_col, prediction_method,prediction_cutoff, f, d, partial)
{

  #v is now a matrix of only the variables of interest (called v_col)
  #conf is now a matrix of only confounders of interest (called conf_col)
  
  #----------perform glm depending on the confounders-----------
  
  ###-----------------no confounders (all multiplication of selected v's)

  if(is.null(conf_col)){  #NO CONFOUNDERS
    model <- edited_speedglm_L1(f, data=d, drop.unused.levels=FALSE)
    partialmodel <- edited_speedglm_L1(y~1, drop.unused.levels=FALSE) #NULL model
  }
  
  ###-------------there ARE confounders-----------
  ##----need to make two models: FULL and partial
  ###-----FULL: prod(v) + prod(conf)
  ###-----PARTIAL: prod(conf)
  
  if(!is.null(conf_col)){
    model <- edited_speedglm_L1(f,data=d, drop.unused.levels=FALSE );# 
    partialmodel <- edited_speedglm_L1(partial, data=d, drop.unused.levels=FALSE); 
  }
  
  q <- 2*abs(model$logLik-partialmodel$logLik)
  df <- abs(model$df-partialmodel$df)
  pval <- pchisq(q, df, lower.tail=FALSE)

  #---------------------------------------------------------------------------


  #prediction
 predict_model <- as.vector(predict.speedglm(model, data.frame(d), type = "response"))
 roc <- pROC::roc(y, predict_model)
 methods <- coords(roc, "best", ret=c("threshold"), best.method="youden")

 #if prediction_method="sensitivity", determining cutoff
if (prediction_method=="sensitivity"){
  ones_indices <- which(y==1)
   predict_model_ones <- predict_model[ones_indices]
  sorted_predict_model <- sort(predict_model_ones, decreasing=FALSE)
prediction_sensitivity_cutoff <- sorted_predict_model[as.integer(.1*length(sorted_predict_model))] #retrieving the indexed value at 90th percentile
predict_model_fit <- as.integer(predict_model>prediction_sensitivity_cutoff)
 } else if (prediction_method=="threshold"){
   prediction_threshold_cutoff <- methods[1]
   predict_model_fit <- as.integer(predict_model>prediction_threshold_cutoff)
 } 
 
  #make a contingency table and calculate odds ratio
y_ <- factor(y, levels=c(0,1))
 predictions <- factor(predict_model_fit, levels = c(0,1))
 
 contingency_table <- xtabs(~y_+predictions, drop.unused.levels=FALSE)

  
  
 TP <- contingency_table[2,2]
FP <- contingency_table[1,2]
FN <- contingency_table[2,1]
TN <- contingency_table[1,1]
  
  oddsratio<- fisher.test(contingency_table)$estimate

  
  #AUC from ROC
  
pr <- ROCR::prediction(predict_model_fit, y) #need ROCR package
 prf <- ROCR::performance(pr, measure="tpr", x.measure="fpr")
 auc <- ROCR::performance(pr, measure="auc")
 auc <- auc@y.values[[1]]

 
 
  #return a row
 data.frame(pval, oddsratio, auc, TP, FP, FN, TN)


  
}
