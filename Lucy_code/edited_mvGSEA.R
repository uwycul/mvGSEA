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
    nullmodel <- edited_speedglm_L1(y~1, drop.unused.levels=FALSE)
  }
  
  ###-------------there ARE confounders-----------
  ##----need to make two models: FULL and partial
  ###-----FULL: prod(v) + prod(conf)
  ###-----PARTIAL: prod(conf)
  
  if(!is.null(conf_col)){
    fullmodel <- edited_speedglm_L1(f,data=d, drop.unused.levels=FALSE )
    partialmodel <- edited_speedglm_L1(partial, data=d, drop.unused.levels=FALSE)
    nullmodel <- edited_speedglm_L1(y~1, drop.unused.levels=FALSE)
    #next, perform likelihood ratio test to pick the full or partial model to use (NOTE: the lrtest function doesn't work with the speedglm function (only glm), so I hacked the code and made my own test that is compatile with speedglm)
    q <- 2*abs(fullmodel$logLik-partialmodel$logLik)
    df <- abs(fullmodel$df-partialmodel$df)
    pval <- pchisq(q, df, lower.tail=FALSE)
  if (pval >0.05){
    model <- fullmodel
  } else{
    model <- partialmodel
  }
  }
 
  
  #---------------------------------------------------------------------------
#comparing model with null model
  #perform chisq----get overall pvalue
  #anova(nullmodel, model, test="Chisq")
  q_2 <- 2*abs(model$logLik-nullmodel$logLik)
  df_2 <- abs(model$df-nullmodel$df)
  pval_ <- pchisq(q_2, df_2, lower.tail=FALSE) #overall pvalue from anova
  
  



  #prediction
 predict_model <- as.vector(predict.speedglm(model, data.frame(d), type = "response"))
ones_indices <- which(y==1)
predict_model_ones <- predict_model[ones_indices]
 
 #if prediction_method="sensitivity", rounding
 if (prediction_method=="sensitivity"){
   sorted_predict_model <- sort(predict_model_ones, decreasing=FALSE)
 prediction_sensitivity_cutoff <- sorted_predict_model[as.integer(9*length(sorted_predict_model)/10)] #retrieving the indexed value at 90th percentile
predict_model_fit <- as.integer(predict_model>prediction_sensitivity_cutoff)
 
   
   #if prediction_method="custom"  , rounding
 } else if (prediction_method=="custom"){
 predict_model_fit <- as.integer(predict_model>prediction_cutoff) #might change the 0.5 cutoff based on the balance of categories
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
data.frame(pval_, oddsratio, auc, TP, FP, FN, TN)

  
}
