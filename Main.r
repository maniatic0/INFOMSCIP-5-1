

# General -----------------------------------------------------------------

library("Hmisc");
library("corrplot");
library("caret");
library("tidyverse");
library("fastICA");
library("dplyr");
library("multcomp");
theme_set(theme_bw());

printf <- function(...) writeLines(sprintf(...));

# Load Data
gt_2011 = read.csv("./pp_gas_emission/gt_2011.csv");
gt_2012 = read.csv("./pp_gas_emission/gt_2012.csv");
gt_2013 = read.csv("./pp_gas_emission/gt_2013.csv");
gt_2014 = read.csv("./pp_gas_emission/gt_2014.csv");
gt_2015 = read.csv("./pp_gas_emission/gt_2015.csv");
gt_total = rbind(gt_2011, gt_2012, gt_2013, gt_2014, gt_2015);



# Part 2 ------------------------------------------------------------------


part2_analysis <- function(data)
{
  cols = colnames(data);
  rows = c("mean", "median", "sd");
  results = data.frame(matrix(ncol = length(cols), nrow = length(rows)));
  colnames(results) = cols;
  rownames(results) = rows;
  results[1,] =  sapply(data, mean);
  results[2,] =  sapply(data, median);
  results[3,] =  sapply(data, sd);
  results = rbind(results, sapply(data, quantile, c(0.01, 0.99)));
  results = rbind(results, sapply(data, range));
  rownames(results)[6] = "min";
  rownames(results)[7] = "max";
  return(results);
};

part2_perform_analysis <- function(data, name)
{
  res = part2_analysis(data);
  csvname = sprintf("./part2/summary_%s.csv", name);
  write.csv(res, csvname);
};

part2_perform_analysis(gt_2011, "2011");
part2_perform_analysis(gt_2012, "2012");
part2_perform_analysis(gt_2013, "2013");
part2_perform_analysis(gt_2014, "2014");
part2_perform_analysis(gt_2015, "2015");
part2_perform_analysis(gt_total, "Total");



# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
# From: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
flattenCorrMatrix <- function(cormat, pmat) {
  ut = upper.tri(cormat);
  res = data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  );
  
  return(res)
}


part2_corr_p_val = 0.01
part2_corr <- function(data, name)
{
  correlation = rcorr(as.matrix(data), type="spearman");
  flattened = flattenCorrMatrix(correlation$r, correlation$P);
  
  csvname = sprintf("./part2/corr_flattened_%s.csv", name);
  write.csv(flattened, csvname);
  
  # Insignificant correlation are crossed
  png(sprintf("./part2/corr_with_significance_%s.png", name));
  corrplot(correlation$r, type="upper", order="alphabet", 
           p.mat = correlation$P, sig.level = part2_corr_p_val,
           title = sprintf("%s's Corr, Significance Level %g", name, part2_corr_p_val),
           
           mar=c(0,0,1,0)
           );
  dev.off();
  
  return (correlation);
}

corr_2011 = part2_corr(gt_2011, "2011");
corr_2012 = part2_corr(gt_2011, "2012");
corr_2013 = part2_corr(gt_2013, "2013");
corr_2014 = part2_corr(gt_2014, "2014");
corr_2015 = part2_corr(gt_2015, "2015");
corr_total = part2_corr(gt_total, "Total");

corr_change = corr_2015$r - corr_2011$r;
corr_change_abs = apply(corr_change, 2, function(x) 1 - abs(x));
png("./part2/corr_change_2011_2015.png");
corrplot(corr_change, type="upper", order="alphabet",
         is.corr = FALSE,
         p.mat = corr_change_abs, sig.level = 0.8,
         title = "Corr Change 2011-2015, Change higher than 0.2 are crossed",
         mar=c(0,0,1,0), insig="label_sig"
);
dev.off();

top_change = which(corr_change_abs == min(corr_change_abs), arr.ind = TRUE);
corr_change_abs_copy = corr_change_abs;
corr_change_abs_copy[top_change] = 2;
top2_change = which(corr_change_abs_copy == min(corr_change_abs_copy), arr.ind = TRUE);
corr_change_abs_copy[top2_change] = 2;
top3_change = which(corr_change_abs_copy == min(corr_change_abs_copy), arr.ind = TRUE);
corr_change_abs_copy[top3_change] = 2;

cols_top_change = c("Var 1", "Var 2", "Value", "Abs Value");
results_top_change = data.frame(matrix(ncol = length(cols_top_change), nrow = 3));
colnames(results_top_change) = cols_top_change;
# 1
results_top_change[1, 1] = colnames(corr_change)[top_change[1, 1]];
results_top_change[1, 2] = colnames(corr_change)[top_change[1, 2]];
results_top_change[1, 3] = corr_change[top_change[1,1], top_change[1,2]];
results_top_change[1, 4] = abs(results_top_change[1, 3]);
# 2
results_top_change[2, 1] = colnames(corr_change)[top2_change[1, 1]];
results_top_change[2, 2] = colnames(corr_change)[top2_change[1, 2]];
results_top_change[2, 3] = corr_change[top2_change[1,1], top2_change[1,2]];
results_top_change[2, 4] = abs(results_top_change[2, 3]);
# 3
results_top_change[3, 1] = colnames(corr_change)[top3_change[1, 1]];
results_top_change[3, 2] = colnames(corr_change)[top3_change[1, 2]];
results_top_change[3, 3] = corr_change[top3_change[1,1], top3_change[1,2]];
results_top_change[3, 4] = abs(results_top_change[3, 3]);

# Save
write.csv(results_top_change, "./part2/top_corr_change_2011_2015.csv");


# Part 3 ------------------------------------------------------------------

predict_model <- function(linear_model, validation)
{
  # Make predictions
  predictions = linear_model %>% predict(validation);
  # Model performance
  # (a) Prediction error, RMSE
  rmse = RMSE(predictions, validation$NOX);
  # (b) R-square
  r2 = R2(predictions, validation$NOX);
  # (c) MAE
  mae = MAE(predictions, validation$NOX);
  # (d) Correlation
  corr = cor(predictions, validation$NOX, method = "spearman");
  
  printf("\nModel Performance");
  printf("RMSE=%g R2=%g MAE=%g Spearman=%g", rmse, r2, mae, corr);
  return (predictions);
}

train_model <- function(formula, training, validation, name, filename)
{
  sink(file = sprintf("./part3/%s_measures.txt", filename));
  linear_model = lm(
    formula, 
    data=training
  );
  # Model Summary
  printf("%s Model Summary", name);
  print(formula);
  print(summary(linear_model));
  predict_model(linear_model, validation);
  sink(file = NULL);
  return(linear_model);
};

# Paper Sets
training = subset(rbind(gt_2011, gt_2012), select = -CO);
validation = subset(gt_2013, select = -CO);
training_validation = rbind(training, validation);
test = subset(rbind(gt_2014, gt_2015), select = -CO);

# Phase 1

total_model_formula = NOX ~ AT + AP + AH + AFDP + GTEP + TIT + TAT + TEY + CDP

# a
train_model(total_model_formula, training, validation, "Phase 1 a", "phase_1_a");

# b
train_model(total_model_formula, training_validation, test, "Phase 1 b", "phase_1_b");

# Phase 2

# center, scale = Standardize data
# YeoJhonson = https://en.wikipedia.org/wiki/Power_transform#Yeo%E2%80%93Johnson_transformation
# ICA = https://en.wikipedia.org/wiki/Independent_component_analysis
# Spatial Sign = Project Everything to a Unit Sphere https://pubs.acs.org/doi/10.1021/ci050498u

nox_col = which(names(training) == "NOX");

our_model_formula = NOX ~ .;
preProcObj = preProcess(training[,-nox_col], method = c("center", "scale", "YeoJohnson", "ica", "spatialSign"), n.comp=8, outcome=training$NOX);

printModel = function(preprocObj, name){
  sink(file = name);
  printf("Our Model\n");
  printf("Summary");
  print(preProcObj);
  printf("\nMeans");
  print(preProcObj$mean);
  printf("\nStandard Deviation");
  print(preProcObj$std);
  printf("\nICA");
  print(preProcObj$ica);
  sink(file = NULL);
}
printModel(preProcObj, "./part3/our_model_info.txt");


prepareData <- function(preprocObj, data)
{
  new_data = cbind(predict(preprocObj, data[,-nox_col]), data$NOX);
  colnames(new_data)[dim(new_data)[2]] = "NOX";
  new_data
};

weight_visualize <- function(preprocObj, linear_model, training, name, filename)
{
  
  png(sprintf("./part3/%s_weights.png", filename));
  barplot(
    t(as.matrix(linear_model$coefficients)), 
    las=2, 
    main=sprintf("Weight Visualization%s", name), 
    sub = "Note that the non-ica parameters are spatial sign processed",
    ylim=c(min(linear_model$coefficients)*1.2, max(linear_model$coefficients)*1.2)
    );
  dev.off();
  
  orig = preprocObj$ica$K %*% preprocObj$ica$W %*% our_model$coefficients[11:18];
  orig = t(orig);
  colnames(orig) = colnames(training)[1:9];
  
  png(sprintf("./part3/%s_ica_weights.png", filename));
  barplot(
    orig, 
    las=2, 
    main=sprintf("ICA Weight Visualization%s", name), 
    sub="Original Features",
    ylim=c(min(orig)*1.2, max(orig)*1.2)
  );
  dev.off();
  
  # See https://topepo.github.io/caret/variable-importance.html
  imp = varImp(linear_model, scale=FALSE);
  png(sprintf("./part3/%s_importance.png", filename));
  barplot(
    t(as.matrix(
      imp[order(imp$Overall,decreasing = TRUE),])
      ),
    names.arg = rownames(imp)[order(imp$Overall,decreasing = TRUE)],
    las=2, 
    main=sprintf("Importance Visualization%s", name),
    ylim=c(0, max(imp$Overall)*1.2),
    sub = "The absolute value of the t-statistic for each model parameter"
    );
  dev.off();
  
};

training_prep = prepareData(preProcObj, training);
validation_prep = prepareData(preProcObj, validation);
test_prep = prepareData(preProcObj, test);
our_model = train_model(our_model_formula, training_prep, validation_prep, "Our Model", "our_model");

# Weight Visualization

weight_visualize(preProcObj, our_model, training, " Our Model", "our_model");

# Probes

probes_model <- function (linear_model, preprocObj, test, n, filename)
{
  sink(file = sprintf("./part3/%s_probes_%d.txt", filename, n));
  probes = sample_n(test, n);
  probes_prep = prepareData(preprocObj, probes);
  predictions = predict_model(linear_model, probes_prep);
  printf("\nProbes");
  data_to_show = cbind(probes, predictions, probes$NOX - predictions);
  colnames(data_to_show)[length(colnames(data_to_show))-1] = "Prediction";
  colnames(data_to_show)[length(colnames(data_to_show))] = "Error";
  print(data_to_show);
  sink(file = NULL);
  write.csv(data_to_show, sprintf("./part3/%s_probes_%d.csv", filename, n));
  
  for (i in 1:n) {
    
    multi = linear_model$coefficients * cbind(1, probes_prep[i, 1:(ncol(probes_prep)-1)]);
    
    multi = cbind(multi, predictions[i], probes$NOX[i]);
    colnames(multi)[1] = "(Intercept)";
    colnames(multi)[19:20] = cbind("Prediction", "Real");
    
    png(sprintf("./part3/%s_probe_%d_weights.png", filename, i));
    barplot(
      as.matrix(multi), 
      las=2, 
      main=sprintf("Probe %d Weight Visualization", i), 
      sub = "Note that the non-ica parameters are spatial sign processed",
      ylim=c(min(multi)*1.2, max(multi)*1.2)
    );
    dev.off();
    
    orig = preprocObj$ica$K %*% preprocObj$ica$W %*% t(as.matrix(multi[11:18]));
    orig = t(orig);
    colnames(orig) = colnames(test)[1:9];
    
    png(sprintf("./part3/%s_probe_%d_ica_weights.png", filename, i));
    barplot(
      orig, 
      las=2, 
      main=sprintf("Probe %d ICA Weight Visualization", i), 
      sub="Original Features",
      ylim=c(min(orig)*1.2, max(orig)*1.2)
    );
    dev.off();
  }
}

probes_model(our_model, preProcObj, test, 2, "our_model");

# Phase 3

measure_model <- function(linear_model, validation)
{
  # Make predictions
  predictions = linear_model %>% predict(validation);
  # Model performance
  # (a) Prediction error, RMSE
  rmse = RMSE(predictions, validation$NOX);
  # (b) R-square
  r2 = R2(predictions, validation$NOX);
  # (c) MAE
  mae = MAE(predictions, validation$NOX);
  # (d) Correlation
  corr = cor(predictions, validation$NOX, method = "spearman");
  
  printf("\nModel Performance");
  printf("RMSE=%g R2=%g MAE=%g Spearman=%g", rmse, r2, mae, corr);
  return (c(rmse, r2, mae, corr));
}


create_blocks <- function(data, block_size)
{
  n = nrow(data);
  r = rep(1:block_size,each=ceiling(n/block_size))[1:n];
  split(data,r)
}


compare_models <- function(training, validation, block_size, filename)
{
  trainData = training;
  blocks = create_blocks(validation, block_size);
  
  cnames = c("Type", "RMSE", "R2", "MAE", "Spearman");
  measures = data.frame(matrix(ncol = length(cnames), nrow = 0));
  colnames(measures) = cnames;
  
  preds = c();
  preds = as.data.frame(preds);
  
  sink(file = sprintf("./part3/%s_phase3_%d_blocks_log.txt", filename, block_size));
  printf("\n\nOnline Tests");
  for (i in 1:length(blocks)) {
    block = blocks[[i]];
    
    preProcObj = preProcess(trainData[,-nox_col], method = c("center", "scale", "YeoJohnson", "ica", "spatialSign"), n.comp=8, outcome=trainData$NOX);
    
    tData = prepareData(preProcObj, trainData);
    vData = prepareData(preProcObj, block);
    
    our_model = lm(
      our_model_formula, 
      data=tData
    );
    
    preds = rbind(preds, cbind("ours_online", our_model %>% predict(vData)));
    printf("\nOur Model Online");
    measures = rbind(measures, c("ours_online", measure_model(our_model, vData)));
    colnames(measures) = cnames;
    
    original_model = lm(
      total_model_formula, 
      data=trainData
    );
    
    preds = rbind(preds, cbind("original_online",original_model %>% predict(block)));
    printf("\nOriginal Model Online");
    measures = rbind(measures, c("original_online", measure_model(original_model, block)));
    colnames(measures) = cnames;
    
    trainData = rbind(trainData, block);
  }
  
  {
    printf("\n\nOffline Tests");
    preProcObj = preProcess(training[,-nox_col], method = c("center", "scale", "YeoJohnson", "ica", "spatialSign"), n.comp=8, outcome=training$NOX);
    
    tData = prepareData(preProcObj, training);
    vData = prepareData(preProcObj, validation);
    
    our_model = lm(
      our_model_formula, 
      data=tData
    );
    
    preds = rbind(preds, cbind("ours_offline", our_model %>% predict(vData)));
    printf("\nOur Model Offline");
    measures = rbind(measures, c("ours_offline", measure_model(our_model, vData)));
    colnames(measures) = cnames;
    
    original_model = lm(
      total_model_formula, 
      data=training
    );
    
    preds = rbind(preds, cbind("original_offline",original_model %>% predict(validation)));
    printf("\nOriginal Model Offline");
    measures = rbind(measures, c("original_offline", measure_model(original_model, validation)));
    colnames(measures) = cnames;
  }
  
  
  sink(file = NULL);
  
  sink(file = sprintf("./part3/%s_phase3_%d_blocks.txt", filename, block_size));
  
  printf("T Test Between Pearson Correlations of Original Model and Our Model Online");
  print(
    t.test(
      as.numeric(measures[which(measures$Type == "original_online"),]$Spearman), 
      as.numeric(measures[which(measures$Type == "ours_online"),]$Spearman)
    )
  );
  
  colnames(preds) = c("Type", "Predictions");
  preds$Predictions = as.numeric(preds$Predictions);
  preds$Type = factor(preds$Type);
  
  printf("\nANOVA Test Between Predictions of the Original Model Online, Our Model Online, Original Model Offline, and Our Model Offilne");
  anova = aov(Predictions ~ Type, data=preds);
  print(anova)
  print(summary(anova))
  # Check https://en.wikipedia.org/wiki/Tukey%27s_range_test
  print(summary(glht(anova, linfct = mcp(Type = "Tukey"))))
  
  sink(file = NULL);
  
  write.csv(measures, sprintf("./part3/%s_phase3_%d_blocks_measurments.csv", filename, block_size));
  
  measures
}

compare_models(training, validation, 10, "1_5");
compare_models(training_validation, test, 20, "6");

