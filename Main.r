

# General -----------------------------------------------------------------

library("Hmisc");
library("corrplot");
library("caret");
library("tidyverse");
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

train_model <- function(training, validation, name)
{
  sink(file = sprintf("./part3/%s.txt", name));
  linear_model = lm(
    NOX ~ AT + AP + AH + AFDP + GTEP + TIT + TAT + TEY + CDP, 
    data=training
  );
  # Model Summary
  printf("%s Model Summary", name);
  print(summary(linear_model));
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
  sink(file = NULL);
};

# Sets
training = rbind(gt_2011, gt_2012);
validation = gt_2013;
training_validation = rbind(training, validation);
test = rbind(gt_2014, gt_2015);

# Phase 1

# a
train_model(training, validation, "Phase 1 a");
train_model(training_validation, test, "Phase 1 b");






