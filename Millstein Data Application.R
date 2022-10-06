library(nnet); library(mcp.project);library(MASS);library(readxl); library(caTools);library(pROC)

# Multinomial logistic regression model function, no additional covariates. Performs DLR likelihood ratio test
MCL<-function(data, blocks, blocksize, ref.grp, plot.param){
  ## data is a nx2 dataset, where the first column is the biomarker data and the second is the binary outcome
  ## blocks is the number of intervals to bin the data into
  ## blocksize is the number of subjects to be in each interval
  ## ref.grp is the interval number to be taken as the reference group in the multinomial logistic regression model
  ## plot.param is a TRUE/FALSE condition where the ROC curve is plotted for the biomarker if 'TRUE'
  
  # Add column names and sort data from smallest to largest
  colnames(data)<-c('marker', 'class')
  ds1sort<-data[order(data$marker),]
  
  
  # Create equally sized intervals of size BLOCKS
  ds1sort$group<-c(rep(1:(blocks-1), each=blocksize), rep(blocks, times=blocksize-(blocks*blocksize-nrow(ds1sort))))
  # Choose reference group
  ds1sort$categoryRef<-ifelse(ds1sort$group==ref.grp, 0, ds1sort$group)
  # Fit multinomial model with Y as main effect
  f1<-multinom(ds1sort$categoryRef ~ ds1sort$class, trace=FALSE)
  # Fit intercept-only multinomial model
  f2<-multinom(ds1sort$categoryRef ~ 1, trace=FALSE)
  # Perform likelihood ratio test and retain test statistic and p-value
  an_f1f2_pval<-anova(f1, f2)[2,7]
  an_f1f2_TS<-anova(f1, f2)[2,6]
  
  
  # Calculate Wald score for individual parameters
  z <- summary(f1)$coefficients/summary(f1)$standard.errors
  # 2-tailed z test for Wald score
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  p<-round(p, digits=4)
  # Summarize model parameters and results into table
  f1sum<-cbind(summary(f1)$coefficients, summary(f1)$standard.errors, exp(coef(f1)), z, p)
  # Calculate referent group DLR
  DLRJ<-1/((1+sum(exp(f1sum[,1]+f1sum[,2])))/(1+sum(f1sum[,5])))
  # Calculate all other DLRj values
  DLRj<-f1sum[,6]*DLRJ
  
  
  f1sum<-cbind(f1sum, DLRj, DLRJ, an_f1f2_TS, an_f1f2_pval)
  f1sum<-round(f1sum, digits=3)
  colnames(f1sum)<-c("Beta0j", "Beta1j", "SE(Beta0j)", "SE(Beta1j)", "exp(Beta0j)", "exp(Beta1j)", 
                     "Beta0j.Z", "Beta1j.Z", "Beta0j.P-val", "Beta1j.P-val", 'DLRj', 'DLRJ', 'LRT TS', 'LRT p-val')
  
  ## Plot ROC curve
  colAUC(ds1sort$group, ds1sort$class, plotROC = plot.param)
  if(plot.param=='T'){
    abline(a=0, b=1)
  }
  return(f1sum)
}
# Multinomial logistic regression model function, with two additional main effect covariates. Performs DLR likelihood ratio test
MCL_age_site<-function(data, blocks, blocksize, ref.grp){
  ## data is a nx4 dataset, where the first column is the biomarker data, the second is the binary outcome, and 3 and 4 are categorical main effect covariates
  ## blocks is the number of intervals to bin the data into
  ## blocksize is the number of subjects to be in each interval
  ## ref.grp is the interval number to be taken as the reference group in the multinomial logistic regression model
  ## plot.param is a TRUE/FALSE condition where the ROC curve is plotted for the biomarker if 'TRUE'
  
  # Add column names and sort data from smallest to largest
  colnames(data)<-c('marker', 'class', 'covar1', 'covar2')
  ds1sort<-data[order(data$marker),]
  
  
  # Create equally sized intervals of size BLOCKS
  ds1sort$group<-c(rep(1:(blocks-1), each=blocksize), rep(blocks, times=blocksize-(blocks*blocksize-nrow(ds1sort))))
  # Choose reference group
  ds1sort$categoryRef<-ifelse(ds1sort$group==ref.grp, 0, ds1sort$group)
  # Fit multinomial model with Y as main effect and Z1 and Z2 as additional covariates
  f03<-multinom(ds1sort$categoryRef ~ ds1sort$class + factor(ds1sort$covar1) + factor(ds1sort$covar2), MaxNWts =2000, trace=F)
  # Fit multinomial model without Y as a main effect
  f_null<-multinom(ds1sort$categoryRef ~ factor(ds1sort$covar1) + factor(ds1sort$covar2), trace=F)
  # Perform likelihood ratio test between the null and full model. Retain test statistic and p-value
  an_f1f2_pval<-anova(f03, f_null)[2,7]
  an_f1f2_TS<-anova(f03, f_null)[2,6]
  
  # Returns p-value from DLR likelihood ratio test
  return(an_f1f2_pval)
}
# Calculates DLR from contingence table, and gives interval ranges of the biomarker
DLR<-function(data, blocks, blocksize){
  ## data is a nx2 dataset, where the first column is the biomarker data and the second is the binary outcome
  ## blocks is the number of intervals to bin the data into
  ## blocksize is the number of subjects to be in each interval
  
  # Add column names to data and sort biomarker data from smallest to largest
  colnames(data)<-c('marker', 'class')
  ds1sort<-data[order(data$marker),]
  
  
  # Create equally sized intervals of size BLOCKS
  ds1sort$group<-c(rep(1:(blocks-1), each=blocksize), rep(blocks, times=blocksize-(blocks*blocksize-nrow(ds1sort))))
  
  # Obtain interval endpoints of each biomarker interval
  mat.all<-NULL
  min.all<-NULL
  max.all<-NULL
  for(i in 1:blocks){
    data_sub<-subset(ds1sort, ds1sort$group==i)
    min.range<-min(data_sub$marker, na.rm=T)
    max.range<-max(data_sub$marker, na.rm=T)
    min.all<-c(min.all, min.range)
    max.all<-c(max.all, max.range)
  }
  
  # Create data matrix
  mat.all<-cbind(t(table(ds1sort$class, ds1sort$group)), min.all, max.all)
  
  # Calculate DLR from contingency table
  N<-sum(mat.all[,1])
  P<-sum(mat.all[,2])
  DLR<-(mat.all[,2]/P)/(mat.all[,1]/N)
  
  # Format data matrix and add column names
  mat.all<-cbind(mat.all[,1:2], DLR, mat.all[,3:4])
  colnames(mat.all)<-c("# Controls", "# Cases", "DLRj", "Interval min.", "Interval max.")
  
  return(mat.all)
}

# Apply Cochran-Armitage test for trend with 10 intervals in the biomarker data
CA_test_10<-function(data, weights, blocksize){
  ## data is a nx2 dataset, where the first column is the biomarker data and the second is the binary outcome
  ## weights is the 10-element vector of weights, traditional or nontraditonal, to use in the conventional, not-modified C-A test for trend
  ## blocksize is the number of observations in each of the 10 intervals. Should be n/10 to the nearest integer.
  
  # Sort biomarker data from smallest to largest
  ds1sort<-data[order(data[,1]),]
  
  # Create equally sized groups using 10 intervals
  ds1sort$group<-c(rep(1:(10-1), each=blocksize), rep(10, times=blocksize-(10*blocksize-nrow(ds1sort))))
  
  # Calculate contingency table on the data, and calculate the marginal values
  t1<-table(ds1sort[,2], ds1sort$group)
  R1<-sum(t1[1,])
  R2<-sum(t1[2,])
  Ci<-c(sum(t1[,1]),sum(t1[,2]),sum(t1[,3]),sum(t1[,4]),sum(t1[,5]),sum(t1[,6]),sum(t1[,7]),sum(t1[,8]),sum(t1[,9]),sum(t1[,10]))

  # Calculate the C-A test statistic  
  CA_TS<-sum(weights*(t1[1,]*R2-t1[2,]*R1))
  
  # Calculate the partial contributions of the variance of the test statistic
  V1<-(weights[1]*Ci[1])*(weights[2]*Ci[2] + weights[3]*Ci[3] + weights[4]*Ci[4] + weights[5]*Ci[5] +
                            weights[6]*Ci[6] + weights[7]*Ci[7] + weights[8]*Ci[8] + weights[9]*Ci[9] + weights[10]*Ci[10])
  V2<-(weights[2]*Ci[2])*(weights[3]*Ci[3] + weights[4]*Ci[4] + weights[5]*Ci[5] + weights[6]*Ci[6] +
                            weights[7]*Ci[7] + weights[8]*Ci[8] + weights[9]*Ci[9] + weights[10]*Ci[10])
  V3<-(weights[3]*Ci[3])*(weights[4]*Ci[4] + weights[5]*Ci[5] + weights[6]*Ci[6] + weights[7]*Ci[7] +
                            weights[8]*Ci[8] + weights[9]*Ci[9] + weights[10]*Ci[10])
  V4<-(weights[4]*Ci[4])*(weights[5]*Ci[5] + weights[6]*Ci[6] + weights[7]*Ci[7] + weights[8]*Ci[8] +
                            weights[9]*Ci[9] + weights[10]*Ci[10])
  V5<-(weights[5]*Ci[5])*(weights[6]*Ci[6] + weights[7]*Ci[7] + weights[8]*Ci[8] + weights[9]*Ci[9] + weights[10]*Ci[10])
  V6<-(weights[6]*Ci[6])*(weights[7]*Ci[7] + weights[8]*Ci[8] + weights[9]*Ci[9] + weights[10]*Ci[10])
  V7<-(weights[7]*Ci[7])*(weights[8]*Ci[8] + weights[9]*Ci[9] + weights[10]*Ci[10])
  V8<-(weights[8]*Ci[8])*(weights[9]*Ci[9] + weights[10]*Ci[10])
  V9<-(weights[9]*Ci[9])*(weights[10]*Ci[10])
  Vi<-sum(V1, V2, V3, V4, V5, V6, V7, V8, V9)
  
  # Calculate variance of the test statistic
  CA_Var<-((R1*R2)/nrow(data))*(sum(((weights**2)*Ci)*(sum(t1[,])-Ci)) - (2*Vi))
  
  # Calculate variance-covariance matrix
  N<-length(ds1sort$group)
  c_vec<-colSums(t1)
  sigma_mat<-diag(c_vec*(N-c_vec))
  sigma_mat[1,2:10]<--c_vec[1]*c_vec[2:10]
  sigma_mat[2:10,1]<--c_vec[1]*c_vec[2:10]
  
  sigma_mat[2,3:10]<--c_vec[2]*c_vec[3:10]
  sigma_mat[3:10,2]<--c_vec[2]*c_vec[3:10]
  
  sigma_mat[3,4:10]<--c_vec[3]*c_vec[4:10]
  sigma_mat[4:10,3]<--c_vec[3]*c_vec[4:10]
  
  sigma_mat[4,5:10]<--c_vec[4]*c_vec[5:10]
  sigma_mat[5:10,4]<--c_vec[4]*c_vec[5:10]
  
  sigma_mat[5,6:10]<--c_vec[5]*c_vec[6:10]
  sigma_mat[6:10,5]<--c_vec[5]*c_vec[6:10]
  
  sigma_mat[6,7:10]<--c_vec[6]*c_vec[7:10]
  sigma_mat[7:10,6]<--c_vec[6]*c_vec[7:10]
  
  sigma_mat[7,8:10]<--c_vec[7]*c_vec[8:10]
  sigma_mat[8:10,7]<--c_vec[7]*c_vec[8:10]
  
  sigma_mat[8,9:10]<--c_vec[8]*c_vec[9:10]
  sigma_mat[9:10,8]<--c_vec[8]*c_vec[9:10]
  
  sigma_mat[9,10]<--c_vec[9]*c_vec[10]
  sigma_mat[10,9]<--c_vec[9]*c_vec[10]
  
  sigma_mat<-((rowSums(t1)[1]*rowSums(t1)[2])/N)*sigma_mat
  
  return(list(CA_TS, sigma_mat))
}

# Import dataset and subjec to training data
GSE132342 <- read_excel("GSE132342.xlsx")

train<-subset(GSE132342, GSE132342$site=='DOV' | GSE132342$site=='AOC' | GSE132342$site=='VAN' |
                GSE132342$site=='SEA' | GSE132342$site=='UKO' | GSE132342$site=='WMH' |
                GSE132342$site=='RTR' | GSE132342$site=='POC' | GSE132342$site=='BRO' |
                GSE132342$site=='AOV' | GSE132342$site=='CNI' | GSE132342$site=='USC' |
                GSE132342$site=='LAX' | GSE132342$site=='GER' | GSE132342$site=='TRI1')
train<-subset(train, train$Stage!=8) #Exclude subjects misssing outcome variable
train<-data.frame(train)

# Recode outcome variable
train$Stage<-ifelse(train$Stage==1, 0, 1)
table(train$Stage)

# Apply DLR likelihood ratio test for all genes, adjusting for age, study, using 10 intervals in the data
LRT_all<-NULL
for(i in 14:526){
  LRT<-MCL_age_site(data.frame(train[,c(i,8, 11, 12)]), blocks = 10, 253, 10)
  LRT_all<-c(LRT_all, LRT)
}
# Apply Benjamini-Hochberg procedure to control FDR at 5%
adj_pval<-fdr(LRT_all, method='BH', q=.05)
table(adj_pval$Pvals$rejected) #249 rejections
sig_genes<-which(adj_pval$Pvals$rejected==T)+13 #gives column number of dataset for significant genes

# Apply modified test for trend to the significant genes
CAtest_all<-NULL
TS_diff_all<-NULL
pval_all<-NULL
for(i in sig_genes){
  
  mean_cont<-mean(train[which(train$Stage==0),i])
  mean_case<-mean(train[which(train$Stage==1),i])
  
  # Keep conditional biomarker distributions such that cases have on average, larger mean than the controls
  if(mean_cont>mean_case){
    train$Stage2<-ifelse(train$Stage==1, 0, 1)
  } else{
    train$Stage2<-train$Stage
  }
  
  # Apply each of the individual C-A tests using monotone and non-monotone weights to obtain T(w) and var-covar matrices
  CAtest_trad<-(CA_test_10(data.frame(train[,i], train$Stage2), 10:1, 253))
  CAtest_nt<-(CA_test_10(data.frame(train[,i], train$Stage2), c(9,7,5,3,1,0,2,4,6,8), 253))
  
  # Calculate numerator of test statistic
  difftruth<-CAtest_trad[[1]]+CAtest_nt[[1]]
  
  # Calculate variance and covariance components of test statistic
  weight_trad<-matrix(data=c(10:1), nrow = 10)
  weight_nt<-c(9,7,5,3,1,0,2,4,6,8)
  var_trad<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_trad
  var_nt<-t(weight_nt)%*%CAtest_nt[[2]]%*%weight_nt
  covar_T<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_nt
  
  # Calculate and retain test statistic
  TS_diff<-difftruth/sqrt(var_nt + var_trad + (2*covar_T))
  TS_diff_all<-c(TS_diff_all, TS_diff)
  
  # Obtain p-value for test statistic
  if(TS_diff<0){
    p_val1<-pnorm(TS_diff, lower.tail = T)*2
  } else{
    p_val1<-pnorm(TS_diff, lower.tail = F)*2
  }
  
  pval_all<-c(pval_all, p_val1)
}

# Split the biomarkers by biomarker type
trad_genes<-sig_genes[which(TS_diff_all > 0 & pval_all < 0.05)] # Column # of traditional significant genes
trad_genes<-trad_genes[-200] # Something happened in switching the class labels of this biomarker making its test statistic positive instead of negative. It should fall into the nontraditional group, so remove it from traditional and add to nontraditional
nt_genes<-c(sig_genes[which(TS_diff_all < 0 & pval_all < 0.05)], 510) # Column # of nontraditional significant genes
ind_genes1<-sig_genes[which(TS_diff_all < 0 & pval_all > 0.05)]
ind_genes2<-sig_genes[which(TS_diff_all > 0 & pval_all > 0.05)]
ind_genes<-c(ind_genes1, ind_genes2) # Column # of indeterminate significant genes


# Conduct the analysis again, but use test for AUC at the step to identify informative biomarkers
auc_test_all<-NULL
AUC_all<-NULL
for(i in 14:526){
  # Use the test for AUC based off the Mann-Whitney U test
  auc_test<-wilcox.test(subset(train, train$Stage==1)[,i], subset(train, train$Stage==0)[,i])$`p.val`
  AUC_all<-c(AUC_all, colAUC(train[,i], train$Stage))
  auc_test_all<-c(auc_test_all, auc_test)
}

# Control FDR
table(ifelse(auc_test_all <0.05, 1, 0)) # 213 rejections
adj_pval_auc<-fdr(auc_test_all, method='BH', q=.05)
table(adj_pval_auc$Pvals$rejected) # 166  rejections
sig_genes_auc<-which(adj_pval_auc$Pvals$rejected==T)+13 # Column #s of significant genes by AUC

# Apply modified test for trend to AUC significant biomarkers
CAtest_all<-NULL
TS_diff_all<-NULL
pval_all<-NULL
for(i in sig_genes_auc){
  
  mean_cont<-mean(train[which(train$Stage==0),i])
  mean_case<-mean(train[which(train$Stage==1),i])
  
  if(mean_cont>mean_case){
    train$Stage2<-ifelse(train$Stage==1, 0, 1)
  } else{
    train$Stage2<-train$Stage
  }
  
  CAtest_trad<-(CA_test_10(data.frame(train[,i], train$Stage2), 10:1, 253))
  CAtest_nt<-(CA_test_10(data.frame(train[,i], train$Stage2), c(9,7,5,3,1,0,2,4,6,8), 253))
  
  difftruth<-CAtest_trad[[1]]+CAtest_nt[[1]]
  
  weight_trad<-matrix(data=c(10:1), nrow = 10)
  weight_nt<-c(9,7,5,3,1,0,2,4,6,8)
  var_trad<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_trad
  var_nt<-t(weight_nt)%*%CAtest_nt[[2]]%*%weight_nt
  covar_T<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_nt
  
  (TS_diff<-difftruth/sqrt(var_nt + var_trad + (2*covar_T)))
  TS_diff_all<-c(TS_diff_all, TS_diff)
  
  if(TS_diff<0){
    p_val1<-pnorm(TS_diff, lower.tail = T)*2
  } else{
    p_val1<-pnorm(TS_diff, lower.tail = F)*2
  }
  
  pval_all<-c(pval_all, p_val1)
}


trad_genes_auc<-sig_genes_auc[which(TS_diff_all > 0 & pval_all < 0.05)] # Column # of traditional significant genes
nt_genes_auc<-sig_genes_auc[which(TS_diff_all < 0 & pval_all < 0.05)] # Column # of nontraditional significant genes
ind_genes1<-sig_genes_auc[which(TS_diff_all < 0 & pval_all > 0.05)]
ind_genes2<-sig_genes_auc[which(TS_diff_all > 0 & pval_all > 0.05)]
ind_genes_auc<-c(ind_genes1, ind_genes2) # Column # of indeterminate significant genes

`%!in%` <- Negate(`%in%`)
sig_genes[which(sig_genes %!in% sig_genes_auc)] #genes that AUC missed


# Conduct the analysis again, but use test for ROC curve length at the step to identify informative biomarkers
# Get case and control sample sizes
length_cont<-length(train[which(train$Stage==0),8])
length_case<-length(train[which(train$Stage==1),8])

# Apply length-based test
pval_all_length<-NULL
for(i in 14:526){
  mu0<-mean(train[which(train$Stage==0),i])
  mu1<-mean(train[which(train$Stage==1),i])
  sd0<-sd(train[which(train$Stage==0),i])
  sd1<-sd(train[which(train$Stage==1),i])
  
  # Length function
  l1<-function(x){
    sqrt(dnorm(x, mu0, sd0)^2 + (dnorm(x, mu1, sd1))^2) #assumes binormal model for markers
  }
  
  # Estimate of length, Equation 2 in Bantis
  l1_est<-integrate(l1, lower=min(train[,i]), upper = max(train[,i]))
  
  
  # Calculate test statistic and p-value
  Z<-((sqrt(2)*4*length_cont*length_case)/(length_cont + length_case))*(l1_est$value - sqrt(2))
  pval_all_length<-c(pval_all_length, pchisq(Z, 2, lower.tail = F))
  
}

# Control FDR
table(ifelse(pval_all_length <0.05, 1, 0)) #140 rejections
adj_pval_length<-fdr(pval_all_length, method='BH', q=.05)
table(adj_pval_length$Pvals$rejected) #73  rejections
sig_genes_length<-which(adj_pval_length$Pvals$rejected==T)+13 # Column #s of significant genes by ROC length

# Table of p-value and adjusted p-values for each search strategy (LR test, AUC test, length test) for each biomarker
View(cbind(1:513, LRT_all, adj_pval$Pvals$adjusted.pvals, AUC_all,auc_test_all, adj_pval_auc$Pvals$adjusted.pvals, pval_all_length, adj_pval_length$Pvals$adjusted.pvals))

# Apply modified test for trend to significant genes by Length test
CAtest_all<-NULL
TS_diff_all<-NULL
pval_all<-NULL
for(i in sig_genes_length){
  
  mean_cont<-mean(train[which(train$Stage==0),i])
  mean_case<-mean(train[which(train$Stage==1),i])
  
  if(mean_cont>mean_case){
    train$Stage2<-ifelse(train$Stage==1, 0, 1)
  } else{
    train$Stage2<-train$Stage
  }
  
  CAtest_trad<-(CA_test_10(data.frame(train[,i], train$Stage2), 10:1, 253))
  CAtest_nt<-(CA_test_10(data.frame(train[,i], train$Stage2), c(9,7,5,3,1,0,2,4,6,8), 253))
  
  difftruth<-CAtest_trad[[1]]+CAtest_nt[[1]]
  
  weight_trad<-matrix(data=c(10:1), nrow = 10)
  weight_nt<-c(9,7,5,3,1,0,2,4,6,8)
  var_trad<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_trad
  var_nt<-t(weight_nt)%*%CAtest_nt[[2]]%*%weight_nt
  covar_T<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_nt
  
  (TS_diff<-difftruth/sqrt(var_nt + var_trad + (2*covar_T)))
  TS_diff_all<-c(TS_diff_all, TS_diff)
  
  if(TS_diff<0){
    p_val1<-pnorm(TS_diff, lower.tail = T)*2
  } else{
    p_val1<-pnorm(TS_diff, lower.tail = F)*2
  }
  
  pval_all<-c(pval_all, p_val1)
}

# Identify traditional, nontraditional, and indeterminate genes
trad_genes_length<-sig_genes_length[which(TS_diff_all > 0 & pval_all < 0.05)] #column # of van with traditional sig. genes
nt_genes_length<-sig_genes_length[which(TS_diff_all < 0 & pval_all < 0.05)] #column # of van with nontraditional sig. genes
ind_genes1<-sig_genes_length[which(TS_diff_all < 0 & pval_all > 0.05)]
ind_genes2<-sig_genes_length[which(TS_diff_all > 0 & pval_all > 0.05)]
ind_genes_length<-c(ind_genes1, ind_genes2)

table(adj_pval_length$Pvals$rejected, adj_pval$Pvals$rejected)
#length as rows, LRT as column
#they agreed on rejecting 195
#54 that LRT thought were important that length did not
#78 that length thought were important that LRT did not
#186 they agreed on failing to reject

# Look at the ROC curves of the 54 that Length did not select, but LR test did
discord_pairs<-ifelse(adj_pval_length$Pvals$rejected==F & adj_pval$Pvals$rejected==T, 1, 0)
length_miss<-which(discord_pairs==1)+13 #column numbers of the biomarkers the AUC test missed compared 
par(mfrow=c(3,3))
for(i in length_miss){
  colAUC(train[,i], train[,8], plotROC = T)
  abline(a=0,b=1)
}


# Identify the most significant nontraditional gene when the likelihood ratio test was used
nt_genes[which(cbind(nt_genes, adj_pval$Pvals$adjusted.pvals[nt_genes - 13])[,2] == min(cbind(nt_genes, adj_pval$Pvals$adjusted.pvals[nt_genes - 13])[,2]))]
colnames(train)[38]
# G25 most nontraditional
# ROC curve of continuous data
marker_roc<-pROC::roc(predictor=train[,38], response=train$Stage, quiet=T)
spec<-1-marker_roc$specificities
sens<-marker_roc$sensitivities
plot(spec, sens, type='s')
abline(a=0,b=1, lty=2)


# Identify the most significant traditional gene when the likelihood ratio test was used
most_sig_trad<-trad_genes[which(cbind(trad_genes, LRT_all[trad_genes - 13])[,2] == min(cbind(trad_genes, LRT_all[trad_genes - 13])[,2]))]
colnames(train)[most_sig_trad]
# There are multiple equally significant traditional genes, pick the one with highest AUC
# Obtain AUC for the most significant traditional genes
most_trad_auc<-NULL
for(i in (most_sig_trad)){
  # ROC curve of continuous data
  marker_roc<-pROC::roc(predictor=train[,i], response=train$Stage, quiet=T)
  spec<-1-marker_roc$specificities
  sens<-marker_roc$sensitivities
  most_trad_auc<-c(most_trad_auc, marker_roc$auc)
}
# G429 strongest traditional by LRT and AUC


# Plot these two biomarkers using the ROC curve and DLR curve
par(mfrow=c(2,2))
# Do traditional biomarker first
# ROC curve
ds1sort<-train[order(train$G429),]
ds1sort$group<-c(rep(1:(10-1), each=253), rep(10, times=253-(10*253-nrow(ds1sort))))
fit1<-glm(ds1sort$Stage ~ ds1sort$group, family = binomial)
pred429<-predict(fit1)
roc429<-roc(ds1sort$Stage, pred429, plot = F)
# Calculate confidence intervals for ROC curve
ci.sp.obj <- ci.sp(roc429, sensitivities=roc429$sensitivities, boot.n=1000)
plot(roc429, ylim=c(0,1), main='ROC Curve', cex.lab=1.2, cex.axis=1.2) # restart a new plot
plot(ci.sp.obj, type="shape", col="grey", ylim=c(0,1))
abline(a=1,b=-1, lwd=1, lty=1)
# Obtain DLR through MLR model
fit1<-MCL(train[,c(442,8)], 10, 253, 10, 'F')
DLR_val<-c(fit1[,6]/fit1[1,11], 1/fit1[1,11])

# DLR curve
# Obtain interval endpoints
bins<-DLR(train[,c(442,8)], 10, 253)
plot(seq(bins[1,4], bins[1,5], length.out = 2), rep(DLR_val[1], times=2), type='l', 
     ylim=c(0,5.5), xlim=c(bins[1,4], bins[10,5]), ylab = 'DLRj', xlab = 'log-standardized marker', main='DLR Curve', sub='Gene: PCDH9', xaxt='n', cex.lab=1.2, cex.axis=1.2)
abline(h=1, lty=2, col='grey')
lines(seq(bins[2,4], bins[2,5], length.out = 2), rep(DLR_val[2], times=2))
lines(seq(bins[3,4], bins[3,5], length.out = 2), rep(DLR_val[3], times=2))
lines(seq(bins[4,4], bins[4,5], length.out = 2), rep(DLR_val[4], times=2))
lines(seq(bins[5,4], bins[5,5], length.out = 2), rep(DLR_val[5], times=2))
lines(seq(bins[6,4], bins[6,5], length.out = 2), rep(DLR_val[6], times=2))
lines(seq(bins[7,4], bins[7,5], length.out = 2), rep(DLR_val[7], times=2))
lines(seq(bins[8,4], bins[8,5], length.out = 2), rep(DLR_val[8], times=2))
lines(seq(bins[9,4], bins[9,5], length.out = 2), rep(DLR_val[9], times=2))
lines(seq(bins[10,4], bins[10,5], length.out = 2), rep(DLR_val[10], times=2))
axis(1, at=c(bins[,4], bins[10,5]), labels=round(c(bins[,4], bins[10,5]), digits=1), cex.axis=1)

# Repeat for nontraditional biomarker
# ROC curve
ds1sort<-train[order(train$G25),]
ds1sort$group<-c(rep(1:(10-1), each=253), rep(10, times=253-(10*253-nrow(ds1sort))))
fit1<-glm(ds1sort$Stage ~ ds1sort$group, family = binomial)
pred364<-predict(fit1)
roc364<-roc(ds1sort$Stage, pred364, plot = F)
ci.sp.obj <- ci.sp(roc364, sensitivities=roc364$sensitivities, boot.n=1000)
plot(roc364, ylim=c(0,1), main='ROC Curve', cex.lab=1.2, cex.axis=1.2) 
plot(ci.sp.obj, type="shape", col="grey", ylim=c(0,1))
abline(a=1,b=-1, lwd=1)

# DLR curve
fit1<-MCL(train[,c(38,8)], 10, 253, 10, 'F')
DLR_val<-c(fit1[,6]/fit1[1,11], 1/fit1[1,11])
bins<-DLR(train[,c(38,8)], 10, 253)
plot(seq(bins[1,4], bins[1,5], length.out = 2), rep(DLR_val[1], times=2), type='l',
     ylim=c(0,5.5), xlim=c(bins[1,4], bins[10,5]), ylab = 'DLRj', xlab = 'log-standardized marker', main='DLR Curve', sub='Gene: PTEN', xaxt='n', cex.lab=1.2, cex.axis=1.2)
abline(h=1, lty=2, col='grey')
lines(seq(bins[2,4], bins[2,5], length.out = 2), rep(DLR_val[2], times=2))
lines(seq(bins[3,4], bins[3,5], length.out = 2), rep(DLR_val[3], times=2))
lines(seq(bins[4,4], bins[4,5], length.out = 2), rep(DLR_val[4], times=2))
lines(seq(bins[5,4], bins[5,5], length.out = 2), rep(DLR_val[5], times=2))
lines(seq(bins[6,4], bins[6,5], length.out = 2), rep(DLR_val[6], times=2))
lines(seq(bins[7,4], bins[7,5], length.out = 2), rep(DLR_val[7], times=2))
lines(seq(bins[8,4], bins[8,5], length.out = 2), rep(DLR_val[8], times=2))
lines(seq(bins[9,4], bins[9,5], length.out = 2), rep(DLR_val[9], times=2))
lines(seq(bins[10,4], bins[10,5], length.out = 2), rep(DLR_val[10], times=2))
axis(1, at=c(bins[,4], bins[10,5]), labels=round(c(bins[,4], bins[10,5]), digits=1), cex.axis=1)




# Validate the 6 nontraditional and top 3 traditional genes from discovery
# Pull the testing dataset
test<-subset(GSE132342, GSE132342$site=='TRI2' | GSE132342$site=='NCO' | GSE132342$site=='MAY' |
               GSE132342$site=='SRF' | GSE132342$site=='HAW' | GSE132342$site=='POL')
test<-subset(test, test$Stage!=8) # Exclude subjects with missing outcome variable
test<-data.frame(test)

# Re-code outcome variable
table(test$Stage)
test$Stage<-ifelse(test$Stage==1, 0, 1)
table(test$Stage)

# Apply modified test for trend to 9 genes of interest
CAtest_all<-NULL
TS_diff_all<-NULL
pval_all<-NULL
auc_all<-NULL
for(i in c(442, 43, 272, 38, 510, 264, 472, 118, 134)){ # <- column numbers of the 9 genes of interest
  
  marker_roc<-pROC::roc(predictor=test[,i], response=test$Stage, quiet=T)
  spec<-1-marker_roc$specificities
  sens<-marker_roc$sensitivities
  auc_all<-c(auc_all, marker_roc$auc)
  
  # Figure 8 (right) plot
  if(i == 264){
    plot(spec, sens, type='s', main='HIST1H2BE\n Validation', xlab='1-Spec.', ylab='Sens.', lwd=2)
    abline(a=0,b=1, lty=2)
  }
  
  mean_cont<-mean(test[which(test$Stage==0),i])
  mean_case<-mean(test[which(test$Stage==1),i])
  
  if(mean_cont>mean_case){
    test$Stage2<-ifelse(test$Stage==1, 0, 1)
  } else{
    test$Stage2<-test$Stage
  }
  
  # When i=264, since ROC starts below and ends above 45-degree line, switch labels
  if(i==264){
    test$Stage2<-1-test$Stage2
  }
  
  CAtest_trad<-(CA_test_10(data.frame(test[,i], test$Stage2), 10:1, 100))
  CAtest_nt<-(CA_test_10(data.frame(test[,i], test$Stage2), c(9,7,5,3,1,0,2,4,6,8), 100))
  
  difftruth<-CAtest_trad[[1]]+CAtest_nt[[1]]
  
  weight_trad<-matrix(data=c(10:1), nrow = 10)
  weight_nt<-c(9,7,5,3,1,0,2,4,6,8)
  var_trad<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_trad
  var_nt<-t(weight_nt)%*%CAtest_nt[[2]]%*%weight_nt
  covar_T<-t(weight_trad)%*%CAtest_trad[[2]]%*%weight_nt
  
  TS_diff<-difftruth/sqrt(var_nt + var_trad + (2*covar_T))
  TS_diff_all<-c(TS_diff_all, TS_diff)
  
  if(TS_diff<0){
    p_val1<-pnorm(TS_diff, lower.tail = T)*2
  } else{
    p_val1<-pnorm(TS_diff, lower.tail = F)*2
  }
  
  pval_all<-c(pval_all, p_val1)
}

sig_genes_test<-c(442, 43, 272, 38, 510, 264, 472, 118, 134)

# Identify the gene classification of the 9 genes of interest
trad_genes_test<-sig_genes_test[which(TS_diff_all > 0 & pval_all < 0.05)] #column # of van with traditional sig. genes
nt_genes_test<-sig_genes_test[which(TS_diff_all < 0 & pval_all < 0.05)] #column # of van with nontraditional sig. genes
ind_genes1<-sig_genes_test[which(TS_diff_all < 0 & pval_all > 0.05)]
ind_genes2<-sig_genes_test[which(TS_diff_all > 0 & pval_all > 0.05)]
ind_genes_test<-c(ind_genes1, ind_genes2)


sum(trad_genes_test %in% trad_genes) #3 of the 3 of traditional validated
sum(nt_genes_test %in% nt_genes) #1 of the 6 nontraditional validated

