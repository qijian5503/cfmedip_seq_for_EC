library(caret)
load("medip_esca_samp146_phen.RData")
table(phen_esca$Group)
valid_sample=createDataPartition(phen_esca$Group, p = 0.2,times = 1)
valid_sample=phen_esca$Sample[valid_sample[[1]]]
save(valid_sample,file="medip_esca_peak_samp146_valid_samp.RData")

load("medip_esca_peak_samp146_valid_samp.RData")
train_test=phen_esca[!phen_esca$Sample %in% valid_sample,]
train_test_samp=train_test$Sample
subsample=createDataPartition(train_test$Group, p = 0.8,times = 100)
for (i in 1:100) {
    train_samp=train_test_samp[subsample[[i]]]
    test_samp=train_test_samp[!train_test_samp %in% train_samp]
    subsample[[i]]=list(train_samp=train_samp,test_samp=test_samp)
}
save(subsample,file="medip_esca_peak_samp146_subsamp.RData")


phen_esca$cohort=ifelse(phen_esca$Sample %in% valid_sample,"Valid_set","Train_set")
table(phen_esca$cohort)

library(gtsummary)
colnames(phen_esca)
t1=tbl_summary(phen_esca[phen_esca$cohort %in% "Train_set",c(2:8,11:12)],by = Group,
               statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
               digits = all_continuous() ~ 2,#确定小数点数
               missing_text = "(Missing)")
t2=tbl_summary(phen_esca[phen_esca$cohort %in% "Valid_set",c(2:8,11:12)],by = Group,
               statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
               digits = all_continuous() ~ 2,#确定小数点数
               missing_text = "(Missing)")
tbl_merge(tbls = list(t1, t2),tab_spanner = c("**Train_set**", "**Valid_set**"))


