## =====================================================================
## ARDS trajectory & omics analysis  
## =====================================================================

setwd('C:/Users/lxqji/OneDrive/R/25SXFX/25ARDS')
setwd('d:/OneDrive/R/25SXFX/25ARDS')

## CMAISE_ards_day1
data <- fread("./data/CMAISE_ards_day1x1.csv")
data <- fread("./data/CMAISE_ards.csv")

library(dplyr)
library(data.table)
library(DESeq2)
library(BiocManager)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(FactoMineR)
library(factoextra)
library(EnhancedVolcano)

### CMAISE1.5v
## d:/OneDrive               pf
dat <- fread("C:/Users/lxqji/OneDrive/R/01sepsis/CMAISE1.5v_x0.csv")

dat <- dat[dat$Days == 1, ]

dat <- dat %>%
  mutate(ARDS = ifelse(pf > 300, 0, 1)) %>%
  filter(copd == 0) %>%
  filter(Hospital_days >= 1) %>%
  filter(age > 17) %>%
  filter(ARDS == 1)

fwrite(dat, "./data/CMAISE_ards_day1.csv")
data <- fread("./data/CMAISE_ards_day1.csv")

ARDS_id <- data %>% select(PtID)

data <- fread("C:/Users/lxqji/OneDrive/R/01sepsis/CMAISE1.5v_x0.csv") %>%
  inner_join(ARDS_id, by = "PtID")

# Convert PtID to factor, then to numeric
data$ID <- as.numeric(as.factor(data$PtID))
data$Days <- as.numeric(data$Days)

fwrite(data, "./data/CMAISE_ards.csv")

day1 <- fread("./data/CMAISE_ards.csv") %>% filter(Days == 1)
fwrite(day1, "./data/CMAISE_ards_day1x1.csv")
data <- fread("./data/CMAISE_ards_day1x1.csv")

data <- fread("./data/CMAISE_ards.csv")

##
set.seed(11)
library(lcmm)

m1 <- hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
           subject = 'ID', ng = 1, data = data)  # 'ID'

m2 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 2, data = data))

m3 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 3, data = data))

m4 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 4, data = data))

m5 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 5, data = data))

m6 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 6, data = data))

table <- summarytable(m1, m2, m3, m4, m5, m6,
                      which = c("G", "loglik", "conv", "npm", "AIC",
                                "BIC", "SABIC", "entropy", "%class"))
flextable::save_as_docx(flextable::flextable(as.data.frame(table)),
                        path = "./table/table_lcmm_full.docx")

plot(m4, which = "fit", var.time = "Days",
     break.times = 3,  ##5
     bty = "l",
     ylab = "PaO2/FiO2",
     xlab = "Days after ICU admission",
     lwd = 2.5,  # bold line (0.5-1.5 pt print requirement)
     marg = TRUE,
     shades = TRUE,
     legend = NULL,
     col = c("#9370DB", "#8B4513", "#FF8C00", "#2E8B57"),  # replace red-yellow-blue: deep green/orange/purple/brown
     lty = rep(1, 4))  # force all trajectories to straight lines

# Add SCI-grade legend
legend("topright",
       legend = c("class 1: Rapid-resolving",
                  "class 2: Delayed-resolving",  # gradually worsening
                  "class 3: Non-resolving",      # low risk
                  "class 4: Mild & Resolving"),  # improving
       col = c("#9370DB", "#8B4513", "#FF8C00", "#2E8B57"),
       lty = rep(1, 4),  # unified straight line
       lwd = 2.5,
       cex = 0.85,  # 8-10 pt font range
       bty = "n",  # no border
       inset = c(0.35, 0.04))  # fine-tune position to avoid overlap

dtclass <- m4$pprob[, 1:2]
fwrite(dtclass, "d:/OneDrive/R/25SXFX/dtclass.csv")

datcom <- merge(data, dtclass, by = 'ID')
fwrite(datcom, "d:/OneDrive/R/25SXFX/25ARDS/datacom.csv")

datcomday1 <- datcom %>% filter(Days == 1)
fwrite(datcomday1, "d:/OneDrive/R/25SXFX/25ARDS/datacomday1.csv")

### Clinical data and expression profile common sample size: 331
### Gene sample
data <- fread("./data/CMAISE_ards.csv")
combined_counts <- fread("C:/Users/lxqji/OneDrive/R/01sepsis/OMIX006457/combined_raw_counts.csv", data.table = FALSE)
# Set first column as row names
rownames(combined_counts) <- combined_counts[, 1]
combined_counts <- combined_counts[, -1]
day1_samples <- colnames(combined_counts)[!grepl("d3$|d5$", colnames(combined_counts))]
countdata_d1 <- combined_counts[, day1_samples]
# Process clinical data and obtain intersecting samples
coldata_raw <- data %>% column_to_rownames("SampleName")
common_samples <- intersect(colnames(countdata_d1), rownames(coldata_raw))

# Final clinical data needed
coldata_final <- coldata_raw[common_samples, ]
### Clinical data and expression profile common sample size: 331
id331 <- coldata_final %>% select(PtID)

data331 <- fread("./data/CMAISE_ards.csv") %>%
  inner_join(id331, by = "PtID")

fwrite(data331, "./data/CMAISE_ards331.csv")
data <- fread("./data/CMAISE_ards331.csv")

library(lcmm)
set.seed(123)
m1 <- hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
           subject = 'ID', ng = 1, data = data)  # 'ID'

m2 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 2, data = data))

m3 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 3, data = data))

m4 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 4, data = data))

m5 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 5, data = data))

m6 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(pf ~ poly(Days, degree = 2, raw = TRUE),
                      mixture = ~ poly(Days, degree = 2, raw = TRUE),
                      subject = 'ID', ng = 6, data = data))

table <- summarytable(m1, m2, m3, m4, m5, m6,
                      which = c("G", "loglik", "conv", "npm", "AIC",
                                "BIC", "SABIC", "entropy", "%class"))
flextable::save_as_docx(flextable::flextable(as.data.frame(table)),
                        path = "./table/table_lcmm331.docx")

plot(m4, which = "fit", var.time = "Days",
     break.times = 3,  ##5
     bty = "l",
     ylab = "PaO2/FiO2",
     xlab = "Days after ICU admission",
     lwd = 2.5,  # bold line (0.5-1.5 pt print requirement)
     marg = TRUE,
     shades = TRUE,
     legend = NULL,
     col = c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513"),  # replace red-yellow-blue: deep green/orange/purple/brown
     lty = rep(1, 4))  # force all trajectories to straight lines

# Add SCI-grade legend
legend("topright",
       legend = c("class 1: Delayed-resolving",
                  "class 2: Non-resolving",  # gradually worsening
                  "class 3: Rapid-resolving",  # low risk
                  "class 4: Mild & Resolving"),  # improving
       col = c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513"),
       lty = rep(1, 4),  # unified straight line
       lwd = 2.5,
       cex = 0.85,  # 8-10 pt font range
       bty = "n",  # no border
       inset = c(0.35, 0.001))  # fine-tune position to avoid overlap

dtclass <- m4$pprob[, 1:2]
fwrite(dtclass, "d:/OneDrive/R/25SXFX/dtclass.csv")

datcom <- merge(data, dtclass, by = 'ID')
fwrite(datcom, "./data/data331.csv")

data331day1 <- datcom %>% filter(Days == 1)
fwrite(data331day1, "./data/data331day1.csv")
data331day1 <- fread("./data/data331day1.csv")

library(mice)
aa <- mice(data331day1, seed = 123)
dat <- complete(aa, action = 3)
# Outliers
dat1 <- dat[, 1:21]
dat2 <- dat[, 22:76]
data <- dat2
# Outlier automated processing 731

data$class <- dat2$class

q99 <- quantile(dat2$CRRT_days, 0.95)
dat2[dat2$CRRT_days > q99, ]$CRRT_days <- q99

data$CRRT_days <- dat2$CRRT_days
dat <- cbind(dat1, data)

dat$infectionSite_SD <- ifelse(dat$infectionSite_SD == 3 | dat$infectionSite_SD == 7 | dat$infectionSite_SD == 8, 6, dat$infectionSite_SD)

summary(dat)

fwrite(dat, "./data/data331day1mice.csv")

dat <- fread("./data/data331day1mice.csv")
library(CBCgrps)
library(flextable)
library(DataExplorer)
# plot_qq(dat, sampled_rows = 1000L)
## Overall non-normal distribution
dat <- fread("./data/data331day1mice.csv") %>%
  select(class, age, sex, diabete, hyperten, myoinfarc, cardiofailure,
         infectionSite_SD, SOFA, fluidin, fluidout, urine, hrmax, hrmin, mapmax, mapmin,
         sapmax, sapmin, rrmax, rrmin, tmax, tmin, wbc, hct, plt, pha,
         paco, pao, lac, pf, cr, crp,
         mort, MV_days, CRRT_days, VASO_days, Hospital_days)
skewvar2 <- c("lac", "paco", "pao", "pha",
              "cr", "crp", "urine", "wbc"
              , "fluidin", "fluidout", "MV_days", "CRRT_days",
              "VASO_days", "Hospital_days")
tab1 <- multigrps(dat, gvar = "class", norm.rd = 1, cat.rd = 1,
                  sk.rd = 1, minfactorlevels = 10,
                  sim = TRUE,
                  workspace = 200000, skewvar = skewvar2)
tab1
colnames(tab1) <- c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "Col7")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1)),
                        path = "./table/tab1.docx")

# table outcome

# table mort
dat <- fread("./data/data331day1mice.csv") %>%
  select(class, age, sex, diabete, hyperten,
         SOFA, fluidin, fluidout, wbc, hct, plt, pha,
         lac, pf, cr, crp, mort, Hospital_days, MV_days, CRRT_days)
skewvar2 <- c("lac", "pha", "cr", "crp", "wbc", "fluidin", "fluidout")
tab1 <- twogrps(dat, gvar = "mort", norm.rd = 1, cat.rd = 1, sk.rd = 1, minfactorlevels = 10, skewvar = skewvar2)

colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "./table/tab1_mort.docx")

library("survival")
library("survminer")
dd <- fread("./data/data331day1mice.csv") %>%
  mutate(Hospital_days = ifelse(Hospital_days > 30, 30, Hospital_days))
fit <- survfit(Surv(Hospital_days, mort) ~ class, dd)
ggsurvplot(fit,
           xlab = 'days',
           pval = T,
           pval.size = 3,
           pval.coord = c(0, 0.6),  # P-value
           legend.labs = c('class1', 'class2', 'class3', 'class4'),
           surv.median.line = "hv",  # show median survival
           ggtheme = theme_bw(),  # ggplot2 theme
           palette = c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513"),  # c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513")
           risk.table = TRUE,
           risk.table.col = "strata",
           ylim = c(0.5, 1))

dd$class <- ifelse(dd$class == 1, "class1", ifelse(dd$class == 2, "class2", ifelse(dd$class == 3, "class3", "class4")))
library(gtsummary)
library(flextable)
cox_mode1 <- coxph(Surv(Hospital_days, mort) ~ class, dd)
tbl_regression(cox_mode1, exponentiate = TRUE)
cox_mode2 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex, dd)
tbl_regression(cox_mode2, exponentiate = TRUE)
cox_mode3 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex + hyperten + diabete +
                     SOFA + copd + lac + pha + pf + cr + crp + wbc + hct + plt + urine, dd)
tbl_regression(cox_mode3, exponentiate = TRUE)

cox_mode3 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex + hyperten + diabete +
                     copd + mv + crrt + SOFA + urine + lac + pha + pf + cr + plt, dd)
tbl_regression(cox_mode3, exponentiate = TRUE)

table1 <- tbl_regression(cox_mode1, exponentiate = TRUE)
save_as_docx(as_flex_table(table1), path = "./table/cox1.docx")
table1 <- tbl_regression(cox_mode2, exponentiate = TRUE)
save_as_docx(as_flex_table(table1), path = "./table/cox2.docx")
table1 <- tbl_regression(cox_mode3, exponentiate = TRUE)
save_as_docx(as_flex_table(table1), path = "./table/cox3.docx")

setwd('C:/Users/lxqji/OneDrive/R/25SXFX/25ARDS')
setwd('d:/OneDrive/R/25SXFX/25ARDS')

dat <- fread("./data/data331day1mice.csv")
library(CBCgrps)
library(flextable)
library(DataExplorer)
library("Ztable")
library(dplyr)
library(data.table)
library(readxl)
library(openxlsx)
library(gtsummary)
library(flextable)

# plot_qq(dat, sampled_rows = 1000L)
## Overall non-normal distribution
dat <- fread("./data/data331day1mice.csv") %>%
  select(class, age, sex, diabete, hyperten, myoinfarc, cardiofailure,
         infectionSite_SD, SOFA, fluidin, fluidout, urine, hrmax, hrmin, mapmax, mapmin,
         sapmax, sapmin, rrmax, rrmin, tmax, tmin, wbc, hct, plt, pha,
         paco, pao, lac, pf, cr, crp,
         mort, MV_days, CRRT_days, VASO_days, Hospital_days)
skewvar2 <- c("lac", "paco", "pao", "pha",
              "cr", "crp", "urine", "wbc"
              , "fluidin", "fluidout", "MV_days", "CRRT_days",
              "VASO_days", "Hospital_days")
tab1 <- multigrps(dat, gvar = "class", norm.rd = 1, cat.rd = 1,
                  sk.rd = 1, minfactorlevels = 10,
                  sim = TRUE,
                  workspace = 200000, skewvar = skewvar2)
tab1
colnames(tab1) <- c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "Col7")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1)),
                        path = "./table/tab1.docx")

# table outcome

# table mort
dat <- fread("./data/data331day1mice.csv") %>%
  select(class, age, sex, diabete, hyperten,
         SOFA, fluidin, fluidout, wbc, hct, plt, pha,
         lac, pf, cr, crp, mort, Hospital_days, MV_days, CRRT_days)
skewvar2 <- c("lac", "pha", "cr", "crp", "wbc", "fluidin", "fluidout")
tab1 <- twogrps(dat, gvar = "mort", norm.rd = 1, cat.rd = 1, sk.rd = 1, minfactorlevels = 10, skewvar = skewvar2)

colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "./table/tab1_mort.docx")

library("survival")
library("survminer")
dd <- fread("./data/data331day1mice.csv") %>%
  mutate(Hospital_days = ifelse(Hospital_days > 30, 30, Hospital_days))
fit <- survfit(Surv(Hospital_days, mort) ~ class, dd)
ggsurvplot(fit,
           xlab = 'days',
           pval = T,
           pval.size = 3,
           pval.coord = c(0, 0.6),  # P-value
           legend.labs = c('class1', 'class2', 'class3', 'class4'),
           surv.median.line = "hv",  # show median survival
           ggtheme = theme_bw(),  # ggplot2 theme
           palette = c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513"),  # c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513")
           risk.table = TRUE,
           risk.table.col = "strata",
           ylim = c(0.5, 1))

### Compare 2 and 3. Class 2: 99 patients, Class 3: 28 patients
dd <- fread("./data/data331day1mice.csv")

dd$class <- ifelse(dd$class == 1, "class1", ifelse(dd$class == 2, "class2", ifelse(dd$class == 3, "class3", "class4")))

dd$class <- ifelse(dd$class == 1, "b", ifelse(dd$class == 2, "a",
                                              ifelse(dd$class == 3, "c", "d")))

# logit_models
model2 <- c("age", "sex")
model3 <- c("age", "sex", "hyperten", "diabete", "SOFA", "lac", "mv")
# categorical variables
cat_var <- c("sex", "class", "hyperten", "diabete", "mv")

tab <- logit_models(df = dd,
                    var_y = "mort",
                    var_x = c("class"),
                    adj_var = list(model1 = NULL,
                                   model2 = model2,
                                   model3 = model3),
                    cat_var = cat_var,  # as independent parameter
                    stadec = "default",
                    Pdec = "default",
                    adj_display = F,
                    path = "./table/mult_log_full.DOCX")

# cox_models
model2 <- c("age", "sex")
model3 <- c("age", "sex", "hyperten", "diabete", "SOFA", "lac", "mv")
# categorical variables
cat_var <- c("sex", "class", "hyperten", "diabete", "mv")

tab <- cox_models(df = dd,
                  var_y = "mort",
                  var_time = "Hospital_days",
                  var_x = c("class"),
                  adj_var = list(
                    model1 = NULL,
                    model2 = model2,
                    model3 = model3),
                  cat_var = cat_var,
                  stadec = "default", Pdec = "default", adj_display = F,  # show adjusted covariates T
                  path = "./table/mult_cox.DOCX")

cox_mode1 <- coxph(Surv(Hospital_days, mort) ~ class, dd)
tbl_regression(cox_mode1, exponentiate = TRUE)
cox_mode2 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex, dd)
tbl_regression(cox_mode2, exponentiate = TRUE)
cox_mode3 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex + hyperten + diabete +
                     SOFA + copd + lac + pha + pf + cr + crp + wbc + hct + plt + urine, dd)
tbl_regression(cox_mode3, exponentiate = TRUE)

cox_mode3 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex + hyperten + diabete +
                     copd + mv + crrt + SOFA + urine + lac + pha + pf + cr + plt, dd)
tbl_regression(cox_mode3, exponentiate = TRUE)

table1 <- tbl_regression(cox_mode1, exponentiate = TRUE)
save_as_docx(as_flex_table(table1), path = "./table/cox1.docx")
table1 <- tbl_regression(cox_mode2, exponentiate = TRUE)
save_as_docx(as_flex_table(table1), path = "./table/cox2.docx")
table1 <- tbl_regression(cox_mode3, exponentiate = TRUE)
save_as_docx(as_flex_table(table1), path = "./table/cox3.docx")

setwd('C:/Users/lxqji/OneDrive/R/25SXFX/25ARDS')
setwd('d:/OneDrive/R/25SXFX/25ARDS')

# All 772
data <- fread("./data/datacomday1mice.csv")
library(CBCgrps)
library(flextable)
library(DataExplorer)
library("Ztable")
library(dplyr)
library(data.table)
library(readxl)
library(openxlsx)

# plot_qq(dat, sampled_rows = 1000L)
## Overall non-normal distribution
dat <- fread("./data/datacomday1mice.csv") %>%
  select(class, age, sex, diabete, hyperten, myoinfarc, cardiofailure,
         infectionSite_SD, SOFA, fluidin, fluidout, urine, hrmax, hrmin, mapmax, mapmin,
         sapmax, sapmin, rrmax, rrmin, tmax, tmin, wbc, hct, plt, pha,
         paco, pao, lac, pf, cr, crp,
         mort, MV_days, CRRT_days, VASO_days, Hospital_days)
skewvar2 <- c("lac", "paco", "pao", "pha",
              "cr", "crp", "urine", "wbc"
              , "fluidin", "fluidout", "MV_days", "CRRT_days",
              "VASO_days", "Hospital_days")
tab1 <- multigrps(dat, gvar = "class", norm.rd = 1, cat.rd = 1,
                  sk.rd = 1, minfactorlevels = 10,
                  sim = TRUE,
                  workspace = 200000, skewvar = skewvar2)
tab1
colnames(tab1) <- c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "Col7")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1)),
                        path = "./table/tab1_full.docx")

# table outcome

# table mort
dat <- fread("./data/datacomday1mice.csv") %>%
  select(class, age, sex, diabete, hyperten,
         SOFA, fluidin, fluidout, wbc, hct, plt, pha,
         lac, pf, cr, crp, mort, Hospital_days, MV_days, CRRT_days)
skewvar2 <- c("lac", "pha", "cr", "crp", "wbc", "fluidin", "fluidout")
tab1 <- twogrps(dat, gvar = "mort", norm.rd = 1, cat.rd = 1, sk.rd = 1, minfactorlevels = 10, skewvar = skewvar2)

colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "./table/tab1_mort_full.docx")

library("survival")
library("survminer")
dd <- fread("./data/datacomday1mice.csv") %>%
  mutate(Hospital_days = ifelse(Hospital_days > 30, 30, Hospital_days))
fit <- survfit(Surv(Hospital_days, mort) ~ class, dd)
ggsurvplot(fit,
           xlab = 'days',
           pval = T,
           pval.size = 3,
           pval.coord = c(0, 0.6),  # P-value
           legend.labs = c('class1', 'class2', 'class3', 'class4'),
           surv.median.line = "hv",  # show median survival
           ggtheme = theme_bw(),  # ggplot2 theme
           palette = c("#9370DB", "#2E8B57", "#FF8C00", "#8B4513"),  # c("#2E8B57", "#FF8C00", "#9370DB", "#8B4513")
           risk.table = TRUE,
           risk.table.col = "strata",
           ylim = c(0.5, 1))

dd <- fread("./data/datacomday1mice.csv")
dd$class <- ifelse(dd$class == 1, "b", ifelse(dd$class == 2, "c",
                                              ifelse(dd$class == 3, "a", "d")))

# logit_models
model2 <- c("age", "sex")
model3 <- c("age", "sex", "hyperten", "diabete", "SOFA", "lac", "mv", "crrt", "pha", "plt")
# categorical variables
cat_var <- c("sex", "class", "hyperten", "diabete", "mv", "crrt")

tab <- logit_models(df = dd,
                    var_y = "mort",
                    var_x = c("class"),
                    adj_var = list(model1 = NULL,
                                   model2 = model2,
                                   model3 = model3),
                    cat_var = cat_var,  # as independent parameter
                    stadec = "default",
                    Pdec = "default",
                    adj_display = F,
                    path = "./table/tab_log_full.DOCX")

library(gtsummary)
library(flextable)
cox_mode1 <- coxph(Surv(Hospital_days, mort) ~ class, dd)
tbl_regression(cox_mode1, exponentiate = TRUE)
cox_mode2 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex, dd)
tbl_regression(cox_mode2, exponentiate = TRUE)

cox_mode3 <- coxph(Surv(Hospital_days, mort) ~ class + age + sex + hyperten + diabete +
                     copd + mv + crrt + SOFA + urine + lac + pha + pf + cr + plt, dd)
tbl_regression(cox_mode3, exponentiate = TRUE)

# cox_models
model2 <- c("age", "sex")
model3 <- c("age", "sex", "hyperten", "diabete", "SOFA", "lac", "mv", "copd", "crrt", "pha", "plt")
# categorical variables
cat_var <- c("sex", "class", "hyperten", "diabete", "mv", "copd", "crrt")

tab <- cox_models(df = dd,
                  var_y = "mort",
                  var_time = "Hospital_days",
                  var_x = c("class"),
                  adj_var = list(
                    model1 = NULL,
                    model2 = model2,
                    model3 = model3),
                  cat_var = cat_var,
                  stadec = "default", Pdec = "default", adj_display = F,  # show adjusted covariates T
                  path = "./table/mult_cox_full.DOCX")