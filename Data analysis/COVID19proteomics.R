library(dplyr)
library(readxl)
library(rstatix)
library(lme4)
library(lmerTest)
library(reshape2)
library(tidyverse)
library(writexl)

# MGH dataset

mydf <- read.csv('clinical_merged.csv')


# KIT and KITLG are our proteins

df <- mydf %>%
  filter(Assay %in% c('KIT', 'KITLG')) %>%
  arrange(subject_id)

# Calculate KIT and KITLG and standardize

df$KIT <- NA
df$KIT[which(df$Assay == 'KIT')] <- df$NPX[which(df$Assay == 'KIT')]
df$KITLG <- NA
df$KITLG[which(df$Assay == 'KITLG')] <- df$NPX[which(df$Assay == 'KITLG')]

# Set labels

### Timepoint

df$Timepoint <-as.factor(df$Timepoint)

levels(df$Timepoint) <- c('D0','D3','D3','D7','D7','D7','DE','DE','DE')


### COVID

df$COVID <- factor(df$COVID, levels = 0:1, labels = c('Negative','Positive'))

### Make the patient with WHO.0 = 1 NA, since it is wrong

df$WHO.0[which(df$WHO.0 == 1)] <- NA

### Make ID factor 
df$subject_id <- as.factor(df$subject_id)


### Age 

df$Age.cat <- as.factor(df$Age.cat)

levels(df$Age.cat) <- c('20-34', '36-49', '50-64', '65-79', '80+')

df$Age_comb <- as.factor(ifelse(df$Age.cat %in% c('65-79','80+'),'65+','65-'))

### Define current WHO score and make it factor

df$WHO <- NA
df$WHO[which(df$Timepoint == 'D0')] <- df$WHO.0[which(df$Timepoint == 'D0')]
df$WHO[which(df$Timepoint == 'D3')] <- df$WHO.3[which(df$Timepoint == 'D3')]
df$WHO[which(df$Timepoint == 'D7')] <- df$WHO.7[which(df$Timepoint == 'D7')]
df$WHO[which(df$Timepoint == 'DE')] <- df$WHO.28[which(df$Timepoint == 'DE')]

df$WHO <- as.factor(df$WHO)

### abs_neut

df$abs_neut <- NA
df$abs_neut[which(df$Timepoint == 'D0')] <- df$abs_neut_0_cat[which(df$Timepoint == 'D0')]
df$abs_neut[which(df$Timepoint == 'D3')] <- df$abs_neut_3_cat[which(df$Timepoint == 'D3')]
df$abs_neut[which(df$Timepoint == 'D7')] <- df$abs_neut_7_cat[which(df$Timepoint == 'D7')]

df$abs_neut_cat <- NA
df$abs_neut_cat[which(df$abs_neut %in% c(2,3))] <- 'normal'
df$abs_neut_cat[which(df$abs_neut %in% c(4,5))] <- 'high'

### abs_lymph

df$abs_lymph <- NA
df$abs_lymph[which(df$Timepoint == 'D0')] <- df$abs_lymph_0_cat[which(df$Timepoint == 'D0')]
df$abs_lymph[which(df$Timepoint == 'D3')] <- df$abs_lymph_3_cat[which(df$Timepoint == 'D3')]
df$abs_lymph[which(df$Timepoint == 'D7')] <- df$abs_lymph_7_cat[which(df$Timepoint == 'D7')]

df$abs_lymph_cat <- NA
df$abs_lymph_cat[which(df$abs_lymph %in% c(3,4))] <- 'normal'
df$abs_lymph_cat[which(df$abs_lymph == 5)] <- 'high'


### abs_mono

df$abs_mono <- NA
df$abs_mono[which(df$Timepoint == 'D0')] <- df$abs_mono_0_cat[which(df$Timepoint == 'D0')]
df$abs_mono[which(df$Timepoint == 'D3')] <- df$abs_mono_3_cat[which(df$Timepoint == 'D3')]
df$abs_mono[which(df$Timepoint == 'D7')] <- df$abs_mono_7_cat[which(df$Timepoint == 'D7')]

df$abs_mono_cat <- NA
df$abs_mono_cat[which(df$abs_mono %in% c(2,3,4))] <- 'normal'
df$abs_mono_cat[which(df$abs_mono == 5)] <- 'high'

### creat

df$creat <- NA
df$creat[which(df$Timepoint == 'D0')] <- df$creat_0_cat[which(df$Timepoint == 'D0')]
df$creat[which(df$Timepoint == 'D3')] <- df$creat_3_cat[which(df$Timepoint == 'D3')]
df$creat[which(df$Timepoint == 'D7')] <- df$creat_7_cat[which(df$Timepoint == 'D7')]

df$creat_cat <- NA
df$creat_cat <- ifelse(df$creat == 1, 'low', 'high')

### crp

df$crp <- NA
df$crp[which(df$Timepoint == 'D0')] <- df$crp_0_cat[which(df$Timepoint == 'D0')]
df$crp[which(df$Timepoint == 'D3')] <- df$crp_3_cat[which(df$Timepoint == 'D3')]
df$crp[which(df$Timepoint == 'D7')] <- df$crp_7_cat[which(df$Timepoint == 'D7')]

df$crp_cat <- NA
df$crp_cat <- ifelse(df$crp == 1, 'normal', 'high')

### ddimer

df$ddimer <- NA
df$ddimer[which(df$Timepoint == 'D0')] <- df$ddimer_0_cat[which(df$Timepoint == 'D0')]
df$ddimer[which(df$Timepoint == 'D3')] <- df$ddimer_3_cat[which(df$Timepoint == 'D3')]
df$ddimer[which(df$Timepoint == 'D7')] <- df$ddimer_7_cat[which(df$Timepoint == 'D7')]

df$ddimer_cat <- NA
df$ddimer_cat <- ifelse(df$ddimer == 1, 'normal', 'high')

### ldh

df$ldh <- NA
df$ldh[which(df$Timepoint == 'D0')] <- df$ldh_0_cat[which(df$Timepoint == 'D0')]
df$ldh[which(df$Timepoint == 'D3')] <- df$ldh_3_cat[which(df$Timepoint == 'D3')]
df$ldh[which(df$Timepoint == 'D7')] <- df$ldh_7_cat[which(df$Timepoint == 'D7')]

df$ldh_cat <- NA
df$ldh_cat[which(df$ldh %in% c(1,2))] <- 'normal'
df$ldh_cat[which(df$ldh %in% c(3,4,5))] <- 'high'

### Severity

df$WHO.max_Comb <- ifelse(df$WHO.max %in% c(1,2,3), 'Serious', 'NonSerious')
df$WHO.max_Comb[which(is.na(df$WHO.max))] <- NA
df$WHO.max_Comb <- as.factor(df$WHO.max_Comb)

df$WHO_Comb <- ifelse(df$WHO %in% c(1,2,3), 'Serious', 'NonSerious')
df$WHO_Comb <- as.factor(df$WHO_Comb)
df$WHO_Comb[which(is.na(df$WHO))] <- NA


### Keep only COVID positive

df_cov_pos <- df %>%
  filter(COVID == 'Positive')

### Remove any duplicated at each timepoint

df_cov_pos <- df_cov_pos %>%
  filter(!SampleID %in% c('321_D7.1', '321_D7.2', '344_DE.1', '344_DE.2', 	'59_D3',	'59_D3.1'))

### Comorbidities

df_cov_pos$KIDNEY <- ifelse(df_cov_pos$KIDNEY == 1, 'YES', 'NO')
df_cov_pos$HTN <- ifelse(df_cov_pos$HTN == 1, 'YES', 'NO')
df_cov_pos$DIABETES <- ifelse(df_cov_pos$DIABETES == 1, 'YES', 'NO')
df_cov_pos$LUNG <- ifelse(df_cov_pos$LUNG == 1, 'YES', 'NO')

df_cov_pos$BMI.cat <- factor(df_cov_pos$BMI.cat, levels = 0:4,
                    labels = c('Underweight', 'normal','overweight', 'obese', 'severely obese'))

df_cov_pos$HEART <- factor(df_cov_pos$HEART, levels = 0:1, labels = c('NO','YES'))
df_cov_pos$IMMUNO <- factor(df_cov_pos$IMMUNO, levels = 0:1, labels = c('NO','YES'))
df_cov_pos$Resp_Symp <- factor(df_cov_pos$Resp_Symp, levels = 0:1, labels = c('NO','YES'))
df_cov_pos$Fever_Sympt <- factor(df_cov_pos$Fever_Sympt, levels = 0:1, labels = c('NO','YES'))
df_cov_pos$GI_Symp <- factor(df_cov_pos$GI_Symp, levels = 0:1, labels = c('NO','YES'))

### Make all characters factors

df_cov_pos <- df_cov_pos %>%
  mutate_if(is.character, factor)


# Differential expression COVID vs non-COVID

kit_D0 <- df %>%
  filter(Assay == 'KIT', Timepoint == 'D0')

kit_D0 %>%
  group_by(COVID) %>%
  get_summary_stats(NPX)

levene_test(data = kit_D0, formula = NPX ~ COVID)

t.test(kit_D0$NPX ~ kit_D0$COVID, var.equal = F)

kitlg_D0 <- df %>%
  filter(Assay == 'KITLG', Timepoint == 'D0')

kitlg_D0 %>%
  group_by(COVID) %>%
  get_summary_stats(NPX)

levene_test(data = kitlg_D0, formula = NPX ~ COVID)

t.test(kitlg_D0$NPX ~ kitlg_D0$COVID, var.equal = T)

################################################################################
################################################################################

## Line-graphs

df_cov_pos %>%
  filter(Assay == 'KIT') %>%
  group_by(WHO_Comb, Timepoint) %>%
  summarise(NPX) -> kit


write.csv(kit, 'KIT.csv')

df_cov_pos %>%
  filter(Assay == 'KITLG') %>%
  group_by(WHO_Comb, Timepoint) %>%
  summarise(NPX) -> kitlg


write.csv(kitlg, 'KITLG.csv')


df_cov_pos %>%
  filter(Assay == 'KIT', Timepoint == 'D0') %>%
  select(NPX, WHO_Comb) %>%
  na.omit() -> kit0


table(kit0$WHO_Comb)

levene_test(data = kit0, formula = NPX ~ WHO_Comb)

t_test(data = kit0, formula = NPX ~ WHO_Comb, var.equal = T)

kit0 %>%
  group_by(WHO_Comb) %>%
  summarise(logFC = mean(NPX))



df_cov_pos %>%
  filter(Assay == 'KIT', Timepoint == 'D3') %>%
  select(NPX, WHO_Comb) %>%
  na.omit() -> kit3


table(kit3$WHO_Comb)

levene_test(data = kit3, formula = NPX ~ WHO_Comb)

t_test(data = kit3, formula = NPX ~ WHO_Comb, var.equal = T)

kit3 %>%
  group_by(WHO_Comb) %>%
  summarise(logFC = mean(NPX))



df_cov_pos %>%
  filter(Assay == 'KIT', Timepoint == 'D3') %>%
  select(NPX, WHO_Comb) %>%
  na.omit() -> kit7


table(kit7$WHO_Comb)

levene_test(data = kit7, formula = NPX ~ WHO_Comb)

t_test(data = kit7, formula = NPX ~ WHO_Comb, var.equal = T)

kit7 %>%
  group_by(WHO_Comb) %>%
  summarise(logFC = mean(NPX))


df_cov_pos %>%
  filter(Assay == 'KITLG', Timepoint == 'D0') %>%
  select(NPX, WHO_Comb) %>%
  na.omit() -> kitlg0


table(kitlg0$WHO_Comb)

levene_test(data = kitlg0, formula = NPX ~ WHO_Comb)

t_test(data = kitlg0, formula = NPX ~ WHO_Comb, var.equal = T)

kitlg0 %>%
  group_by(WHO_Comb) %>%
  summarise(logFC = mean(NPX))



df_cov_pos %>%
  filter(Assay == 'KITLG', Timepoint == 'D3') %>%
  select(NPX, WHO_Comb) %>%
  na.omit() -> kitlg3


table(kitlg3$WHO_Comb)

levene_test(data = kitlg3, formula = NPX ~ WHO_Comb)

t_test(data = kitlg3, formula = NPX ~ WHO_Comb, var.equal = T)

kitlg3 %>%
  group_by(WHO_Comb) %>%
  summarise(logFC = mean(NPX))



df_cov_pos %>%
  filter(Assay == 'KITLG', Timepoint == 'D7') %>%
  select(NPX, WHO_Comb) %>%
  na.omit() -> kitlg7


table(kitlg7$WHO_Comb)

levene_test(data = kitlg7, formula = NPX ~ WHO_Comb)

t_test(data = kitlg7, formula = NPX ~ WHO_Comb, var.equal = T)

kitlg7 %>%
  group_by(WHO_Comb) %>%
  summarise(logFC = mean(NPX))

################################################################################
################################################################################

# Linear mixed models - Repeated measures (without day 28)

df_cp_subE <- df_cov_pos %>%
  filter(Timepoint != 'DE') %>%
  droplevels()

kit_ancova <- lmer(KIT ~ Timepoint * WHO_Comb + (1|subject_id), data = df_cp_subE)
kitlg_ancova <- lmer(KITLG ~ Timepoint * WHO_Comb + (1|subject_id), data = df_cp_subE)

summary(kitlg_ancova)
summary(kit_ancova)

## linear mixed models for biochemical markers and comorbidities

lmer(KIT ~ abs_neut_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ abs_lymph_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ abs_mono_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ creat_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ crp_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ ddimer_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ ldh_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()


lmer(KITLG ~ abs_neut_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ abs_lymph_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ abs_mono_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ creat_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ crp_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ ddimer_cat*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ ldh_cat * Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ HEART*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ LUNG*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ KIDNEY*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ DIABETES*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ HTN*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ IMMUNO*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ Resp_Symp*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ Fever_Sympt*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ GI_Symp*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()


lmer(KITLG ~ HEART*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ LUNG*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ KIDNEY*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ DIABETES*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ HTN*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ IMMUNO*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ Resp_Symp*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ Fever_Sympt*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ GI_Symp*Timepoint + (1|subject_id), data = df_cp_subE) %>%
  summary()

# Confounding with biochemical markers, comorbidities and age
lmer(KIT ~ WHO.max_Comb + Age_comb + abs_neut_cat + abs_mono_cat + crp_cat + GI_Symp +
       (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KITLG ~ WHO.max_Comb + Age_comb + abs_neut_cat + abs_mono_cat + ldh_cat + creat_cat +
       HEART + KIDNEY + LUNG + HTN + (1|subject_id), data = df_cp_subE) %>%
  summary()

lmer(KIT ~ WHO.max_Comb + Age_comb + abs_neut_cat + abs_mono_cat + crp_cat + GI_Symp +
       (1|subject_id), data = df_cp_subE) %>%
  step() %>%
  print()

lmer(KIT ~ WHO.max_Comb + abs_mono_cat + crp_cat + GI_Symp +
       (1|subject_id), data = df_cp_subE) %>%
  summary()


d <- df_cp_subE %>%
  dplyr::select(KITLG, WHO.max_Comb, Age_comb, abs_neut_cat, abs_mono_cat, creat_cat, ldh_cat,
                HEART, KIDNEY, LUNG, HTN, subject_id) %>%
  na.omit()

lmer(KITLG ~ WHO.max_Comb + Age_comb + abs_neut_cat + abs_mono_cat + creat_cat + ldh_cat +
       HEART + KIDNEY + LUNG + HTN + (1|subject_id), data = d) %>%
  step() %>%
  print()

lmer(KITLG ~ WHO.max_Comb + abs_neut_cat + abs_mono_cat + creat_cat +
       KIDNEY  + HTN + (1|subject_id), data = d) %>%
  summary()

# Logistic Regression

## Keep only day 0
df0 <- df_cp_subE %>%
  filter(Timepoint == 'D0') %>%
  droplevels()
## Standardize KIT and KITLG of day 0

df0$KIT <- (df0$KIT - mean(df0$KIT, na.rm = T))/sd(df0$KIT, na.rm = T)
df0$KITLG <- (df0$KITLG - mean(df0$KITLG, na.rm = T))/sd(df0$KITLG, na.rm = T)


kitlg0 <- df0 %>%
  dplyr::select(KITLG, WHO.max_Comb, Age_comb) %>%
  na.omit()

kit0 <- df0 %>%
  dplyr::select(KIT, WHO.max_Comb, Age_comb) %>%
  na.omit()

kit_logist <- glm(WHO.max_Comb ~ KIT + Age_comb, family = 'binomial', data = kit0)
kitlg_logist <- glm(WHO.max_Comb ~ KITLG + Age_comb, family = 'binomial', data = kitlg0)

summary(kitlg_logist)
summary(kit_logist)

exp(coef(kit_logist))
exp(coef(kitlg_logist))

# Data for ROC curves
m_kit = kit_logist$fitted.values
d_kit = kit0$WHO.max_Comb

actuals_kit <- d_kit == 'Serious'

actuals_kit <- actuals_kit[order(m_kit)]

thresholds_kit <- sort(m_kit)

pos_kit <- sum(actuals_kit)
neg_kit <- sum(!actuals_kit)

tn_kit <- cumsum(!actuals_kit)
spec_kit <- tn_kit/neg_kit

tp_kit <- pos_kit - cumsum(actuals_kit)
sens_kit <- tp_kit/pos_kit

spec_ant_kit <- 1 - spec_kit
plot(spec_ant_kit, sens_kit, type = "l", col = "red", 
     ylab = "Sensitivity", xlab = "1 - Specificity")
abline(c(0,0),c(1,1))

mydata <- data.frame(cbind(spec_ant_kit, sens_kit))
library(writexl)
write_xlsx(mydata, 'kit_roc.xlsx')

m_kitlg = kitlg_logist$fitted.values
d_kitlg = kitlg0$WHO.max_Comb

actuals_kitlg <- d_kitlg == 'Serious'

actuals_kitlg <- actuals_kitlg[order(m_kitlg)]

thresholds_kitlg <- sort(m_kitlg)

pos_kitlg <- sum(actuals_kitlg)
neg_kitlg <- sum(!actuals_kitlg)

tn_kitlg <- cumsum(!actuals_kitlg)
spec_kitlg <- tn_kitlg/neg_kitlg

tp_kitlg <- pos_kitlg - cumsum(actuals_kitlg)
sens_kitlg <- tp_kitlg/pos_kitlg

spec_ant_kitlg <- 1 - spec_kitlg

plot(spec_ant_kitlg, sens_kitlg, type = "l", col = "red", 
     ylab = "Sensitivity", xlab = "1 - Specificity")
abline(c(0,0),c(1,1))

mydata <- data.frame(cbind(spec_ant_kitlg, sens_kitlg))
library(writexl)
write_xlsx(mydata, 'kitlg_roc.xlsx')


# Linear mixed models for volcano plots

mydf$Age <- mydf$Age.cat

mydf$Age.cat <- as.factor(mydf$Age.cat)

levels(mydf$Age.cat) <- c('20-34', '36-49', '50-64', '65-79', '80+')

mydf$Age_comb <- as.factor(ifelse(mydf$Age.cat %in% c('65-79','80+'),'65+','65-'))
levels(as.factor(mydf$Timepoint))

mydf <- mydf %>%
  filter(!Timepoint %in% c('DE', 'DE.1', 'DE.2'))

assays <- unique(mydf$Assay)

n <- length(assays)

## COVID vs non-COVID
myp <- 0
logFC <- 0
myf <- 0
mycon <- 0



for (i in 1:n) {
  mydf %>%
    filter(Assay == assays[i]) -> data
  
  data$NPX <-(data$NPX - mean(data$NPX, na.rm = T))/sd(data$NPX, na.rm = T)
  
  model <- lmer(NPX ~ COVID + (1|subject_id),
                data = data)
  
  logFC <- c(logFC, summary(model)$coefficients[2,1])
  
  myp <- c(myp,summary(model)$coefficients[2,5])
  
  if (any( grepl("failed to converge", model@optinfo$conv$lme4$messages))){
    myf <- c(myf, Assay_all[i])
  }
  
  if (any( grepl("nearly unidentifiable", model@optinfo$conv$lme4$messages))){
    mycon <- c(mycon, Assay_all[i])
  }
}

results <- data.frame(assays, logFC[-1],myp[-1])

results <- results %>%
  arrange(myp..1.)

colnames(results) <- c('Assay', 'logFC', 'p')

results$p_ad <- p.adjust(results$p, method = 'BH')

results$logp <- -log10(results$p_ad)

write.csv(results, 'LMM_COVID.txt')

## Severity

mydf$WHO.0[which(mydf$WHO.0 == 1)] <- NA

mydf$WHO <- NA
mydf$WHO[which(mydf$Timepoint == 'D0')] <- mydf$WHO.0[which(mydf$Timepoint == 'D0')]
mydf$WHO[which(mydf$Timepoint == 'D3')] <- mydf$WHO.3[which(mydf$Timepoint == 'D3')]
mydf$WHO[which(mydf$Timepoint == 'D7')] <- mydf$WHO.7[which(mydf$Timepoint == 'D7')]
mydf$WHO[which(mydf$Timepoint == 'DE')] <- mydf$WHO.28[which(mydf$Timepoint == 'DE')]

mydf$WHO <- as.factor(mydf$WHO)

mydf$WHO_Comb <- ifelse(mydf$WHO %in% c(1,2,3), 'Serious', 'NonSerious')
mydf$WHO_Comb[which(is.na(mydf$WHO))] <- NA
mydf$WHO_Comb <- as.factor(mydf$WHO_Comb)

myp <- 0
logFC <- 0
myf <- 0
mycon <- 0

for (i in 1:n) {
  mydf %>%
    filter(Assay == assays[i]) -> data
  
  data$NPX <- (data$NPX - mean(data$NPX, na.rm = T))/sd(data$NPX, na.rm = T)
  
  model <- lmer(NPX ~ WHO_Comb + (1|subject_id),
                data = data)
  
  logFC <- c(logFC, summary(model)$coefficients[2,1])
  
  myp <- c(myp,summary(model)$coefficients[2,5])
  
  if (any( grepl("failed to converge", model@optinfo$conv$lme4$messages))){
    myf <- c(myf, Assay_all[i])
  }
  
  if (any( grepl("nearly unidentifiable", model@optinfo$conv$lme4$messages))){
    mycon <- c(mycon, Assay_all[i])
  }
}


results <- data.frame(assays, logFC[-1],myp[-1])

results <- results %>%
  arrange(myp..1.)

colnames(results) <- c('Assay', 'logFC', 'p')

results$p_ad <- p.adjust(results$p, method = 'BH')

results$logp <- -log10(results$p_ad)

write.csv(results, 'LMM_Severity.txt')


### Correlations for heatmaps


wide_data <- read_excel('lmm_MGH.xlsx') %>%
  as.data.frame()

wide_data <- wide_data[,-1]

covid_data <- wide_data %>%
  filter(Timepoint == 'D0', COVID == 1)

noncovid_data <- wide_data %>%
  filter(Timepoint == 'D0', COVID == 0)


assays <- colnames(wide_data)[-c(1,1422,1423,1424)]


n <- length(assays)

correlations_noncov <- matrix(1, n, n)
pvalue_noncov <- matrix('', n, n)

correlations_cov <- matrix(1, n, n)
pvalue_cov <- matrix('', n, n)

for (i in 1:n) {
  for (j in (i+1):n){
    a <- cor.test(noncovid_data[, assays[i]],
                  noncovid_data[, assays[j]])
    
    correlations_noncov[i,j] <- a$estimate
    pvalue_noncov[i,j] <- a$p.value
    
    correlations_noncov[j,i] <- a$estimate
    pvalue_noncov[j,i] <- a$p.value
    
    b <- cor.test(covid_data[, assays[i]],
                  covid_data[, assays[j]])
    
    correlations_cov[i,j] <- b$estimate
    pvalue_cov[i,j] <- b$p.value
    
    correlations_cov[j,i] <- b$estimate
    pvalue_cov[j,i] <- b$p.value

  }

}


correlations_noncov <- as.data.frame(correlations_noncov)
pvalue_noncov <- as.data.frame(pvalue_noncov)

colnames(correlations_noncov) <- assays
row.names(correlations_noncov) <- assays

colnames(pvalue_noncov) <- assays
row.names(pvalue_noncov) <- assays



correlations_cov <- as.data.frame(correlations_cov)
pvalue_cov <- as.data.frame(pvalue_cov)

colnames(correlations_cov) <- assays
row.names(correlations_cov) <- assays

colnames(pvalue_cov) <- assays
row.names(pvalue_cov) <- assays

write.csv(correlations_noncov, 'Non_COVID_correlations.csv')
write.csv(pvalue_noncov, 'Non_COVID_pvalue.csv')

write.csv(correlations_cov, 'COVID_correlations.csv')
write.csv(pvalue_cov, 'COVID_pvalue.csv')

write_xlsx(correlations_noncov, 'Non_COVID_correlations.xlsx')
write_xlsx(pvalue_noncov, 'Non_COVID_pvalue.xlsx')

write_xlsx(correlations_cov, 'COVID_correlations.xlsx')
write_xlsx(pvalue_cov, 'COVID_pvalue.xlsx')

################################################################################

# ICL dataset



med2 <- read.csv('media-2.txt')
med3 <- read.csv('media-3.txt')

med3 <- med3 %>%
  arrange(SampleID)

med2 <- med2 %>%
  arrange(SampleID)

# Joint the two datasets
med23 <- left_join(med2,med3) %>%
  arrange(Individual_ID)

# Keep only COVID+ cases
med23_pos <- med23 %>%
  filter(Case_Control == 'POSITIVE')


# Convert to data frame
med23_pos <- as.data.frame(med23_pos)

# Classify: Critical as "serious", severe, moderate, mild as "non-serious"
med23_pos$WHO.max_Comb <- as.factor(
  ifelse(med23_pos$WHO_Severity_Peak == 'critical',
         'Serious',
         'NonSerious')
)

med23_pos$WHO_Comb <- as.factor(
  ifelse(med23_pos$WHO_Severity_Contemporaneous == 'critical',
         'Serious',
         'NonSerious')
)

# Timepoints

med23_pos$time <- NA

med23_pos$time[which(med23_pos$Time_From_First_Symptoms < 8)] <- 'D0-8'
med23_pos$time[which(med23_pos$Time_From_First_Symptoms >= 8 &
                 med23_pos$Time_From_First_Symptoms < 12)] <- 'D8-12'

med23_pos$time[which(med23_pos$Time_From_First_Symptoms >= 12 &
                 med23_pos$Time_From_First_Symptoms < 19)] <- 'D12-19'

med23_pos$time[which(med23_pos$Time_From_First_Symptoms >= 19)] <- 'D19+'


# Create the datasets for KITLG and for KIT

kitlg <- med23_pos %>%
  filter(Assay == 'SCF')

kitlg %>%
  arrange(Time_From_First_Symptoms) -> kitlg

## Keep only the non-duplicated

## Define age categories
kitlg$Age_cat <- ifelse(kitlg$Age %in% c('(60,80]', '(80,100]'), '60+','60-')



kit <- med23_pos %>%
  filter(Assay == 'KIT')

kit %>%
  arrange(Time_From_First_Symptoms) -> kit



## Define age categories
kit$Age_cat <- ifelse(kit$Age %in% c('(60,80]', '(80,100]'), '60+','60-')



## Linear mixed models for differential expression

Assay_all <- unique(med23_pos$GeneID)

n <- length(Assay_all)

myp <- 0
logFC <- 0
myf <- 0
mycon <- 0


for (i in 1:n) {
  med23 %>%
    filter(GeneID == Assay_all[i], Case_Control == 'POSITIVE') -> data
  
  data$WHO.max_Comb <- as.factor(
    ifelse(data$WHO_Severity_Peak == 'critical',
           'Serious',
           'NonSerious')
  )
  
  model <- lmer(NPX ~ WHO.max_Comb + (1|Individual_ID),
                data = data)
  
  logFC <- c(logFC, summary(model)$coefficients[2,1])
  
  myp <- c(myp,summary(model)$coefficients[2,5])
  
  if (any( grepl("failed to converge", model@optinfo$conv$lme4$messages))){
    myf <- c(myf, Assay_all[i])
  }
  
  if (any( grepl("nearly unidentifiable", model@optinfo$conv$lme4$messages))){
    mycon <- c(mycon, Assay_all[i])
  }
}

ICL_res <- data.frame(Assay_all, logFC[-1],myp[-1])

ICL_res_final  <- ICL_res %>%
  filter(!(Assay_all %in% myf))

ICL_res_final <- ICL_res_final %>%
  arrange(myp..1.)

ICL_res_final$p_ad <- p.adjust(ICL_res_final$myp..1., method = 'BH')

ICL_res_final <- ICL_res_final[,-3]

low <- ICL_res_final[which(ICL_res_final$logFC<0),]

ICL_res_final$logp <- -log10(ICL_res_final$p_ad)

write.csv(ICL_res_final, 'ICL_LMM.csv')


lmer(NPX ~ WHO.max_Comb * time + Age_cat +
       (1|Individual_ID), data = kit) %>%
  summary()

lmer(NPX ~ WHO.max_Comb * time + Age_cat +
       (1|Individual_ID), data = kitlg) %>%
  summary() -> model

model

model$coefficients
-log10(model$coefficients[,5])


### Linegraphs

med23_pos %>%
  dplyr::select(Individual_ID, Assay, time, NPX) -> data_for_lineplot

# Remove patients with more than one counts inside the same time interval
a <- paste(data_for_lineplot$Individual_ID,
           data_for_lineplot$Assay,
           data_for_lineplot$time)

data_for_lineplot <- data_for_lineplot[-which(duplicated(a)),]

spread(data_for_lineplot, key = 'time', value = 'NPX') -> wide_data_lineplot

write.csv(wide_data_lineplot, 'ICL_Data_Lineplot.csv')

## ICL correlations

## Keep only the first measurement
mydat <- med23 %>%
  arrange(Time_From_First_Symptoms)

dat1 <- mydat %>%
  dplyr::select(SampleID, UniProt, NPX)

df_wide <- dat1 %>%
  pivot_wider(names_from = UniProt, values_from = NPX)

mymerge <- med23 %>%
  dplyr::select(SampleID, Individual_ID, Case_Control, WHO_Severity_Peak, Time_From_First_Symptoms)

df_merged <- merge(df_wide, mymerge, by = 'SampleID')

df_final <- df_merged[which(!duplicated(df_merged$SampleID)),]

df_final <- df_final %>%
  arrange(Time_From_First_Symptoms)

df_final <- df_final[-which(duplicated(df_final$Individual_ID)),]

uniprots <- colnames(df_final)[2:437]

n <- length(uniprots)

correlations <- matrix(1, n, n)
pvalue <- matrix('', n, n)

df_final <- as.data.frame(df_final)

df_final <- df_final %>%
  mutate_if(is.character, factor)

for (i in 1:n) {
  for (j in (i+1):n){
    a <- cor.test(df_final[, uniprots[i]], df_final[, uniprots[j]])
    
    correlations[i,j] <- a$estimate
    pvalue[i,j] <- a$p.value
    correlations[j,i] <- a$estimate
    pvalue[j,i] <- a$p.value
  }

}


correlations <- as.data.frame(correlations)
pvalue <- as.data.frame(pvalue)

colnames(correlations) <- uniprots
row.names(correlations) <- uniprots

colnames(pvalue) <- uniprots
row.names(pvalue) <- uniprots

write.csv(correlations, 'ICL_correlations_First_Measure.txt')
write.csv(pvalue, 'ICL_pvalue_First_Measure.txt')



################################################################################

# KCL

##  Correlations

df <- read_excel('COVID19bloodbiomarkerprofilewithphenotypeModified.xlsx', sheet = 'Data')

df <- df %>%
  filter(!is.na(`Pat Code`))

df <- df %>%
  mutate_at('Pat Code', factor)

df <- df %>%
  rename(ID = 'Pat Code')


df$Group <- as.factor(df$Group)
levels(df$Group) <- c('Control', 'Mild', 'Severe', 'Critical')
df$WHO_max <- as.factor(ifelse(df$WHO == 'Severe','Serious','Non-Serious'))


## Keep only the first measurement

df <- df %>%
  arrange(`Days since onset`)

df <- df[which(!duplicated(df$ID)),]

df <- as.data.frame(df)


assays <- colnames(df)[8:362]

n <- length(assays)
logfc <- 0
p <- 0
for (i in 8:362) {
  model <- lm(df[,i] ~ WHO_max, data = df)
  logfc <- c(logfc, summary(model)$coefficients[2,1])
  p <- c(p, summary(model)$coefficients[2,4])
  
}

logfc <- logfc[-1]
p <- p[-1]

pval <- p.adjust(p, method = 'BH')

mypval <- -log10(pval)

results <- data.frame(cbind(assays, logfc, mypval))

write.csv(results, 'KCL_Critical.csv')

correlations <- matrix(1, n, n)
pvalue <- matrix('', n, n)

df <- as.data.frame(df)

df <- df %>%
  mutate_if(is.character, factor)


for (i in 1:n) {
  for (j in (i+1):n){
    a <- cor.test(df[, assays[i]],df[, assays[j]])
    
    correlations[i,j] <- a$estimate
    pvalue[i,j] <- a$p.value
    correlations[j,i] <- a$estimate
    pvalue[j,i] <- a$p.value
  }
}

correlations <- as.data.frame(correlations)
pvalue <- as.data.frame(pvalue)

colnames(correlations) <- assays
row.names(correlations) <- assays

colnames(pvalue) <- assays
row.names(pvalue) <- assays

write.csv(correlations, 'KCL_correlations.txt')
write.csv(pvalue, 'KCL_pvalue.txt')
