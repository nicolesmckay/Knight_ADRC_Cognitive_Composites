#### Knight ADRC PACC and Domain Specific Composites Project 
#### Authors: Nicole S. McKay, Peter R. Millar, Andrew Aschenbrenner, Jessica Nicosia, Jason Hassenstab
#### Finalized: November 2023
#### Last Changes: January 2024 #added extra comments 

#### This code will create domain-specific composites, and a Knight-PACC, and compares these to the ADCS-PACC and a Global cognitive composite 

#### Please see manuscript "Comparing domain-specific and general cognitive composites in Alzheimer disease research" Accepted 2024 at Neuropsychology 

#### Packages ####
# If you do not have pacman already on your computer, you must first run: install.packages("pacman") 
pacman::p_load(data.table, corrplot, psych, jtools, ggpubr, effectsize, ppcor, tidyverse, gdata, RVAideMemoire)
scaleFUN <- function(x) sprintf("%.1f", x)

#### Load Data and Set Positivity Variables ####
# For these lines of code to work, you need to adjust the file location and name to match where you have saved the data files. Note, you may also need to convert these to .csv
# If there are issues, please also check variable names match the names in your spreadsheet
demo <- read.csv("Data/mod_demographics.csv")%>% 
  mutate(birth_Date = as.Date(BIRTH, "%m/%d/%Y"))%>%
  distinct(ID, birth_Date, EDUC, race, sex)%>% mutate(sex = as.numeric(str_replace_all(sex, c("M" = "1", "F" = "2"))))

apoe <- read.csv("Data/apoe.csv")%>% 
  mutate(ID = id)%>% distinct(ID, apoe)%>% mutate(apoe = as.numeric(str_replace_all(apoe, c("22" = "0", "24" = "1", "44" = "1", "34" = "1", "23" = "0", "33" = "0"))))

clin <- read.csv("Data/mod_b4_cdr.csv")%>% 
  mutate(clin_Date = as.Date(TESTDATE, "%m/%d/%Y"), cdr_bin = case_when(cdr == 0 ~ 0, cdr > 0 ~ 1), cdr_bin = as.numeric(cdr_bin))%>% 
  distinct(ID, clin_Date, cdr_bin, cdr, MMSE)

cog <- read.csv("Data/mod_psychometrics.csv")%>% 
  select(ID, psy_date, mocatots, craftdre, MEMUNITS, srtfree, asscmem, minttots, ANIMALS, VEG, BOSTON, tma, tmb, trailb, TRAILA, digsym, switchmixed, DIGIF, DIGIB, digforct, digbacct, lettnum)%>% 
  mutate(psy_Date = as.Date(psy_date, "%m/%d/%Y"))%>% 
  select (-psy_date)%>% 
  mutate(MEMUNITS = case_when(MEMUNITS <= 25 ~ MEMUNITS),switchmixed = case_when(switchmixed <= 52 ~ switchmixed))

tau_pos <- read.csv("Data/tau_fs53.csv")%>% 
  mutate(tau_Date = as.Date(PET_Date, "%m/%d/%Y"))%>% 
  select(ID, tau_Date, Tauopathy)%>% 
  drop_na()%>% 
  mutate(taupos = case_when (Tauopathy >= 1.22 ~ 1, Tauopathy < 1.22 ~ 0))%>% 
  distinct(ID, tau_Date, taupos)

pib <- read.csv("Data/pib_fs53.csv")%>% 
  drop_na(pib_fsuvr_rsf_tot_cortmean)%>% 
  mutate(Date = as.Date(PET_Date, "%m/%d/%Y"), pos = case_when( pib_fsuvr_rsf_tot_cortmean >= 1.42 ~ 1, pib_fsuvr_rsf_tot_cortmean < 1.42 ~ 0))%>% 
  select(ID, Date, pos)
av <- read.csv("Data/av45_fs53.csv")%>% 
  drop_na(a45_fsuvr_rsf_tot_cortmean)%>% 
  mutate(Date = as.Date(PET_Date, "%m/%d/%Y"), pos = case_when(a45_fsuvr_rsf_tot_cortmean >= 1.16 ~ 1, a45_fsuvr_rsf_tot_cortmean < 1.16 ~ 0))%>% 
  select(ID, Date, pos) 
csf_pos <- read.csv("Data/csf_markers.csv")%>% 
  mutate(Date = as.Date(CSF_LP_DATE, "%m/%d/%Y"), AB4240 = LUMIPULSE_CSF_AB42/LUMIPULSE_CSF_AB40, pos = case_when (AB4240 <= 0.0673 ~ 1, AB4240 > 0.0673 ~ 0))%>% 
  drop_na(AB4240)%>% 
  select(ID, Date, pos)

# combine the two amyloid-PET dataframes and then add the CSF dataframe to create one combined amyloid positivity dataframe
pet_pos <- rbind(pib, av)%>% 
  arrange(ID)
ab_pos <- rbind(pet_pos, csf_pos)%>% 
  arrange(ID, Date)%>% 
  group_by(ID)%>% 
  mutate(poslag = lag(pos), abpos = case_when(is.na(poslag) ~ pos, poslag == pos ~ pos, poslag!= pos ~ 1))%>% 
  distinct(ID, Date, abpos)%>% 
  ungroup()

rm(av, csf_pos, pet_pos, pib)

#### Merge the Data to have one cog data form ####
# Merging, but only keeping if biomarker is within three years of cognition 
full_data <- merge(cog, apoe, by = "ID", all.x = TRUE, sort = TRUE) #apoe
full_data <- merge(full_data, demo, by = "ID", all.x = TRUE, sort = TRUE) #demographics
full_data <- merge(full_data, clin, by = "ID", all.x = TRUE)%>% #clinical
  mutate(abs_diff=as.numeric(abs(clin_Date-psy_Date)), abs_diff = replace_na(abs_diff, 0))%>% 
  group_by(ID, psy_Date)%>%  
  slice_min(abs_diff)%>% 
  ungroup()%>% 
  filter(abs_diff <= (365.25*3))%>% 
  select(-abs_diff)

full_data <- merge(full_data, tau_pos, by = "ID", all.x = TRUE)%>% #tau_pos
  mutate(abs_diff=as.numeric(abs(tau_Date-psy_Date)), abs_diff = replace_na(abs_diff, 0))%>% 
  group_by(ID, psy_Date)%>%  
  slice_min(abs_diff)%>% 
  ungroup()%>% 
  filter(abs_diff <= (365.25*3))%>% 
  select(-abs_diff)

full_data <- merge(full_data, ab_pos, by = "ID", all.x = TRUE)%>% #ab_pos
  mutate(abs_diff=as.numeric(abs(Date-psy_Date)), abs_diff = replace_na(abs_diff, 0))%>% 
  group_by(ID, psy_Date)%>%  
  slice_min(abs_diff)%>% 
  ungroup()%>% 
  filter(abs_diff <= (365.25*3))%>% 
  select(-abs_diff)

# fill in variables missing based on other visits (will fill top to bottom), also restrict to 2005 and later and ignore anything collected during the pandemic
full_data <- full_data%>% 
  group_by(ID)%>% 
  arrange(psy_Date)%>% 
  fill(taupos)%>% 
  fill(abpos)%>% 
  fill(apoe, .direction = "downup")%>% 
  fill(birth_Date, .direction = "downup")%>% 
  fill(EDUC, .direction = "downup")%>% 
  fill(race, .direction = "downup")%>% 
  fill(sex, .direction = "downup")%>% 
  ungroup()%>% 
  select(-tau_Date, -Date, -clin_Date)%>% 
  mutate(year = format(psy_Date, format = "%Y"))%>% 
  filter(year > 2003)%>% 
  filter(year < 2020)%>% 
  mutate(Age = as.numeric(round(psy_Date - birth_Date)/365))

rm(ab_pos, apoe, clin, cog, demo, tau_pos)

#### EQUATE VARIABLES ####
# remove any data values that are invalid
full_data <- full_data %>%
  mutate(MMSE = as.numeric(MMSE), mocatots = as.numeric(mocatots), eqMMSE = case_when (
    mocatots == 30 ~ 30, mocatots == 29 ~ 30, mocatots == 28 ~ 30,
    mocatots == 27 ~ 30, mocatots == 26 ~ 29, mocatots == 25 ~ 29,
    mocatots == 24 ~ 29, mocatots == 23 ~ 28, mocatots == 22 ~ 28, 
    mocatots == 21 ~ 27, mocatots == 20 ~ 26, mocatots == 19 ~ 25,
    mocatots == 18 ~ 24, mocatots == 17 ~ 23, mocatots == 16 ~ 22,
    mocatots == 15 ~ 21, mocatots == 14 ~ 20, mocatots == 13 ~ 19,
    mocatots == 12 ~ 18, mocatots == 11 ~ 17, mocatots == 10 ~ 16,
    mocatots == 9 ~ 15, mocatots == 8 ~15, mocatots == 7 ~ 14,
    mocatots == 6 ~ 13, mocatots == 5 ~ 12, mocatots == 4 ~ 12,
    mocatots == 3 ~ 11, mocatots == 2 ~ 10, mocatots == 1 ~ 9,
    mocatots == 0 ~ 6, is.na(mocatots) ~ MMSE), 
    BOSTON = as.numeric(BOSTON), minttots = as.numeric(minttots), eqMINT = case_when (
      minttots == 32 ~ 30, minttots == 31 ~ 29, minttots == 30 ~ 28,
      minttots == 29 ~ 27, minttots == 28 ~ 26, minttots == 27 ~ 25,
      minttots == 26 ~ 24, minttots == 25 ~ 22, minttots == 24 ~ 21,
      minttots == 23 ~ 20, minttots == 22 ~ 18, minttots == 21 ~ 17,
      minttots == 20 ~ 16, minttots == 19 ~ 15, minttots == 18 ~ 14,
      minttots == 17 ~ 13, minttots == 16 ~ 12, minttots == 15 ~ 11,
      minttots == 14 ~ 11, minttots == 13 ~ 10, minttots == 12 ~ 9,
      minttots == 11 ~ 9, minttots == 10 ~ 8, minttots == 9 ~ 8,
      minttots == 8 ~ 7, minttots == 7 ~ 7, minttots == 6 ~ 6,
      minttots == 5 ~ 6, minttots == 4 ~ 5, minttots == 3 ~ 4,
      minttots == 2 ~ 3, minttots == 1 ~ 2, minttots == 0 ~ 1,
      is.na(minttots) ~ BOSTON), 
    MEMUNITS = as.numeric(MEMUNITS), craftdre = as.numeric(craftdre),eqCRFTDel = case_when (
      craftdre == 25 ~ 24, craftdre == 24 ~ 23, craftdre == 23 ~ 22,     
      craftdre == 22 ~ 21, craftdre == 21 ~ 20, craftdre == 20 ~ 19, 
      craftdre == 19 ~ 18, craftdre == 18 ~ 16, craftdre == 17 ~ 15, 
      craftdre == 16 ~ 14, craftdre == 15 ~ 13, craftdre == 14 ~ 12, 
      craftdre == 13 ~ 11, craftdre == 12 ~ 10, craftdre == 11 ~ 9, 
      craftdre == 10 ~ 8, craftdre == 9 ~ 8, craftdre == 8 ~ 7, 
      craftdre == 7 ~ 6, craftdre == 6 ~ 5, craftdre == 5 ~ 5, 
      craftdre == 4 ~ 4, craftdre == 3 ~ 3, craftdre == 2 ~ 3, 
      craftdre == 1 ~ 1, craftdre == 0 ~ 0, is.na(craftdre) ~ MEMUNITS), 
    DIGIF = as.numeric(DIGIF), digforct = as.numeric(digforct), eqDIGFCT = case_when (
      digforct == 14 ~ 12, digforct == 13 ~ 12, 
      digforct == 12 ~ 11, digforct == 11 ~ 11, digforct == 10 ~ 10, 
      digforct == 9 ~ 9, digforct == 8 ~ 9, digforct == 7 ~ 8, 
      digforct == 6 ~ 7, digforct == 5 ~ 6, digforct == 4 ~ 5, 
      digforct == 3 ~ 3, digforct == 2 ~ 2, digforct == 1 ~ 1, 
      digforct == 0 ~ 0, is.na(digforct) ~ DIGIF), DIGIB = as.numeric(DIGIB), digbacct = as.numeric(digbacct), eqDIGBCT = case_when (
        digbacct == 14 ~ 12, digbacct == 13 ~ 12, digbacct == 12 ~ 12, 
        digbacct == 11 ~ 11, digbacct == 10 ~ 10, digbacct == 9 ~ 9, 
        digbacct == 8 ~ 8, digbacct == 7 ~ 7, digbacct == 6 ~ 6, 
        digbacct == 5 ~ 5, digbacct == 4 ~ 4, digbacct == 3 ~ 3, 
        digbacct == 2 ~ 2, digbacct == 1 ~ 1, digbacct == 0 ~ 0, 
        is.na(digbacct) ~ DIGIB), tma = as.numeric(tma), tmb = as.numeric(tmb), TRAILA = as.numeric(TRAILA), trailb = as.numeric(trailb), tmax = case_when(
          tma > 149 ~ 150, tma < 150 ~ tma), TRAILAx = case_when(
            TRAILA > 149 ~ 150, TRAILA < 150 ~ TRAILA), tmbx = case_when(
              tmb > 300 ~ 300, tmb < 300 ~ tmb), TRAILBx = case_when(
                trailb > 299 ~ 300, trailb < 300 ~ trailb),TRAILA_Final = case_when (
                  is.na(TRAILAx) ~ tmax, TRAILAx > 0 ~ TRAILAx), TRAILB_Final = case_when (
                    is.na(TRAILBx) ~ tmbx, TRAILBx > 0 ~ TRAILBx))%>% 
  select(ID, Age, abpos, taupos, psy_Date, year, cdr_bin, cdr, apoe, EDUC, race, sex, birth_Date, eqMMSE, eqCRFTDel, srtfree, asscmem, VEG, ANIMALS, eqMINT, TRAILA_Final, TRAILB_Final, digsym, switchmixed, eqDIGBCT, eqDIGFCT, lettnum)

full_data$na_count <- apply(full_data [, c(13:26)], 1, function(x) sum(is.na(x)))


#### EXPLORATORY FACTOR ANALYSIS + PLOTS #### 
# Create dataset for the EFA by selecting cdr_bin = 0 and latest visit 
data_efa <- full_data%>% filter(cdr_bin == 0)%>% arrange(psy_Date)%>% group_by(ID)%>% slice_max(psy_Date)%>% ungroup()

#subset to tasks only (for EFA) and get correlation matrix 
cog <- data_efa%>% select(-ID, -abpos, -taupos, -psy_Date, -cdr_bin, -cdr, -apoe, -EDUC, -race, -sex, -birth_Date, -eqMMSE, -year, -na_count, -Age)
cog_cor <- round(cor(cog, use = "pairwise"),digits =3)

# visualize heatmap and screeplots 
corrplot(cog_cor, type = "lower" , method = "color", order = "hclust", addCoef.col = "black", number.cex=0.6, tl.cex = 0.8, cl.cex = 0.5, insig ='blank')
plot(eigen(cog_cor)$values,xlab = "Factors",ylab = "Eigen Values", type ="b")

# run EFA and make path diagram 
efa <- fa(cog_cor, nfactors = 3, rotate="varimax")
fa.diagram(efa, digits = 2)

# Store demographic data for later 
data_efa <- data_efa%>% 
  select(ID, abpos, taupos, year, cdr_bin, cdr, apoe, EDUC, race, sex, eqMMSE, Age)
write.csv(data_efa, 'Manuscript/demo_efa.csv', row.names = FALSE)

# calculate the proportion of variance explained by each factor in the scree plot
proportions_explained <- data.frame(eigen(cog_cor)$values/sum(eigen(cog_cor)$values))
colnames(proportions_explained) <- "variance"
proportions_explained <- proportions_explained %>%
  mutate(Factor = 1:13, eigen = eigen(cog_cor)$values)
write.csv(proportions_explained, 'Manuscript/proportions_explained_EFA.csv', row.names = FALSE)

# Create nicer scree plot
ggplot(aes (x = Factor, y = eigen), data = proportions_explained)+
  geom_line()+
  geom_point()+ 
  scale_y_continuous(
    "Eigenvalues", 
    sec.axis = sec_axis(~ . /13 * 100, name = "Percentage of Variance Explained"), labels = scaleFUN)+
  labs(x = "Factors")+
  scale_x_continuous(breaks = 1:13)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), axis.text = element_text(size = 10), axis.title = element_text (size = 12))
ggsave("Manuscript/final_scree_plot.png", width=10,height =4)  
#keep(full_data, scaleFUN, sure = TRUE)

#### COMPUTE DERIVED FACTORS IN FULL DATASET ####
# Create control population to create z-scores relative to
data_cdr_bin0 <- full_data%>% filter(cdr_bin == 0)

#Calculate the means and sd of cdr_bin = 0, then use these to create z-scores for each relevant column for everyone's latest (pre 2020) visit 
m_eqMMSE  <- mean(data_cdr_bin0$eqMMSE, na.rm = T); sd_eqMMSE <- sd(data_cdr_bin0$eqMMSE, na.rm = T); full_data$z_eqMMSE <- (full_data$eqMMSE - m_eqMMSE)/sd_eqMMSE
m_eqMINT <- mean(data_cdr_bin0$eqMINT, na.rm = T); sd_eqMINT <- sd(data_cdr_bin0$eqMINT, na.rm = T); full_data$z_eqMINT <- (full_data$eqMINT - m_eqMINT)/sd_eqMINT
m_eqCRFTDel <- mean(data_cdr_bin0$eqCRFTDel, na.rm = T); sd_eqCRFTDel <- sd(data_cdr_bin0$eqCRFTDel, na.rm = T); full_data$z_eqCRFTDel <- (full_data$eqCRFTDel - m_eqCRFTDel)/sd_eqCRFTDel
m_eqDIGFCT  <- mean(data_cdr_bin0$eqDIGFCT, na.rm = T); sd_eqDIGFCT <- sd(data_cdr_bin0$eqDIGFCT, na.rm = T); full_data$z_eqDIGFCT <- (full_data$eqDIGFCT - m_eqDIGFCT)/sd_eqDIGFCT
m_eqDIGBCT  <- mean(data_cdr_bin0$eqDIGBCT, na.rm = T); sd_eqDIGBCT <- sd(data_cdr_bin0$eqDIGBCT, na.rm = T); full_data$z_eqDIGBCT <- (full_data$eqDIGBCT - m_eqDIGBCT)/sd_eqDIGBCT
m_switchmixed <- mean(data_cdr_bin0$switchmixed, na.rm = T); sd_switchmixed <- sd(data_cdr_bin0$switchmixed, na.rm = T); full_data$z_switchmixed <- (full_data$switchmixed- m_switchmixed)/sd_switchmixed
m_asscmem <- mean(data_cdr_bin0$asscmem, na.rm = T); sd_asscmem <- sd(data_cdr_bin0$asscmem, na.rm = T); full_data$z_asscmem <- (full_data$asscmem - m_asscmem)/sd_asscmem
m_TRAILA_Final  <- mean(data_cdr_bin0$TRAILA_Final, na.rm = T); sd_TRAILA_Final <- sd(data_cdr_bin0$TRAILA_Final, na.rm = T); full_data$z_TRAILA_Final <- -1*(full_data$TRAILA_Final - m_TRAILA_Final)/sd_TRAILA_Final
m_TRAILB_Final <- mean(data_cdr_bin0$TRAILB_Final, na.rm = T); sd_TRAILB_Final <- sd(data_cdr_bin0$TRAILB_Final, na.rm = T); full_data$z_TRAILB_Final <- -1*(full_data$TRAILB_Final - m_TRAILB_Final)/sd_TRAILB_Final
m_lettnum <- mean(data_cdr_bin0$lettnum, na.rm = T); sd_lettnum <- sd(data_cdr_bin0$lettnum, na.rm = T); full_data$z_lettnum <- (full_data$lettnum - m_lettnum)/sd_lettnum
m_ANIMALS <- mean(data_cdr_bin0$ANIMALS, na.rm = T); sd_ANIMALS <- sd(data_cdr_bin0$ANIMALS, na.rm = T); full_data$z_ANIMALS <- (full_data$ANIMALS - m_ANIMALS)/sd_ANIMALS
m_VEG  <- mean(data_cdr_bin0$VEG, na.rm = T); sd_VEG <- sd(data_cdr_bin0$VEG, na.rm = T); full_data$z_VEG<- (full_data$VEG - m_VEG)/sd_VEG
m_srtfree  <- mean(data_cdr_bin0$srtfree, na.rm = T); sd_srtfree <- sd(data_cdr_bin0$srtfree, na.rm = T); full_data$z_srtfree <- (full_data$srtfree - m_srtfree)/sd_srtfree
m_digsym  <- mean(data_cdr_bin0$digsym, na.rm = T); sd_digsym <- sd(data_cdr_bin0$digsym, na.rm = T); full_data$z_digsym <- (full_data$digsym - m_digsym)/sd_digsym

# Count NA for each factor for each person
full_data$F1_na_count <- apply(full_data [, c("z_eqCRFTDel", "z_asscmem", "z_srtfree")], 1, function(x) sum(is.na(x)))
full_data$F2_na_count <- apply(full_data [, c("z_ANIMALS", "z_eqMINT", "z_VEG")], 1, function(x) sum(is.na(x)))
full_data$F3_na_count <- apply(full_data [, c("z_TRAILA_Final", "z_TRAILB_Final", "z_digsym", "z_switchmixed")], 1, function(x) sum(is.na(x)))
full_data$F4_na_count <- apply(full_data [, c("z_eqDIGBCT", "z_eqDIGFCT", "z_lettnum")], 1, function(x) sum(is.na(x)))
full_data$K1_na_count <- apply(full_data [, c("z_srtfree", "z_digsym", "z_ANIMALS", "z_TRAILB_Final")], 1, function(x) sum(is.na(x)))
full_data$G1_na_count <- apply(full_data [, c("z_eqCRFTDel", "z_asscmem", "z_srtfree","z_ANIMALS", "z_eqMINT", "z_VEG","z_TRAILA_Final", "z_TRAILB_Final", "z_digsym", "z_switchmixed","z_eqDIGBCT", "z_eqDIGFCT", "z_lettnum")], 1, function(x) sum(is.na(x)))
full_data$A1_na_count <- apply(full_data [, c("z_eqCRFTDel", "z_digsym", "z_eqMMSE", "z_srtfree")], 1, function(x) sum(is.na(x)))

# Using the na_counts from above, calculate the scores, if possible, if too many NAs for any specific factor, then give an NA (also divide by adjusted denominator)
index_craftdel <- grep("z_eqCRFTDel", colnames(full_data))
index_asscmem <- grep("z_asscmem", colnames(full_data))
index_srtfree <- grep("z_srtfree", colnames(full_data))
index_animals <- grep("z_ANIMALS", colnames(full_data))
index_mint <- grep("z_eqMINT", colnames(full_data))
index_veg <- grep("z_VEG", colnames(full_data))
index_tma <- grep("z_TRAILA_Final", colnames(full_data))
index_tmb <- grep("z_TRAILB_Final", colnames(full_data))
index_digsym <- grep("z_digsym", colnames(full_data))
index_switch <- grep("z_switchmixed", colnames(full_data))
index_digb <- grep("z_eqDIGBCT", colnames(full_data))
index_digf <- grep("z_eqDIGFCT", colnames(full_data))
index_lettnum <- grep("z_lettnum", colnames(full_data))
index_mmse <- grep("z_eqMMSE", colnames(full_data))
full_data$G1 <- rowMeans(full_data[, c(index_craftdel,index_asscmem,index_srtfree,
                                           index_animals,index_mint, index_veg,
                                           index_tma, index_tmb, index_digsym, index_switch,
                                           index_digb, index_digf, index_lettnum)], na.rm = TRUE)
full_data$F1 <- rowMeans(full_data[, c(index_craftdel,index_asscmem,index_srtfree)], na.rm = TRUE)
full_data$F2 <- rowMeans(full_data[, c(index_animals,index_mint,index_veg)], na.rm = TRUE)
full_data$F3 <- rowMeans(full_data[, c(index_tma,index_tmb,index_digsym, index_switch)], na.rm = TRUE)
full_data$F4 <- rowMeans(full_data[, c(index_digb,index_digf,index_lettnum)], na.rm = TRUE)
full_data$K1 <- rowMeans(full_data[, c(index_srtfree,index_animals,index_digsym, index_tmb)], na.rm = TRUE)
full_data$A1 <- rowMeans(full_data[, c(index_srtfree,index_craftdel,index_digsym, index_mmse)], na.rm = TRUE)

## if too many components were missing then set to missing/NA
full_data$F1 = ifelse(full_data$F1_na_count <2, full_data$F1, NA)
full_data$F2 = ifelse(full_data$F2_na_count <2, full_data$F2, NA)
full_data$F3 = ifelse(full_data$F3_na_count <3, full_data$F3, NA)
full_data$F4 = ifelse(full_data$F4_na_count <2, full_data$F4, NA)
full_data$K1 = ifelse(full_data$K1_na_count <3, full_data$K1, NA)
full_data$G1 = ifelse(full_data$G1_na_count <8, full_data$G1, NA)
full_data$A1 = ifelse(full_data$A1_na_count <4, full_data$A1, NA)

full_data = full_data %>% select(ID, abpos, taupos, psy_Date, cdr, cdr_bin, Age, apoe, EDUC, race, sex, birth_Date, F1, F2, F3, F4, K1, G1, A1, eqMMSE, na_count, F1_na_count, F2_na_count, F3_na_count, F4_na_count, K1_na_count, G1_na_count, A1_na_count) 

#remove outliers, anyone with values exceeding +/- 3.5 SD from mean for group
desc <- data.frame(t(describe(full_data[c("F1","F2","F3","F4","K1","G1","A1")])))

full_data <- full_data %>%
  mutate(F1 = case_when(F1 <= (desc$F1[3]+3.5*desc$F1[4]) & F1 >= (desc$F1[3]-3.5*desc$F1[4]) ~ F1))%>%
  mutate(F2 = case_when(F2 <= (desc$F2[3]+3.5*desc$F2[4]) & F2 >= (desc$F2[3]-3.5*desc$F2[4]) ~ F2))%>%
  mutate(F3 = case_when(F3 <= (desc$F3[3]+3.5*desc$F3[4]) & F3 >= (desc$F3[3]-3.5*desc$F3[4]) ~ F3))%>%
  mutate(F4 = case_when(F4 <= (desc$F4[3]+3.5*desc$F4[4]) & F4 >= (desc$F4[3]-3.5*desc$F4[4]) ~ F4))%>%
  mutate(K1 = case_when(K1 <= (desc$K1[3]+3.5*desc$K1[4]) & K1 >= (desc$K1[3]-3.5*desc$K1[4]) ~ K1))%>%
  mutate(A1 = case_when(A1 <= (desc$A1[3]+3.5*desc$A1[4]) & A1 >= (desc$A1[3]-3.5*desc$A1[4]) ~ A1))%>%
  mutate(G1 = case_when(G1 <= (desc$G1[3]+3.5*desc$G1[4]) & G1 >= (desc$G1[3]-3.5*desc$G1[4]) ~ G1))

#remove excess variables etc 
#keep(full_data, scaleFUN, sure = TRUE)

#### CREATE EXTRA VARIABLES ####
full_data <- full_data%>% 
  drop_na(cdr_bin)%>%
  group_by(ID)%>% 
  arrange(psy_Date)%>%
  mutate(min_Date = psy_Date[1], Time = as.numeric(psy_Date - min_Date), BLAge = Age[1], BLcdr_bin = cdr_bin[1], maxcdr_bin = max(cdr_bin), con = case_when(maxcdr_bin - BLcdr_bin == 0 ~ 0, maxcdr_bin - BLcdr_bin != 0 ~ 1), maxDate = max(psy_Date))%>% 
  ungroup()%>% mutate(Year = format(psy_Date, format = "%Y"))

#### FIGURE 2: HISTOGRAMS OF COGNITION BY CDR #### 
data_hist <- full_data%>% arrange(psy_Date)%>% group_by(ID)%>% slice_max(psy_Date)%>% ungroup()%>% mutate(cdr_bin = as.factor(cdr_bin))

#get mean values to place onto histograms in plotting step
means <- data_hist%>% gather (Factor, Cognition, F1:A1)%>% mutate(Factor = as.factor(Factor))%>% select(Factor, Cognition, cdr_bin)%>% group_by(Factor, cdr_bin)%>% drop_na()%>% summarise_at(vars(Cognition), list(MeanF = mean))

data_hist%>% gather (Factor, Cognition, F1:A1)%>% 
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x=Cognition, color = cdr_bin, fill = cdr_bin)) + 
  geom_histogram(alpha = 0.3, position = "identity", bins = 100)+
  scale_color_manual(values = c("#EBB261","#5A4A6F"), labels=c('Cognitively Normal', 'Cognitively Impaired'))+
  scale_fill_manual(values = c("#EBB261","#5A4A6F"), labels=c( 'Cognitively Normal', 'Cognitively Impaired'))+
  scale_y_continuous(labels = scaleFUN)+
  labs (y="Frequency", x = "Standardized Cognitive Score", color = "Cognitive Status", fill = "Cognitive Status") +
  geom_vline(aes(xintercept=MeanF, color = cdr_bin),  means, linetype = "dashed")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggsave("Manuscript/cdr_bin_hist_lastvisit.png", width=10,height =4)  

data_hist <- data_hist%>% select(ID, abpos, taupos, cdr_bin, cdr, Age, apoe, EDUC, race, sex, eqMMSE, F1, F2, F3, F4, K1, A1, G1)
write.csv(data_hist, "Manuscript/demographics_histograms.csv")

#keep(full_data, scaleFUN, sure = TRUE)
#### INTRODUCE AND ORGANIZE BIOMARKER DATA ####
# Read in MRI
mri <- read.csv("Data/mri_t3_fs53.csv")%>% mutate(mri_Date = as.Date(MR_Date, "%m/%d/%Y"), CortSig = LOAD_CorticalSignature_Thickness)%>% select(ID, mri_Date, CortSig)%>% drop_na(CortSig)

# Read in TAU
tau <- read.csv("Data/tau_fs53.csv")%>% mutate(tau_Date = as.Date(PET_Date, "%m/%d/%Y"))%>% select(ID, tau_Date, Tauopathy)%>% drop_na(Tauopathy)

# Read in PIB
pib <- read.csv("Data/pib_fs53.csv")%>% mutate(ab_Date = as.Date(PET_Date, "%m/%d/%Y"), suvr = pib_fsuvr_rsf_tot_cortmean, Centiloid = 45.0*suvr - 47.5, Tracer = "pib")%>% select(ID, ab_Date, suvr, Centiloid, Tracer)%>% drop_na(suvr)

av <- read.csv("Data/av45_fs53.csv")%>% mutate(ab_Date = as.Date(PET_Date, "%m/%d/%Y"), suvr = a45_fsuvr_rsf_tot_cortmean, Centiloid = 53.6*suvr - 43.2, Tracer = "av")%>% select(ID, ab_Date, suvr, Centiloid, Tracer)%>% drop_na(suvr)

# merge the two types of amyloid-PET
ab <- rbind(pib, av)

# Remove Outliers 
desc <- data.frame(t(describe(mri$CortSig))) 
mri <- mri%>% mutate(CortSig = case_when(CortSig >= (desc$X1[3]-3*desc$X1[4]) & CortSig <= (desc$X1[3]+3*desc$X1[4])~ CortSig))

desc <- data.frame(t(describe(tau$Tauopathy)))
tau <- tau%>% mutate(Tauopathy = case_when(Tauopathy <= (desc$X1[3]+3*desc$X1[4]) & Tauopathy >= (desc$X1[3]-3*desc$X1[4]) ~ Tauopathy))

desc <- data.frame(t(describe(ab$Centiloid)))
ab <- ab%>% mutate(Centiloid = case_when(Centiloid <= (desc$X1[3]+3*desc$X1[4]) & Centiloid >= (desc$X1[3]-3*desc$X1[4]) ~ Centiloid))
rm(desc, pib, av) 

# add closest cognition to each biomarker and then subset to the most complete first visit based on cognition values missing (and time must be within three years of biomarker)
ab <- merge(ab, full_data, by = "ID", all.x = TRUE)%>% #ab_pos
  drop_na(Centiloid)%>% 
  mutate(abs_diff=as.numeric(abs(ab_Date-psy_Date)), abs_diff = replace_na(abs_diff, 0))%>% 
  group_by(ID, ab_Date)%>% 
  slice_min(abs_diff)%>% 
  ungroup()%>% 
  filter(abs_diff <= (365.25*3))%>% 
  select(-abs_diff)

ab$na_count <- apply(ab [, c(16:22)], 1, function(x) sum(is.na(x)))
ab <- ab%>% 
  arrange(na_count)%>% 
  group_by(ID, psy_Date)%>% 
  slice_min(na_count)%>% 
  ungroup()%>% 
  arrange(ab_Date)%>% 
  group_by(ID)%>% 
  slice_min(ab_Date)%>% ungroup()%>% mutate(sex = as.numeric(str_replace_all(sex, c("M" = "1", "F" = "2"))))

tau <- merge(tau, full_data, by = "ID", all.x = TRUE)%>% #ab_pos
  drop_na(Tauopathy)%>% 
  mutate(abs_diff=as.numeric(abs(tau_Date-psy_Date)), abs_diff = replace_na(abs_diff, 0))%>% 
  group_by(ID, tau_Date)%>% 
  slice_min(abs_diff)%>% 
  ungroup()%>% 
  filter(abs_diff <= (365.25*3))%>% 
  select(-abs_diff)

tau$na_count <- apply(tau [, c(16:22)], 1, function(x) sum(is.na(x)))
tau <- tau%>% 
  arrange(na_count)%>% 
  group_by(ID, psy_Date)%>% 
  slice_min(na_count)%>% 
  ungroup()%>% 
  arrange(tau_Date)%>% 
  group_by(ID)%>% 
  slice_min(tau_Date)%>% ungroup()%>% mutate(sex = as.numeric(str_replace_all(sex, c("M" = "1", "F" = "2"))))

mri <- merge(mri, full_data, by = "ID", all.x = TRUE)%>% #ab_pos
  drop_na(CortSig)%>% 
  mutate(abs_diff=as.numeric(abs(mri_Date-psy_Date)), abs_diff = replace_na(abs_diff, 0))%>% 
  group_by(ID, mri_Date)%>% 
  slice_min(abs_diff)%>% 
  ungroup()%>% 
  filter(abs_diff <= (365.25*3))%>% 
  select(-abs_diff)

mri$na_count <- apply(mri [, c(14:20)], 1, function(x) sum(is.na(x)))
mri <- mri%>% 
  arrange(na_count)%>% 
  group_by(ID, psy_Date)%>% 
  slice_min(na_count)%>% 
  ungroup()%>% 
  arrange(mri_Date)%>% 
  group_by(ID)%>% 
  slice_min(mri_Date)%>% ungroup()%>% mutate(sex = as.numeric(str_replace_all(sex, c("M" = "1", "F" = "2"))))

#keep(full_data, scaleFUN, ab, mri, tau, sure = TRUE)

#### FIGURE 3: CORRELATIONS w AGE####
# data is most recent visit, cdr_bin == 0 only
data_corr <- full_data%>% 
  filter (cdr_bin == 0)%>% 
  group_by(ID)%>% arrange(psy_Date)%>%
  slice_min(psy_Date)%>% ungroup()%>% mutate(AgeS = scale(Age), EDUCS = scale (EDUC, scale = FALSE))%>% 
  drop_na(AgeS, sex, EDUCS, apoe)
write.csv(data_corr, "Manuscript/demographics_age_correlation.csv")

# set up the partial correlations to get the values for the line plots 
data_corr <- data_corr %>%
  gather (Factor, Cognition, F1:A1)%>%
  select(ID, Age, AgeS, sex, EDUCS, apoe, Factor, Cognition)

F1 <- data_corr %>% 
  filter(Factor == "F1")%>% 
  drop_na(Cognition)
F1 <- pcor.test(F1$AgeS, F1$Cognition, list(F1$EDUCS, F1$sex, F1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F1$p.value; r <- F1$estimate; CI <- F1$conf.int[1:2]; F <- "F1"
F1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F1 <- F1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F2 <- data_corr %>% 
  filter(Factor == "F2")%>% 
  drop_na(Cognition)
F2 <- pcor.test(F2$AgeS, F2$Cognition, list(F2$EDUCS, F2$sex, F2$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F2$p.value; r <- F2$estimate; CI <- F2$conf.int[1:2]; F <- "F2"
F2 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F2) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F2 <- F2%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F3 <- data_corr %>% 
  filter(Factor == "F3")%>% 
  drop_na(Cognition)
F3 <- pcor.test(F3$AgeS, F3$Cognition, list(F3$EDUCS, F3$sex, F3$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F3$p.value; r <- F3$estimate; CI <- F3$conf.int[1:2]; F <- "F3"
F3 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F3) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F3 <- F3%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F4 <- data_corr %>% 
  filter(Factor == "F4")%>% 
  drop_na(Cognition)
F4 <- pcor.test(F4$AgeS, F4$Cognition, list(F4$EDUCS, F4$sex, F4$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F4$p.value; r <- F4$estimate; CI <- F4$conf.int[1:2]; F <- "F4"
F4 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F4) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F4 <- F4%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

K1 <- data_corr %>% 
  filter(Factor == "K1")%>% 
  drop_na(Cognition)
K1 <- pcor.test(K1$AgeS, K1$Cognition, list(K1$EDUCS, K1$sex, K1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- K1$p.value; r <- K1$estimate; CI <- K1$conf.int[1:2]; F <- "K1"
K1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(K1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
K1 <- K1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

A1 <- data_corr %>% 
  filter(Factor == "A1")%>% 
  drop_na(Cognition)
A1 <- pcor.test(A1$AgeS, A1$Cognition, list(A1$EDUCS, A1$sex, A1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- A1$p.value; r <- A1$estimate; CI <- A1$conf.int[1:2]; F <- "A1"
A1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(A1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
A1 <- A1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

G1 <- data_corr %>% 
  filter(Factor == "G1")%>% 
  drop_na(Cognition)
G1 <- pcor.test(G1$AgeS, G1$Cognition, list(G1$EDUCS, G1$sex, G1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- G1$p.value; r <- G1$estimate; CI <- G1$conf.int[1:2]; F <- "G1"
G1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(G1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
G1 <- G1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

part_corrs <- rbind(F1, F2, F3, F4, K1, G1, A1)%>% 
  mutate_at(2:5, round, 3)%>% 
  `colnames<-`(c("Factor", "r", "p", "cil", "cih"))%>% 
  mutate(p = case_when (p > 0 ~ p, p == 0 ~ 0.001), text = if_else(p > 0.001, sprintf("r = %s, p = %s, \n95CI = %s:%s",r, p, cil, cih), sprintf("r = %s, p < %s, \n95CI = %s-%s",r, p, cil, cih)))
part_corrs$Factor <- factor(part_corrs$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
data_corr$Factor <- factor(data_corr$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))

data_corr %>%
  ggplot(., aes(x = Age, y= Cognition))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm", color = "red", fill = "red", linewidth = 0.7)+ 
  labs(y = "Standardized Cognitive Score", x = "Age (years)")+
  ylim(-3, 4)+
  geom_text(data = part_corrs, mapping = aes(x = -Inf, y = Inf, label = text), hjust = -0.2, vjust = 1.2, size = 3.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
 facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggsave("Manuscript/corr_plot_age.png", width=10,height =4)  
#keep(full_data, scaleFUN, ab, tau, mri, sure = TRUE)

#### FIGURE 5: CLASSIFICATION VIOLINS ####
violins_cdr_bin <- full_data%>% mutate(cdr_bin = as.factor(cdr_bin))%>% 
  group_by(ID)%>% 
  slice_min(psy_Date)%>% 
  ungroup()%>% 
  mutate(AgeS = scale(Age), EDUCS = scale (EDUC, scale = FALSE))%>% 
  gather (Factor, Cognition, F1:A1)%>% 
  select(ID, Age, AgeS, sex, EDUCS, apoe, Factor, Cognition, cdr_bin)%>% drop_na(AgeS, sex, EDUCS, apoe, cdr_bin)
write.csv(violins_cdr_bin, "Manuscript/demographics_cdr_bin_violin.csv")

F1 <- violins_cdr_bin %>% 
  filter(Factor == "F1")%>% drop_na(Cognition)
F1 <- aov(Cognition ~ cdr_bin + AgeS + sex + EDUCS + apoe, data = F1)
F1_store <- summary(F1)[[1]][1,4:5]; F1_store <- F1_store%>% mutate(Factor = "F1", df = "1")
F1_store$dfresid <- F1$df.residual

F2 <- violins_cdr_bin %>% 
  filter(Factor == "F2")%>% drop_na(Cognition)
F2 <- aov(Cognition ~ cdr_bin + AgeS + sex + EDUCS + apoe, data = F2)
F2_store <- summary(F2)[[1]][1,4:5]; F2_store <- F2_store%>% mutate(Factor = "F2", df = "1")
F2_store$dfresid <- F2$df.residual

F3 <- violins_cdr_bin %>% 
  filter(Factor == "F3")%>% drop_na(Cognition)
F3 <- aov(Cognition ~ cdr_bin + AgeS + sex + EDUCS + apoe, data = F3)
F3_store <- summary(F3)[[1]][1,4:5]; F3_store <- F3_store%>% mutate(Factor = "F3", df = "1")
F3_store$dfresid <- F3$df.residual

F4 <- violins_cdr_bin %>% 
  filter(Factor == "F4")%>% drop_na(Cognition)
F4 <- aov(Cognition ~ cdr_bin + AgeS + sex + EDUCS + apoe, data = F4)
F4_store <- summary(F4)[[1]][1,4:5]; F4_store <- F4_store%>% mutate(Factor = "F4", df = "1")
F4_store$dfresid <- F4$df.residual

K1 <- violins_cdr_bin %>% 
  filter(Factor == "K1")%>% drop_na(Cognition)
K1 <- aov(Cognition ~ cdr_bin + AgeS + sex + EDUCS + apoe, data = K1)
K1_store <- summary(K1)[[1]][1,4:5]; K1_store <- K1_store%>% mutate(Factor = "K1", df = "1")
K1_store$dfresid <- K1$df.residual

G1 <- violins_cdr_bin %>% 
  filter(Factor == "G1")%>% drop_na(Cognition)
G1 <- aov(Cognition ~ cdr_bin + AgeS + sex + EDUCS + apoe, data = G1)
G1_store <- summary(G1)[[1]][1,4:5]; G1_store <- G1_store%>% mutate(Factor = "G1", df = "1")
G1_store$dfresid <- G1$df.residual

A1 <- violins_cdr_bin %>% 
  filter(Factor == "A1")%>% drop_na(Cognition)
A1 <- aov(Cognition ~ cdr_bin + AgeS + sex + EDUCS + apoe, data = A1)
A1_store <- summary(A1)[[1]][1,4:5]; A1_store <- A1_store%>% mutate(Factor = "A1", df = "1")
A1_store$dfresid <- A1$df.residual

all_anovas <- rbind(F1_store, F2_store, F3_store, F4_store, K1_store, G1_store, A1_store)%>% 
  mutate_at(1:2, round, 3)%>% 
  `colnames<-`(c("F", "p", "Factor", "df", "dfr"))%>% 
  mutate(p = case_when (p > 0 ~ p, p == 0 ~ 0.001), text = if_else(p > 0.001, sprintf("F = %s, p = %s",F, p), sprintf("F = %s, p < %s",F, p)))%>% 
  mutate(cdr_bin = "1")
all_anovas$Factor <- factor(all_anovas$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
violins_cdr_bin$Factor <- factor(violins_cdr_bin$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))


cdr_bin_violin <- violins_cdr_bin %>% 
  ggplot(., aes(x = cdr_bin, y = Cognition, group=cdr_bin))+
  geom_violin(trim = FALSE , alpha = 0.5, aes(fill = cdr_bin, color=cdr_bin), linewidth = 1)+
  geom_point(aes(y=Cognition),size =1,  show.legend = F, position = position_jitter(width = .03), alpha = 0.03)+
  geom_boxplot(width = .15, size = 0.5,  show.legend = F, outlier.shape = NA, aes(fill = cdr_bin), alpha = 0.5)+
  scale_color_manual(values = c("#EBB261","#5A4A6F"), labels=c('Cognitively Normal', 'Cognitively Impaired'))+
  scale_fill_manual(values = c("#EBB261","#5A4A6F"), labels=c('Cognitively Normal', 'Cognitively Impaired'))+
  scale_x_discrete(labels = c("", "")) +
  geom_text(data = all_anovas, mapping = aes(x = -Inf, y = Inf, label = text), hjust = -0.2, vjust = 1.2)+
  ylim(-3, 3.5)+
  labs (x ="", y = "", color = "Clinical Dementia Rating", fill = "Clinical Dementia Rating") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
cdr_bin_violin

violins_abpos <- full_data%>% mutate(abpos = as.factor(abpos))%>% drop_na(abpos)%>%
  group_by(ID)%>% 
  slice_min(psy_Date)%>% 
  ungroup()%>% 
  mutate(AgeS = scale(Age), EDUCS = scale (EDUC, scale = FALSE))%>% 
  gather (Factor, Cognition, F1:A1)%>% 
  select(ID, Age, AgeS, sex, EDUCS, apoe, Factor, Cognition, cdr_bin, abpos)%>% drop_na(AgeS, sex, EDUCS, apoe, cdr_bin)
write.csv(violins_abpos, "Manuscript/demographics_abpos_violin.csv")

F1 <- violins_abpos %>% 
  filter(Factor == "F1")%>% drop_na(Cognition)
F1 <- aov(Cognition ~ abpos + AgeS + sex + EDUCS + apoe, data = F1)
F1_store <- summary(F1)[[1]][1,4:5]; F1_store <- F1_store%>% mutate(Factor = "F1", df = "1")
F1_store$dfresid <- F1$df.residual

F2 <- violins_abpos %>% 
  filter(Factor == "F2")%>% drop_na(Cognition)
F2 <- aov(Cognition ~ abpos + AgeS + sex + EDUCS + apoe, data = F2)
F2_store <- summary(F2)[[1]][1,4:5]; F2_store <- F2_store%>% mutate(Factor = "F2", df = "1")
F2_store$dfresid <- F2$df.residual

F3 <- violins_abpos %>% 
  filter(Factor == "F3")%>% drop_na(Cognition)
F3 <- aov(Cognition ~ abpos + AgeS + sex + EDUCS + apoe, data = F3)
F3_store <- summary(F3)[[1]][1,4:5]; F3_store <- F3_store%>% mutate(Factor = "F3", df = "1")
F3_store$dfresid <- F3$df.residual

F4 <- violins_abpos %>% 
  filter(Factor == "F4")%>% drop_na(Cognition)
F4 <- aov(Cognition ~ abpos + AgeS + sex + EDUCS + apoe, data = F4)
F4_store <- summary(F4)[[1]][1,4:5]; F4_store <- F4_store%>% mutate(Factor = "F4", df = "1")
F4_store$dfresid <- F4$df.residual

K1 <- violins_abpos %>% 
  filter(Factor == "K1")%>% drop_na(Cognition)
K1 <- aov(Cognition ~ abpos + AgeS + sex + EDUCS + apoe, data = K1)
K1_store <- summary(K1)[[1]][1,4:5]; K1_store <- K1_store%>% mutate(Factor = "K1", df = "1")
K1_store$dfresid <- K1$df.residual

G1 <- violins_abpos %>% 
  filter(Factor == "G1")%>% drop_na(Cognition)
G1 <- aov(Cognition ~ abpos + AgeS + sex + EDUCS + apoe, data = G1)
G1_store <- summary(G1)[[1]][1,4:5]; G1_store <- G1_store%>% mutate(Factor = "G1", df = "1")
G1_store$dfresid <- G1$df.residual

A1 <- violins_abpos %>% 
  filter(Factor == "A1")%>% drop_na(Cognition)
A1 <- aov(Cognition ~ abpos + AgeS + sex + EDUCS + apoe, data = A1)
A1_store <- summary(A1)[[1]][1,4:5]; A1_store <- A1_store%>% mutate(Factor = "A1", df = "1")
A1_store$dfresid <- A1$df.residual

all_anovas <- rbind(F1_store, F2_store, F3_store, F4_store, K1_store, G1_store, A1_store)%>% 
  mutate_at(1:2, round, 3)%>% 
  `colnames<-`(c("F", "p", "Factor", "df", "dfr"))%>% 
  mutate(p = case_when (p > 0 ~ p, p == 0 ~ 0.001), text = if_else(p > 0.001, sprintf("F = %s, p = %s",F, p), sprintf("F = %s, p < %s",F, p)))%>% 
  mutate(abpos = "1")
all_anovas$Factor <- factor(all_anovas$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
violins_abpos$Factor <- factor(violins_abpos$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))

abpos_violin <- violins_abpos %>% 
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x = abpos, y = Cognition, group=abpos))+
  geom_violin(trim = FALSE , alpha = 0.5, aes(fill = abpos, color=abpos), linewidth = 1)+
  geom_point(aes(y=Cognition),size =1,  show.legend = F, position = position_jitter(width = .03), alpha = 0.03)+
  geom_boxplot(width = .15, size = 0.5,  show.legend = F, outlier.shape = NA, aes(fill = abpos), alpha = 0.5)+
  scale_color_manual(values = c("#FC4E07", "#00AFBB"), labels=c('Amyloid Negative', 'Amyloid Positive'))+
  scale_fill_manual(values = c( "#FC4E07", "#00AFBB"), labels=c('Amyloid Negative', 'Amyloid Positive'))+
  scale_x_discrete(labels = c("", "")) +
  geom_text(data = all_anovas, mapping = aes(x = -Inf, y = Inf, label = text), hjust = -0.2, vjust = 1.2)+
  ylim(-3, 3.5)+
  labs (x ="", y = "Standardized Cognitive Score", color = "Amyloid-Positivity", fill = "Amyloid-Positivity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
abpos_violin

violins_taupos <- full_data%>% mutate(taupos = as.factor(taupos))%>% drop_na(taupos)%>%
  group_by(ID)%>% 
  slice_min(psy_Date)%>% 
  ungroup()%>% 
  mutate(AgeS = scale(Age), EDUCS = scale (EDUC, scale = FALSE))%>% 
  gather (Factor, Cognition, F1:A1)%>% 
  select(ID, Age, AgeS, sex, EDUCS, apoe, Factor, Cognition, cdr_bin, taupos)%>% drop_na(AgeS, sex, EDUCS, apoe, cdr_bin)
write.csv(violins_taupos, "Manuscript/demographics_taupos_violin.csv")

F1 <- violins_taupos %>% 
  filter(Factor == "F1")%>% drop_na(Cognition)
F1 <- aov(Cognition ~ taupos + AgeS + sex + EDUCS + apoe, data = F1)
F1_store <- summary(F1)[[1]][1,4:5]; F1_store <- F1_store%>% mutate(Factor = "F1", df = "1")
F1_store$dfresid <- F1$df.residual

F2 <- violins_taupos %>% 
  filter(Factor == "F2")%>% drop_na(Cognition)
F2 <- aov(Cognition ~ taupos + AgeS + sex + EDUCS + apoe, data = F2)
F2_store <- summary(F2)[[1]][1,4:5]; F2_store <- F2_store%>% mutate(Factor = "F2", df = "1")
F2_store$dfresid <- F2$df.residual

F3 <- violins_taupos %>% 
  filter(Factor == "F3")%>% drop_na(Cognition)
F3 <- aov(Cognition ~ taupos + AgeS + sex + EDUCS + apoe, data = F3)
F3_store <- summary(F3)[[1]][1,4:5]; F3_store <- F3_store%>% mutate(Factor = "F3", df = "1")
F3_store$dfresid <- F3$df.residual

F4 <- violins_taupos %>% 
  filter(Factor == "F4")%>% drop_na(Cognition)
F4 <- aov(Cognition ~ taupos + AgeS + sex + EDUCS + apoe, data = F4)
F4_store <- summary(F4)[[1]][1,4:5]; F4_store <- F4_store%>% mutate(Factor = "F4", df = "1")
F4_store$dfresid <- F4$df.residual

K1 <- violins_taupos %>% 
  filter(Factor == "K1")%>% drop_na(Cognition)
K1 <- aov(Cognition ~ taupos + AgeS + sex + EDUCS + apoe, data = K1)
K1_store <- summary(K1)[[1]][1,4:5]; K1_store <- K1_store%>% mutate(Factor = "K1", df = "1")
K1_store$dfresid <- K1$df.residual

G1 <- violins_taupos %>% 
  filter(Factor == "G1")%>% drop_na(Cognition)
G1 <- aov(Cognition ~ taupos + AgeS + sex + EDUCS + apoe, data = G1)
G1_store <- summary(G1)[[1]][1,4:5]; G1_store <- G1_store%>% mutate(Factor = "G1", df = "1")
G1_store$dfresid <- G1$df.residual

A1 <- violins_taupos %>% 
  filter(Factor == "A1")%>% drop_na(Cognition)
A1 <- aov(Cognition ~ taupos + AgeS + sex + EDUCS + apoe, data = A1)
A1_store <- summary(A1)[[1]][1,4:5]; A1_store <- A1_store%>% mutate(Factor = "A1", df = "1")
A1_store$dfresid <- A1$df.residual

all_anovas <- rbind(F1_store, F2_store, F3_store, F4_store, K1_store, G1_store, A1_store)%>% 
  mutate_at(1:2, round, 3)%>% 
  `colnames<-`(c("F", "p", "Factor", "df", "dfr"))%>% 
  mutate(p = case_when (p > 0 ~ p, p == 0 ~ 0.001), text = if_else(p > 0.001, sprintf("F = %s, p = %s",F, p), sprintf("F = %s, p < %s",F, p)))%>% 
  mutate(taupos = "1")
all_anovas$Factor <- factor(all_anovas$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
violins_taupos$Factor <- factor(violins_taupos$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))

taupos_violin <- violins_taupos %>%  
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x = taupos, y = Cognition, group=taupos))+
  geom_violin(trim = FALSE , alpha = 0.5, aes(fill = taupos, color=taupos), linewidth = 1)+
  geom_point(aes(y=Cognition),size =1,  show.legend = F, position = position_jitter(width = .03), alpha = 0.03)+
  geom_boxplot(width = .15, size = 0.5,  show.legend = F, outlier.shape = NA, aes(fill = taupos), alpha = 0.5)+
  scale_color_manual(values = c("#FF9900","#339966"), labels=c('Tau Negative', 'Tau Positive'))+
  scale_fill_manual(values = c( "#FF9900","#339966"), labels=c('Tau Negative', 'Tau Positive'))+
  scale_x_discrete(labels = c("", "")) +
  geom_text(data = all_anovas, mapping = aes(x = -Inf, y = Inf, label = text), hjust = -0.2, vjust = 1.2)+
  ylim(-3, 3.5)+
  labs (x ="", y = "", color = "Tau-Positivity", fill = "Tau-Positivity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
taupos_violin

ggarrange(cdr_bin_violin, abpos_violin, taupos_violin, ncol = 1)
ggsave("Manuscript/all_violins.png", width=10,height =12)  

#keep(full_data, scaleFUN, ab, tau, mri, sure = TRUE)
#### FIGURE 4: CORRELATIONS w BIOMARKERS #### 
# set data aside 
data_corr_ab <- ab%>% drop_na(Centiloid)%>% filter (cdr_bin == 0)%>% mutate(AgeS = as.numeric(scale(Age)), EDUCS = as.numeric(scale (EDUC, scale = FALSE)))%>% drop_na(AgeS, EDUCS, sex, apoe)
write.csv(data_corr_ab, "Manuscript/demographics_corr_ab.csv")
data_corr_tau <- tau%>% drop_na(Tauopathy)%>% filter (cdr_bin == 0)%>% mutate(AgeS = scale(Age), EDUCS = scale (EDUC, scale = FALSE))%>% drop_na(AgeS, EDUCS, sex, apoe)
write.csv(data_corr_tau, "Manuscript/demographics_corr_tau.csv")
data_corr_mri <- mri%>% drop_na(CortSig)%>% filter (cdr_bin == 0)%>% mutate(AgeS = scale(Age), EDUCS = scale (EDUC, scale = FALSE))%>% drop_na(AgeS, EDUCS, sex, apoe)
write.csv(data_corr_mri, "Manuscript/demographics_corr_mri.csv")

# line plots with partial correlations added 
data_corr_ab <- data_corr_ab %>%
  gather (Factor, Cognition, F1:A1)%>%
  select(ID, AgeS, sex, EDUCS, apoe, Factor, Cognition, Centiloid)

F1 <- data_corr_ab %>% filter(Factor == "F1")%>% drop_na(Cognition)
F1 <- pcor.test(F1$Centiloid, F1$Cognition, list(F1$AgeS, F1$EDUCS, F1$sex, F1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F1$p.value; r <- F1$estimate; CI <- F1$conf.int[1:2]; F <- "F1"
F1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F1 <- F1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F2 <- data_corr_ab %>% filter(Factor == "F2")%>% drop_na(Cognition)
F2 <- pcor.test(F2$Centiloid, F2$Cognition, list(F2$AgeS, F2$EDUCS, F2$sex, F2$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F2$p.value; r <- F2$estimate; CI <- F2$conf.int[1:2]; F <- "F2"
F2 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F2) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F2 <- F2%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F3 <- data_corr_ab %>% filter(Factor == "F3")%>% drop_na(Cognition)
F3 <- pcor.test(F3$Centiloid, F3$Cognition, list(F3$AgeS, F3$EDUCS, F3$sex, F3$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F3$p.value; r <- F3$estimate; CI <- F3$conf.int[1:2]; F <- "F3"
F3 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F3) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F3 <- F3%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F4 <- data_corr_ab %>% filter(Factor == "F4")%>% drop_na(Cognition)
F4 <- pcor.test(F4$Centiloid, F4$Cognition, list(F4$AgeS, F4$EDUCS, F4$sex, F4$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F4$p.value; r <- F4$estimate; CI <- F4$conf.int[1:2]; F <- "F4"
F4 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F4) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F4 <- F4%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

K1 <- data_corr_ab %>% filter(Factor == "K1")%>% drop_na(Cognition)
K1 <- pcor.test(K1$Centiloid, K1$Cognition, list(K1$AgeS, K1$EDUCS, K1$sex, K1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- K1$p.value; r <- K1$estimate; CI <- K1$conf.int[1:2]; F <- "K1"
K1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(K1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
K1 <- K1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

G1 <- data_corr_ab %>% filter(Factor == "G1")%>% drop_na(Cognition)
G1 <- pcor.test(G1$Centiloid, G1$Cognition, list(G1$AgeS, G1$EDUCS, G1$sex, G1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- G1$p.value; r <- G1$estimate; CI <- G1$conf.int[1:2]; F <- "G1"
G1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(G1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
G1 <- G1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

A1 <- data_corr_ab %>% filter(Factor == "A1")%>% drop_na(Cognition)
A1 <- pcor.test(A1$Centiloid, A1$Cognition, list(A1$AgeS, A1$EDUCS, A1$sex, A1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- A1$p.value; r <- A1$estimate; CI <- A1$conf.int[1:2]; F <- "A1"
A1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(A1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
A1 <- A1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

part_corrs <- rbind(F1, F2, F3, F4, K1, G1, A1)%>% 
  mutate_at(2:5, round, 3)%>% 
  `colnames<-`(c("Factor", "r", "p", "cil", "cih"))%>% 
  mutate(p = case_when (p > 0 ~ p, p == 0 ~ 0.001), text = if_else(p > 0.001, sprintf("r = %s, p = %s, \n95CI = %s:%s",r, p, cil, cih), sprintf("r = %s, p < %s, \n95CI = %s-%s",r, p, cil, cih)))
part_corrs$Factor <- factor(part_corrs$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
data_corr_ab$Factor <- factor(data_corr_ab$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))

corr_plot_ab <- data_corr_ab %>%
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x = Centiloid, y= Cognition))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm", color = "red", fill = "red", linewidth = 0.7)+ 
  labs(y = "Cog Score", x = "Amyloid PET (Centiloid)")+
  ylim(-3, 4)+
  geom_text(data = part_corrs, mapping = aes(x = -Inf, y = Inf, label = text), hjust = -0.2, vjust = 1.2, size = 3.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
corr_plot_ab 

data_corr_tau <- data_corr_tau %>%
  gather (Factor, Cognition, F1:A1)%>%
  select(ID, AgeS, sex, EDUCS, apoe, Factor, Cognition, Tauopathy)

F1 <- data_corr_tau %>% filter(Factor == "F1")%>% drop_na(Cognition)
F1 <- pcor.test(F1$Tauopathy, F1$Cognition, list(F1$AgeS, F1$EDUCS, F1$sex, F1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F1$p.value; r <- F1$estimate; CI <- F1$conf.int[1:2]; F <- "F1"
F1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F1 <- F1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F2 <- data_corr_tau %>% filter(Factor == "F2")%>% drop_na(Cognition)
F2 <- pcor.test(F2$Tauopathy, F2$Cognition, list(F2$AgeS, F2$EDUCS, F2$sex, F2$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F2$p.value; r <- F2$estimate; CI <- F2$conf.int[1:2]; F <- "F2"
F2 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F2) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F2 <- F2%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F3 <- data_corr_tau %>% filter(Factor == "F3")%>% drop_na(Cognition)
F3 <- pcor.test(F3$Tauopathy, F3$Cognition, list(F3$AgeS, F3$EDUCS, F3$sex, F3$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F3$p.value; r <- F3$estimate; CI <- F3$conf.int[1:2]; F <- "F3"
F3 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F3) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F3 <- F3%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F4 <- data_corr_tau %>% filter(Factor == "F4")%>% drop_na(Cognition)
F4 <- pcor.test(F4$Tauopathy, F4$Cognition, list(F4$AgeS, F4$EDUCS, F4$sex, F4$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F4$p.value; r <- F4$estimate; CI <- F4$conf.int[1:2]; F <- "F4"
F4 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F4) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F4 <- F4%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

K1 <- data_corr_tau %>% filter(Factor == "K1")%>% drop_na(Cognition)
K1 <- pcor.test(K1$Tauopathy, K1$Cognition, list(K1$AgeS, K1$EDUCS, K1$sex, K1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- K1$p.value; r <- K1$estimate; CI <- K1$conf.int[1:2]; F <- "K1"
K1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(K1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
K1 <- K1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

G1 <- data_corr_tau %>% filter(Factor == "G1")%>% drop_na(Cognition)
G1 <- pcor.test(G1$Tauopathy, G1$Cognition, list(G1$AgeS, G1$EDUCS, G1$sex, G1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- G1$p.value; r <- G1$estimate; CI <- G1$conf.int[1:2]; F <- "G1"
G1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(G1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
G1 <- G1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

A1 <- data_corr_tau %>% filter(Factor == "A1")%>% drop_na(Cognition)
A1 <- pcor.test(A1$Tauopathy, A1$Cognition, list(A1$AgeS, A1$EDUCS, A1$sex, A1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- A1$p.value; r <- A1$estimate; CI <- A1$conf.int[1:2]; F <- "A1"
A1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(A1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
A1 <- A1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

part_corrs <- rbind(F1, F2, F3, F4, K1, G1, A1)%>% 
  mutate_at(2:5, round, 3)%>% 
  `colnames<-`(c("Factor", "r", "p", "cil", "cih"))%>% 
  mutate(p = case_when (p > 0 ~ p, p == 0 ~ 0.001), text = if_else(p > 0.001, sprintf("r = %s, p = %s, \n95CI = %s:%s",r, p, cil, cih), sprintf("r = %s, p < %s, \n95CI = %s-%s",r, p, cil, cih)))
part_corrs$Factor <- factor(part_corrs$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
data_corr_tau$Factor <- factor(data_corr_tau$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))

corr_plot_tau <- data_corr_tau %>%
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x = Tauopathy, y= Cognition))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm", color = "red", fill = "red", linewidth = 0.7)+ 
  labs(y = "Standardized Cognitive Score", x = "Tau PET (SUVR)")+
  ylim(-3, 4)+
  geom_text(data = part_corrs, mapping = aes(x = -Inf, y = Inf, label = text), hjust = -0.2, vjust = 1.2, size = 3.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
corr_plot_tau 

data_corr_mri <- data_corr_mri %>%
  gather (Factor, Cognition, F1:A1)%>%
  select(ID, AgeS, sex, EDUCS, apoe, Factor, Cognition, CortSig)

F1 <- data_corr_mri %>% filter(Factor == "F1")%>% drop_na(Cognition)
F1 <- pcor.test(F1$CortSig, F1$Cognition, list(F1$AgeS, F1$EDUCS, F1$sex, F1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F1$p.value; r <- F1$estimate; CI <- F1$conf.int[1:2]; F <- "F1"
F1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F1 <- F1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F2 <- data_corr_mri %>% filter(Factor == "F2")%>% drop_na(Cognition)
F2 <- pcor.test(F2$CortSig, F2$Cognition, list(F2$AgeS, F2$EDUCS, F2$sex, F2$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F2$p.value; r <- F2$estimate; CI <- F2$conf.int[1:2]; F <- "F2"
F2 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F2) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F2 <- F2%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F3 <- data_corr_mri %>% filter(Factor == "F3")%>% drop_na(Cognition)
F3 <- pcor.test(F3$CortSig, F3$Cognition, list(F3$AgeS, F3$EDUCS, F3$sex, F3$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F3$p.value; r <- F3$estimate; CI <- F3$conf.int[1:2]; F <- "F3"
F3 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F3) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F3 <- F3%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

F4 <- data_corr_mri %>% filter(Factor == "F4")%>% drop_na(Cognition)
F4 <- pcor.test(F4$CortSig, F4$Cognition, list(F4$AgeS, F4$EDUCS, F4$sex, F4$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- F4$p.value; r <- F4$estimate; CI <- F4$conf.int[1:2]; F <- "F4"
F4 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(F4) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
F4 <- F4%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

G1 <- data_corr_mri %>% filter(Factor == "G1")%>% drop_na(Cognition)
G1 <- pcor.test(G1$CortSig, G1$Cognition, list(G1$AgeS, G1$EDUCS, G1$sex, G1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- G1$p.value; r <- G1$estimate; CI <- G1$conf.int[1:2]; F <- "G1"
G1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(G1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
G1 <- G1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

K1 <- data_corr_mri %>% filter(Factor == "K1")%>% drop_na(Cognition)
K1 <- pcor.test(K1$CortSig, K1$Cognition, list(K1$AgeS, K1$EDUCS, K1$sex, K1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- K1$p.value; r <- K1$estimate; CI <- K1$conf.int[1:2]; F <- "K1"
K1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(K1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
K1 <- K1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

A1 <- data_corr_mri %>% filter(Factor == "A1")%>% drop_na(Cognition)
A1 <- pcor.test(A1$CortSig, A1$Cognition, list(A1$AgeS, A1$EDUCS, A1$sex, A1$apoe), semi = FALSE, conf.level = 0.95, nrep = 1000, method = "pearson")
p <- A1$p.value; r <- A1$estimate; CI <- A1$conf.int[1:2]; F <- "A1"
A1 <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(A1) <- c("Factor", "estimate", "p.value", "CI.Low", "CI.High")
A1 <- A1%>% mutate(Factor = F, CI.Low = CI[1], CI.High = CI[2], p.value = p, estimate = r)

part_corrs <- rbind(F1, F2, F3, F4, K1, G1, A1)%>% 
  mutate_at(2:5, round, 3)%>% 
  `colnames<-`(c("Factor", "r", "p", "cil", "cih"))%>% 
  mutate(p = case_when (p > 0 ~ p, p == 0 ~ 0.001), text = if_else(p > 0.001, sprintf("r = %s, p = %s, \n95CI = %s:%s",r, p, cil, cih), sprintf("r = %s, p < %s, \n95CI = %s-%s",r, p, cil, cih)))
part_corrs$Factor <- factor(part_corrs$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
data_corr_mri$Factor <- factor(data_corr_mri$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))

corr_plot_mri <- data_corr_mri %>%
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x = CortSig, y= Cognition))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm", color = "red", fill = "red", linewidth = 0.7)+ 
  labs(y = "Cog Score", x = "Cortical Thickness (mm)")+
  ylim(-3, 4)+
  geom_text(data = part_corrs, mapping = aes(x = -Inf, y = Inf, label = text), hjust = -0.2, vjust = 1.2, size = 3.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
corr_plot_mri

corr_fig <- ggarrange(corr_plot_ab + rremove("ylab"), corr_plot_tau + rremove("ylab"), corr_plot_mri + rremove("ylab"), ncol = 1)

annotate_figure(corr_fig, left = text_grob("Standardized Cognitive Score", rot = 90))

ggsave("Manuscript/all_corrs.png", width=10,height =12)  
#### FIGURE 7: TEST RETEST ####
retest <- full_data%>% filter(cdr_bin==0)%>% group_by(ID)%>% arrange(psy_Date)%>% mutate(maxtime = max(Time), count = "1", visit = cumsum(count), maxvisit = max(visit))%>% ungroup()%>% filter(maxtime >= 4*365.25 & maxvisit >= 4)
write.csv(retest, 'Manuscript/reliability_demographics.csv')

F1_5 <- retest%>% drop_na(F1)%>% select(ID, F1, visit)%>% pivot_wider(names_from = visit, values_from = F1)
F2_5 <- retest%>% select(ID, F2, visit)%>% pivot_wider(names_from = visit, values_from = F2)%>% filter(ID %in% F1_5$ID)%>% select(-ID)
F3_5 <- retest%>% select(ID, F3, visit)%>% pivot_wider(names_from = visit, values_from = F3)%>% filter(ID %in% F1_5$ID)%>% select(-ID)
F4_5 <- retest%>% select(ID, F4, visit)%>% pivot_wider(names_from = visit, values_from = F4)%>% filter(ID %in% F1_5$ID)%>% select(-ID)
K1_5 <- retest%>% select(ID, K1, visit)%>% pivot_wider(names_from = visit, values_from = K1)%>% filter(ID %in% F1_5$ID)%>% select(-ID)
G1_5 <- retest%>% select(ID, G1, visit)%>% pivot_wider(names_from = visit, values_from = G1)%>% filter(ID %in% F1_5$ID)%>% select(-ID)
A1_5 <- retest%>% select(ID, A1, visit)%>% pivot_wider(names_from = visit, values_from = A1)%>% filter(ID %in% F1_5$ID)%>% select(-ID)

F1_5 <- F1_5%>% select(-ID)

ICC(F1_5)
ICC(F2_5)
ICC(F3_5)
ICC(F4_5)
ICC(K1_5)
ICC(G1_5)
ICC(A1_5)

# Make reliability figure V1 vs V2, 3, 4, and max
V1V2 <- retest%>%
  filter(visit == 1 | visit == 2)%>%
  gather (Factor, Cognition, F1:A1)%>%
  select(ID, visit, Cognition, Factor)%>% drop_na(Cognition)%>%
  pivot_wider(names_from = visit, values_from = Cognition)
colnames(V1V2) <- c("ID", "Factor", "V1", "V2")

plot1 <- V1V2 %>%
  drop_na(V1, V2)%>%
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x = V1, y= V2))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm", color = "red", fill = "red", size = 0.7)+ 
  labs(y = "Visit 2 Standardized Cognitive Score \n~one year delay", x = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., nrow = 2, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
plot1
ggsave("Manuscript/cog_reliability.png", width=10,height =4)  

V1V4 <- retest%>%
  filter(visit == 1 | visit == 4)%>%
  gather (Factor, Cognition, F1:A1)%>%
  select(ID, visit, Cognition, Factor)%>% drop_na(Cognition)%>%
  pivot_wider(names_from = visit, values_from = Cognition)
colnames(V1V4) <- c("ID", "Factor", "V1", "V4")

plot2 <- V1V4 %>%
  drop_na(V1, V4)%>%
  mutate(across(Factor, ~factor(., levels=c("F1","F2","F3","F4","K1","A1","G1"))))%>%
  ggplot(., aes(x = V1, y= V4))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm", color = "red", fill = "red", size = 0.7)+ 
  labs(y = "Visit 4 Standardized Cognitive Score \n~three year delay", x = "Visit 1 Standardized Cognitive Score")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., nrow = 2, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggarrange(plot1, plot2, nrow = 2)
ggsave("Manuscript/cog_reliability2b.png", width=10,height =8)  

#keep(full_data, scaleFUN, ab, tau, mri, sure = TRUE)


#### FIGURE 6: BIOMARKER LONGITUDINAL COGNITION ####
pacman::p_load(lme4)
lmerCont <- lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun=1000000))

ab <- ab%>% distinct(ID, ab_Date, Centiloid, BLAge, BLcdr_bin)%>% mutate(CentiloidS = as.numeric(scale(Centiloid)), baseline_ageS = as.numeric(scale(BLAge)))
tau <- tau%>% distinct(ID, tau_Date, Tauopathy, BLAge, BLcdr_bin)%>% mutate(TauopathyS = as.numeric(scale(Tauopathy)), baseline_ageS = as.numeric(scale(BLAge)))
mri <- mri%>% distinct(ID, mri_Date, CortSig, BLAge, BLcdr_bin)%>% mutate(CortSigS = as.numeric(scale(CortSig)), baseline_ageS = as.numeric(scale(BLAge)))

cog_ab <- merge(ab, full_data, by = "ID", all = TRUE)%>% 
  filter(psy_Date > ab_Date)%>% 
  group_by(ID)%>% 
  arrange(psy_Date)%>% 
  mutate(count = "1", visit = cumsum(count), max_visit = max(visit), Time = as.numeric(psy_Date-ab_Date), Time = Time / 365.25)%>% 
  ungroup()%>% mutate(EDUCS = as.numeric(scale(EDUC)))%>% filter(max_visit > 1)%>% filter(Time < 6)
write.csv(cog_ab, "Manuscript/demographics_ab_cog.csv")

cog_tau <- merge(tau, full_data, by = "ID", all = TRUE)%>% 
  filter(psy_Date > tau_Date)%>% 
  group_by(ID)%>% 
  arrange(psy_Date)%>% 
  mutate(count = "1", visit = cumsum(count), max_visit = max(visit), Time = as.numeric(psy_Date-tau_Date), Time = Time / 365.25)%>% 
  ungroup()%>% mutate(EDUCS = as.numeric(scale(EDUC)))%>% filter(max_visit > 1)%>% filter(Time < 6)
write.csv(cog_tau, "Manuscript/demographics_tau_cog.csv")

cog_mri <- merge(mri, full_data, by = "ID", all = TRUE)%>% 
  filter(psy_Date > mri_Date)%>% 
  group_by(ID)%>% 
  arrange(psy_Date)%>% 
  mutate(count = "1", visit = cumsum(count), max_visit = max(visit), Time = as.numeric(psy_Date-mri_Date), Time = Time / 365.25)%>% 
  ungroup()%>% mutate(EDUCS = as.numeric(scale(EDUC)))%>% filter(max_visit > 1)%>% filter(Time < 6)
write.csv(cog_mri, "Manuscript/demographics_mri_cog.csv")

F1 <- lmer(F1 ~ 1 + Time * CentiloidS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_ab, na.action = na.omit, control = lmerCont)
F2 <- lmer(F2 ~ 1 + Time * CentiloidS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_ab, na.action = na.omit, control = lmerCont)
F3 <- lmer(F3 ~ 1 + Time * CentiloidS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_ab, na.action = na.omit, control = lmerCont)
F4 <- lmer(F4 ~ 1 + Time * CentiloidS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_ab, na.action = na.omit, control = lmerCont)
K1 <- lmer(K1 ~ 1 + Time * CentiloidS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_ab, na.action = na.omit, control = lmerCont)
A1 <- lmer(A1 ~ 1 + Time * CentiloidS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_ab, na.action = na.omit, control = lmerCont)
G1 <- lmer(G1 ~ 1 + Time * CentiloidS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_ab, na.action = na.omit, control = lmerCont)

sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T, title = "Baseline Centiloid and Longitudinal Cognition", file = "Manuscript/baselineABPET_predicts_longcog.doc")
sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T)
sjPlot::get_model_data(F1, type = "est", terms = "Time:CentiloidS")[9]
sjPlot::get_model_data(F2, type = "est", terms = "Time:CentiloidS")[9]
sjPlot::get_model_data(F3, type = "est", terms = "Time:CentiloidS")[9]
sjPlot::get_model_data(F4, type = "est", terms = "Time:CentiloidS")[9]
sjPlot::get_model_data(K1, type = "est", terms = "Time:CentiloidS")[9]
sjPlot::get_model_data(G1, type = "est", terms = "Time:CentiloidS")[9]
sjPlot::get_model_data(A1, type = "est", terms = "Time:CentiloidS")[9]

# Collect data 
F1mod <- data.frame(sjPlot::get_model_data(F1, type = "pred", terms = c("Time", "CentiloidS"), mdrt.values = "meansd"))%>% mutate(F1pred = predicted, F1stderr = std.error, F1CIH = conf.high, F1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F1pred, F1stderr, F1CIL, F1CIH, Year, Group)

F2mod <- data.frame(sjPlot::get_model_data(F2, type = "pred", terms = c("Time", "CentiloidS"), mdrt.values = "meansd"))%>% mutate(F2pred = predicted, F2stderr = std.error, F2CIH = conf.high, F2CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F2pred, F2stderr, F2CIL, F2CIH, Year, Group)

F3mod <- data.frame(sjPlot::get_model_data(F3, type = "pred", terms = c("Time", "CentiloidS"), mdrt.values = "meansd"))%>% mutate(F3pred = predicted, F3stderr = std.error, F3CIH = conf.high, F3CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F3pred, F3stderr, F3CIL, F3CIH, Year, Group)

F4mod <- data.frame(sjPlot::get_model_data(F4, type = "pred", terms = c("Time", "CentiloidS"), mdrt.values = "meansd"))%>% mutate(F4pred = predicted, F4stderr = std.error, F4CIH = conf.high, F4CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F4pred, F4stderr, F4CIL, F4CIH, Year, Group)

K1mod <- data.frame(sjPlot::get_model_data(K1, type = "pred", terms = c("Time", "CentiloidS"), mdrt.values = "meansd"))%>% mutate(K1pred = predicted, K1stderr = std.error, K1CIH = conf.high, K1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(K1pred, K1stderr, K1CIL, K1CIH, Year, Group)

G1mod <- data.frame(sjPlot::get_model_data(G1, type = "pred", terms = c("Time", "CentiloidS"), mdrt.values = "meansd"))%>% mutate(G1pred = predicted, G1stderr = std.error, G1CIH = conf.high, G1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(G1pred, G1stderr, G1CIL, G1CIH, Year, Group)

A1mod <- data.frame(sjPlot::get_model_data(A1, type = "pred", terms = c("Time", "CentiloidS"), mdrt.values = "meansd"))%>% mutate(A1pred = predicted, A1stderr = std.error, A1CIH = conf.high, A1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(A1pred, A1stderr, A1CIL, A1CIH, Year, Group)

# Add all the data together
ab <- merge(F1mod, F2mod, by = c("Group", "Year")); ab <- merge(ab, F3mod, by = c("Group", "Year")); ab <- merge(ab, F4mod, by = c("Group", "Year")); ab <- merge(ab, G1mod, by = c("Group", "Year")); ab <- merge(ab, K1mod, by = c("Group", "Year")) ; ab <- merge(ab, A1mod, by = c("Group", "Year"))

abpred <- ab %>% select(Group, Year, F1pred, F2pred, F3pred, F4pred, G1pred, K1pred, A1pred)%>% gather (Factor, Cognition, F1pred, F2pred, F3pred, F4pred, G1pred, K1pred, A1pred)

abCIH <- ab %>% select(Group, Year, F1CIH, F2CIH, F3CIH, F4CIH, G1CIH, K1CIH, A1CIH)%>% gather (Factor, CIH, F1CIH, F2CIH, F3CIH, F4CIH, G1CIH, K1CIH, A1CIH)%>% mutate(Factor = str_replace_all(Factor, c("F1CIH" = "F1pred", "F2CIH" = "F2pred", "F3CIH" = "F3pred", "F4CIH" = "F4pred", "K1CIH" = "K1pred", "G1CIH" = "G1pred", "A1CIH" = "A1pred")))

abCIL <- ab %>% select(Group, Year, F1CIL, F2CIL, F3CIL, F4CIL, G1CIL, K1CIL, A1CIL)%>% gather (Factor, CIL, F1CIL, F2CIL, F3CIL, F4CIL, G1CIL, K1CIL, A1CIL)%>% mutate(Factor = str_replace_all(Factor, c("F1CIL" = "F1pred", "F2CIL" = "F2pred", "F3CIL" = "F3pred", "F4CIL" = "F4pred", "K1CIL" = "K1pred", "G1CIL" = "G1pred", "A1CIL" = "A1pred")))

abplot <- merge(abpred, abCIH, by = c("Group", "Year", "Factor")); abplot <- merge(abplot, abCIL,  by = c("Group", "Year", "Factor"))
abplot$Factor <- factor(abplot$Factor, levels=c("F1pred","F2pred","F3pred","F4pred","K1pred","A1pred","G1pred"))

#ggplot visualize  
abplot <- abplot %>%
  ggplot(aes(x = Year, y = Cognition,  color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin = CIL, ymax = CIH, color = Group, fill = Group), linetype = "dashed", size = 0.2,alpha = 0.1)+
  scale_color_manual(values = c("#edae49", "#d1495b", "#00798c"), 
                     labels=c("Average", "High", "Low"))+
  scale_fill_manual(values = c("#edae49", "#d1495b", "#00798c"), 
                    labels=c("Average", "High", "Low"))+
  scale_y_continuous(labels=scaleFUN)+
  labs (y="", x = "Time From Imaging Visit (years)", color = "Baseline Amyloid Pathology", fill = "Baseline Amyloid Pathology") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1pred = "Episodic Memory", F2pred = "Semantic Memory", F3pred = "Attention and Processing Speed", F4pred = "Working Memory", K1pred = "Knight ADRC PACC", G1pred = "Global", A1pred = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
abplot

F1 <- lmer(F1 ~ 1 + Time * TauopathyS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_tau, na.action = na.omit, control = lmerCont)
F2 <- lmer(F2 ~ 1 + Time * TauopathyS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_tau, na.action = na.omit, control = lmerCont)
F3 <- lmer(F3 ~ 1 + Time * TauopathyS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_tau, na.action = na.omit, control = lmerCont)
F4 <- lmer(F4 ~ 1 + Time * TauopathyS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_tau, na.action = na.omit, control = lmerCont)
K1 <- lmer(K1 ~ 1 + Time * TauopathyS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_tau, na.action = na.omit, control = lmerCont)
A1 <- lmer(A1 ~ 1 + Time * TauopathyS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_tau, na.action = na.omit, control = lmerCont)
G1 <- lmer(G1 ~ 1 + Time * TauopathyS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_tau, na.action = na.omit, control = lmerCont)

sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T, title = "Baseline Tau and Longitudinal Cognition", file = "Manuscript/baselineTAU_predicts_longcog.doc")
sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T)
performance::r2_nakagawa(F2, tolerance = 1e-08)
performance::r2_nakagawa(F3, tolerance = 1e-08)

sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T)
sjPlot::get_model_data(F1, type = "est", terms = "Time:TauopathyS")[9]
sjPlot::get_model_data(F2, type = "est", terms = "Time:TauopathyS")[9]
sjPlot::get_model_data(F3, type = "est", terms = "Time:TauopathyS")[9]
sjPlot::get_model_data(F4, type = "est", terms = "Time:TauopathyS")[9]
sjPlot::get_model_data(K1, type = "est", terms = "Time:TauopathyS")[9]
sjPlot::get_model_data(G1, type = "est", terms = "Time:TauopathyS")[9]
sjPlot::get_model_data(A1, type = "est", terms = "Time:TauopathyS")[9]

# Collect data 
F1mod <- data.frame(sjPlot::get_model_data(F1, type = "pred", terms = c("Time", "TauopathyS"), mdrt.values = "meansd"))%>% mutate(F1pred = predicted, F1stderr = std.error, F1CIH = conf.high, F1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F1pred, F1stderr, F1CIL, F1CIH, Year, Group)

F2mod <- data.frame(sjPlot::get_model_data(F2, type = "pred", terms = c("Time", "TauopathyS"), mdrt.values = "meansd"))%>% mutate(F2pred = predicted, F2stderr = std.error, F2CIH = conf.high, F2CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F2pred, F2stderr, F2CIL, F2CIH, Year, Group)

F3mod <- data.frame(sjPlot::get_model_data(F3, type = "pred", terms = c("Time", "TauopathyS"), mdrt.values = "meansd"))%>% mutate(F3pred = predicted, F3stderr = std.error, F3CIH = conf.high, F3CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F3pred, F3stderr, F3CIL, F3CIH, Year, Group)

F4mod <- data.frame(sjPlot::get_model_data(F4, type = "pred", terms = c("Time", "TauopathyS"), mdrt.values = "meansd"))%>% mutate(F4pred = predicted, F4stderr = std.error, F4CIH = conf.high, F4CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F4pred, F4stderr, F4CIL, F4CIH, Year, Group)

K1mod <- data.frame(sjPlot::get_model_data(K1, type = "pred", terms = c("Time", "TauopathyS"), mdrt.values = "meansd"))%>% mutate(K1pred = predicted, K1stderr = std.error, K1CIH = conf.high, K1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(K1pred, K1stderr, K1CIL, K1CIH, Year, Group)

G1mod <- data.frame(sjPlot::get_model_data(G1, type = "pred", terms = c("Time", "TauopathyS"), mdrt.values = "meansd"))%>% mutate(G1pred = predicted, G1stderr = std.error, G1CIH = conf.high, G1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(G1pred, G1stderr, G1CIL, G1CIH, Year, Group)

A1mod <- data.frame(sjPlot::get_model_data(A1, type = "pred", terms = c("Time", "TauopathyS"), mdrt.values = "meansd"))%>% mutate(A1pred = predicted, A1stderr = std.error, A1CIH = conf.high, A1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(A1pred, A1stderr, A1CIL, A1CIH, Year, Group)

# Add all the data together
tau <- merge(F1mod, F2mod, by = c("Group", "Year")); tau <- merge(tau, F3mod, by = c("Group", "Year")); tau <- merge(tau, F4mod, by = c("Group", "Year")); tau <- merge(tau, G1mod, by = c("Group", "Year")); tau <- merge(tau, K1mod, by = c("Group", "Year")); tau <- merge(tau, A1mod, by = c("Group", "Year"))

taupred <- tau %>% select(Group, Year, F1pred, F2pred, F3pred, F4pred, G1pred, K1pred, A1pred)%>% gather (Factor, Cognition, F1pred, F2pred, F3pred, F4pred, G1pred, K1pred, A1pred)

tauCIH <- tau %>% select(Group, Year, F1CIH, F2CIH, F3CIH, F4CIH, G1CIH, K1CIH, A1CIH)%>% gather (Factor, CIH, F1CIH, F2CIH, F3CIH, F4CIH, G1CIH, K1CIH, A1CIH)%>% mutate(Factor = str_replace_all(Factor, c("F1CIH" = "F1pred", "F2CIH" = "F2pred", "F3CIH" = "F3pred", "F4CIH" = "F4pred", "K1CIH" = "K1pred", "G1CIH" = "G1pred", "A1CIH" = "A1pred")))

tauCIL <- tau %>% select(Group, Year, F1CIL, F2CIL, F3CIL, F4CIL, G1CIL, K1CIL, A1CIL)%>% gather (Factor, CIL, F1CIL, F2CIL, F3CIL, F4CIL, G1CIL, K1CIL, A1CIL)%>% mutate(Factor = str_replace_all(Factor, c("F1CIL" = "F1pred", "F2CIL" = "F2pred", "F3CIL" = "F3pred", "F4CIL" = "F4pred", "K1CIL" = "K1pred", "G1CIL" = "G1pred", "A1CIL" = "A1pred")))

tauplot <- merge(taupred, tauCIH, by = c("Group", "Year", "Factor")); tauplot <- merge(tauplot, tauCIL,  by = c("Group", "Year", "Factor"))
tauplot$Factor <- factor(tauplot$Factor, levels=c("F1pred","F2pred","F3pred","F4pred","K1pred","A1pred","G1pred"))

#ggplot visualize  
tauplot <- tauplot %>%
  ggplot(aes(x = Year, y = Cognition,  color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin = CIL, ymax = CIH, color = Group, fill = Group), linetype = "dashed", size = 0.2,alpha = 0.1)+
  scale_color_manual(values = c("#edae49", "#d1495b", "#00798c"), 
                     labels=c("Average", "High", "Low"))+
  scale_fill_manual(values = c("#edae49", "#d1495b", "#00798c"), 
                    labels=c("Average", "High", "Low"))+
  scale_y_continuous(labels=scaleFUN)+
  labs (y="Predicted Standardized Cognitive Scores", x = "Time From Imaging Visit (years)", color = "Baseline Tau Pathology", fill = "Baseline Tau Pathology") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1pred = "Episodic Memory", F2pred = "Semantic Memory", F3pred = "Attention and Processing Speed", F4pred = "Working Memory", K1pred = "Knight ADRC PACC", G1pred = "Global", A1pred = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
tauplot

F1 <- lmer(F1 ~ 1 + Time * CortSigS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_mri, na.action = na.omit, control = lmerCont)
F2 <- lmer(F2 ~ 1 + Time * CortSigS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_mri, na.action = na.omit, control = lmerCont)
F3 <- lmer(F3 ~ 1 + Time * CortSigS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_mri, na.action = na.omit, control = lmerCont)
F4 <- lmer(F4 ~ 1 + Time * CortSigS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_mri, na.action = na.omit, control = lmerCont)
K1 <- lmer(K1 ~ 1 + Time * CortSigS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_mri, na.action = na.omit, control = lmerCont)
A1 <- lmer(A1 ~ 1 + Time * CortSigS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_mri, na.action = na.omit, control = lmerCont)
G1 <- lmer(G1 ~ 1 + Time * CortSigS + (Time | ID) +baseline_ageS + EDUCS + sex + apoe, data = cog_mri, na.action = na.omit, control = lmerCont)

sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T, title = "Baseline mri and Longitudinal Cognition", file = "Manuscript/baselinemri_predicts_longcog.doc")
sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T)
sjPlot:::tab_model(F1, F2, F3, F4, K1, G1, A1, show.aic = T)
sjPlot::get_model_data(F1, type = "est", terms = "Time:CortSigS")[9]
sjPlot::get_model_data(F2, type = "est", terms = "Time:CortSigS")[9]
sjPlot::get_model_data(F3, type = "est", terms = "Time:CortSigS")[9]
sjPlot::get_model_data(F4, type = "est", terms = "Time:CortSigS")[9]
sjPlot::get_model_data(K1, type = "est", terms = "Time:CortSigS")[9]
sjPlot::get_model_data(G1, type = "est", terms = "Time:CortSigS")[9]
sjPlot::get_model_data(A1, type = "est", terms = "Time:CortSigS")[9]


# Collect data 
F1mod <- data.frame(sjPlot::get_model_data(F1, type = "pred", terms = c("Time", "CortSigS"), mdrt.values = "meansd"))%>% mutate(F1pred = predicted, F1stderr = std.error, F1CIH = conf.high, F1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F1pred, F1stderr, F1CIL, F1CIH, Year, Group)

F2mod <- data.frame(sjPlot::get_model_data(F2, type = "pred", terms = c("Time", "CortSigS"), mdrt.values = "meansd"))%>% mutate(F2pred = predicted, F2stderr = std.error, F2CIH = conf.high, F2CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F2pred, F2stderr, F2CIL, F2CIH, Year, Group)

F3mod <- data.frame(sjPlot::get_model_data(F3, type = "pred", terms = c("Time", "CortSigS"), mdrt.values = "meansd"))%>% mutate(F3pred = predicted, F3stderr = std.error, F3CIH = conf.high, F3CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F3pred, F3stderr, F3CIL, F3CIH, Year, Group)

F4mod <- data.frame(sjPlot::get_model_data(F4, type = "pred", terms = c("Time", "CortSigS"), mdrt.values = "meansd"))%>% mutate(F4pred = predicted, F4stderr = std.error, F4CIH = conf.high, F4CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(F4pred, F4stderr, F4CIL, F4CIH, Year, Group)

K1mod <- data.frame(sjPlot::get_model_data(K1, type = "pred", terms = c("Time", "CortSigS"), mdrt.values = "meansd"))%>% mutate(K1pred = predicted, K1stderr = std.error, K1CIH = conf.high, K1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(K1pred, K1stderr, K1CIL, K1CIH, Year, Group)

G1mod <- data.frame(sjPlot::get_model_data(G1, type = "pred", terms = c("Time", "CortSigS"), mdrt.values = "meansd"))%>% mutate(G1pred = predicted, G1stderr = std.error, G1CIH = conf.high, G1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(G1pred, G1stderr, G1CIL, G1CIH, Year, Group)

A1mod <- data.frame(sjPlot::get_model_data(A1, type = "pred", terms = c("Time", "CortSigS"), mdrt.values = "meansd"))%>% mutate(A1pred = predicted, A1stderr = std.error, A1CIH = conf.high, A1CIL = conf.low, Year = x, group = as.numeric(group), Group = case_when(group == 3 ~ "High", group == 2 ~ "Average", group == 1 ~ "Low"))%>% select(A1pred, A1stderr, A1CIL, A1CIH, Year, Group)

# Add all the data together
mri <- merge(F1mod, F2mod, by = c("Group", "Year")); mri <- merge(mri, F3mod, by = c("Group", "Year")); mri <- merge(mri, F4mod, by = c("Group", "Year")); mri <- merge(mri, G1mod, by = c("Group", "Year")); mri <- merge(mri, K1mod, by = c("Group", "Year")); mri <- merge(mri, A1mod, by = c("Group", "Year"))

mripred <- mri %>% select(Group, Year, F1pred, F2pred, F3pred, F4pred, G1pred, K1pred, A1pred)%>% gather (Factor, Cognition, F1pred, F2pred, F3pred, F4pred, G1pred, K1pred, A1pred)

mriCIH <- mri %>% select(Group, Year, F1CIH, F2CIH, F3CIH, F4CIH, G1CIH, K1CIH, A1CIH)%>% gather (Factor, CIH, F1CIH, F2CIH, F3CIH, F4CIH, G1CIH, K1CIH, A1CIH)%>% mutate(Factor = str_replace_all(Factor, c("F1CIH" = "F1pred", "F2CIH" = "F2pred", "F3CIH" = "F3pred", "F4CIH" = "F4pred", "K1CIH" = "K1pred", "G1CIH" = "G1pred", "A1CIH" = "A1pred")))

mriCIL <- mri %>% select(Group, Year, F1CIL, F2CIL, F3CIL, F4CIL, G1CIL, K1CIL, A1CIL)%>% gather (Factor, CIL, F1CIL, F2CIL, F3CIL, F4CIL, G1CIL, K1CIL, A1CIL)%>% mutate(Factor = str_replace_all(Factor, c("F1CIL" = "F1pred", "F2CIL" = "F2pred", "F3CIL" = "F3pred", "F4CIL" = "F4pred", "K1CIL" = "K1pred", "G1CIL" = "G1pred", "A1CIL" = "A1pred")))

mriplot <- merge(mripred, mriCIH, by = c("Group", "Year", "Factor")); mriplot <- merge(mriplot, mriCIL,  by = c("Group", "Year", "Factor"))
mriplot$Factor <- factor(mriplot$Factor, levels=c("F1pred","F2pred","F3pred","F4pred","K1pred","A1pred","G1pred"))

#ggplot visualize  
mriplot <- mriplot %>%
  ggplot(aes(x = Year, y = Cognition,  color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin = CIL, ymax = CIH, color = Group, fill = Group), linetype = "dashed", size = 0.2,alpha = 0.1)+
  scale_color_manual(values = c("#edae49", "#d1495b", "#00798c"), 
                     labels=c("Average", "High", "Low"))+
  scale_fill_manual(values = c("#edae49", "#d1495b", "#00798c"), 
                    labels=c("Average", "High", "Low"))+
  scale_y_continuous(labels=scaleFUN)+
  labs (y="", x = "Time From Imaging Visit (years)", color = "Baseline Cortical Atrophy", fill = "Baseline Cortical Atrophy") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1pred = "Episodic Memory", F2pred = "Semantic Memory", F3pred = "Attention and Processing Speed", F4pred = "Working Memory", K1pred = "Knight ADRC PACC", G1pred = "Global", A1pred = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
mriplot

bio_fig <- ggarrange(abplot + rremove("ylab"), tauplot +rremove("ylab"), mriplot +rremove("ylab"), ncol = 1)

annotate_figure(bio_fig, left = text_grob("Predicted Standardized Cognitive Score", rot = 90))

ggsave("Manuscript/all_biomarkers_cog.png", width=10,height =8)  

#keep(full_data, scaleFUN, ab, tau, mri, sure = TRUE)
#### FIGURE 8: SURVIVAL CURVES ####
pacman::p_load(survival, sjPlot, survminer, longpower, hrbrthemes)
full_data$na_count <- apply(full_data [, c(12:15)], 1, function(x) sum(is.na(x)))

survival <- full_data%>% filter(na_count < 1 & BLcdr_bin == 0)%>% arrange(psy_Date)%>% group_by(ID)%>% mutate(count = "1", visit = cumsum(count), maxvisit = max(visit), min_Date = psy_Date[1], maxDate = max(psy_Date), maxTime = max(Time), BLAge =Age[1], BLcdr_bin = cdr_bin[1], maxcdr_bin = max(cdr_bin), con = case_when(maxcdr_bin == 0 ~ 0, maxcdr_bin != 0 ~ 1))%>% ungroup()%>% filter(maxvisit >= 4 & maxTime >= 4*365.25)%>% mutate(BLAgeS = as.numeric(scale(BLAge)), EDUCS = as.numeric(scale(EDUC, scale = FALSE)))
  
survival1 <- survival%>% filter(con == 0)%>% group_by(ID)%>% arrange(psy_Date)%>% mutate(conDate = max(psy_Date), conTime = as.numeric(conDate - min_Date))%>% ungroup()%>% distinct(ID, conDate, conTime)

survival2 <- survival%>% filter(maxcdr_bin > 0 & cdr_bin > 0)%>% group_by(ID)%>% arrange(psy_Date)%>% mutate(conDate = psy_Date[1], conTime = as.numeric(conDate - min_Date))%>% ungroup()%>% distinct(ID, conDate, conTime)
survival3 <- rbind(survival1, survival2)
survival <- merge(survival, survival3, by = "ID")

#IDs to run models for 
MapIDs <- unique(survival$ID)

# Get cog slopes
F1_coefs <- rep(NA, length(MapIDs))
j = 1
for (i in MapIDs) {
  # Slice by individual participant 
  tmp = survival[survival$ID == i,] 
  # Run regression
  results_F1 = lm(F1 ~ Time, data = tmp)  
  # Take coefficient 
  tmp_coef_cog_F1 = coef(results_F1) 
  # Move to correct place 
  F1_coefs[j] = tmp_coef_cog_F1[2] 
  
  j = j + 1
}
F1_slope <- as.data.frame(F1_coefs)

F2_coefs <- rep(NA, length(MapIDs))
j = 1
for (i in MapIDs) {
  # Slice by individual participant 
  tmp = survival[survival$ID == i,] 
  # Run regression
  results_F2 = lm(F2 ~ Time, data = tmp)  
  # Take coefficient 
  tmp_coef_cog_F2 = coef(results_F2) 
  # Move to correct place 
  F2_coefs[j] = tmp_coef_cog_F2[2] 
  
  j = j + 1
}
F2_slope <- as.data.frame(F2_coefs)

F3_coefs <- rep(NA, length(MapIDs))
j = 1
for (i in MapIDs) {
  # Slice by individual participant 
  tmp = survival[survival$ID == i,] 
  # Run regression
  results_F3 = lm(F3 ~ Time, data = tmp)  
  # Take coefficient 
  tmp_coef_cog_F3 = coef(results_F3) 
  # Move to correct place 
  F3_coefs[j] = tmp_coef_cog_F3[2] 
  
  j = j + 1
}
F3_slope <- as.data.frame(F3_coefs)

F4_coefs <- rep(NA, length(MapIDs))
j = 1
for (i in MapIDs) {
  # Slice by individual participant 
  tmp = survival[survival$ID == i,] 
  # Run regression
  results_F4 = lm(F4 ~ Time, data = tmp)  
  # Take coefficient 
  tmp_coef_cog_F4 = coef(results_F4) 
  # Move to correct place 
  F4_coefs[j] = tmp_coef_cog_F4[2] 
  
  j = j + 1
}
F4_slope <- as.data.frame(F4_coefs)

K1_coefs <- rep(NA, length(MapIDs))
j = 1
for (i in MapIDs) {
  # Slice by individual participant 
  tmp = survival[survival$ID == i,] 
  # Run regression
  results_K1 = lm(K1 ~ Time, data = tmp)  
  # Take coefficient 
  tmp_coef_cog_K1 = coef(results_K1) 
  # Move to correct place 
  K1_coefs[j] = tmp_coef_cog_K1[2] 
  
  j = j + 1
}
K1_slope <- as.data.frame(K1_coefs)

G1_coefs <- rep(NA, length(MapIDs))
j = 1
for (i in MapIDs) {
  # Slice by individual participant 
  tmp = survival[survival$ID == i,] 
  # Run regression
  results_G1 = lm(G1 ~ Time, data = tmp)  
  # Take coefficient 
  tmp_coef_cog_G1 = coef(results_G1) 
  # Move to correct place 
  G1_coefs[j] = tmp_coef_cog_G1[2] 
  
  j = j + 1
}
G1_slope <- as.data.frame(G1_coefs)

A1_coefs <- rep(NA, length(MapIDs))
j = 1
for (i in MapIDs) {
  # Slice by individual participant 
  tmp = survival[survival$ID == i,] 
  # Run regression
  results_A1 = lm(A1 ~ Time, data = tmp)  
  # Take coefficient 
  tmp_coef_cog_A1 = coef(results_A1) 
  # Move to correct place 
  A1_coefs[j] = tmp_coef_cog_A1[2] 
  
  j = j + 1
}
A1_slope <- as.data.frame(A1_coefs)

Results_cog <- cbind(MapIDs, F1_slope, F2_slope, F3_slope, F4_slope, K1_slope, G1_slope, A1_slope)%>% drop_na()%>% mutate(ID = MapIDs)%>% select(-MapIDs)

survival <- survival%>% select(ID, psy_Date, cdr_bin, BLAgeS, apoe, EDUCS, sex, race, BLcdr_bin, maxcdr_bin, maxvisit, min_Date, Time, con, conDate, conTime)%>% merge(., Results_cog, by = "ID")%>% group_by(ID)%>% arrange(psy_Date)%>% slice_min(psy_Date)%>% ungroup()

# calculate the median for each factor 
med_F1 <- median(survival$F1_coefs)
med_F2 <- median(survival$F2_coefs)
med_F3 <- median(survival$F3_coefs)
med_F4 <- median(survival$F4_coefs)
med_K1 <- median(survival$K1_coefs)
med_G1 <- median(survival$G1_coefs)
med_A1 <- median(survival$A1_coefs)

#then bin factors as high or lw based on median vals 
surv_cog <- survival%>% mutate(F1cut = med_F1, F1bin = case_when(F1_coefs < F1cut ~ "Low", F1_coefs >= F1cut ~ "High"))
surv_cog <- surv_cog%>% mutate(F2cut = med_F2, F2bin = case_when(F2_coefs < F2cut ~ "Low", F2_coefs >= F2cut ~ "High"))
surv_cog <- surv_cog%>% mutate(F3cut = med_F3, F3bin = case_when(F3_coefs < F3cut ~ "Low", F3_coefs >= F3cut ~ "High"))
surv_cog <- surv_cog%>% mutate(F4cut = med_F4, F4bin = case_when(F4_coefs < F4cut ~ "Low", F4_coefs >= F4cut ~ "High"))
surv_cog <- surv_cog%>% mutate(K1cut = med_K1, K1bin = case_when(K1_coefs < K1cut ~ "Low", K1_coefs >= K1cut ~ "High"))
surv_cog <- surv_cog%>% mutate(G1cut = med_G1, G1bin = case_when(G1_coefs < G1cut ~ "Low", G1_coefs >= G1cut ~ "High"))
surv_cog <- surv_cog%>% mutate(A1cut = med_A1, A1bin = case_when(A1_coefs < A1cut ~ "Low", A1_coefs >= A1cut ~ "High"))%>% mutate(con = as.numeric(con))

cog_Surv <- Surv(time = surv_cog$conTime, event = surv_cog$con)

cog_Survfit_F1 <- survfit(cog_Surv ~ surv_cog$F1bin, data = surv_cog)
cog_Survfit_F2 <- survfit(cog_Surv ~ surv_cog$F2bin, data = surv_cog)
cog_Survfit_F3 <- survfit(cog_Surv ~ surv_cog$F3bin, data = surv_cog)
cog_Survfit_F4 <- survfit(cog_Surv ~ surv_cog$F4bin, data = surv_cog)
cog_Survfit_K1 <- survfit(cog_Surv ~ surv_cog$K1bin, data = surv_cog)
cog_Survfit_G1 <- survfit(cog_Surv ~ surv_cog$G1bin, data = surv_cog)
cog_Survfit_A1 <- survfit(cog_Surv ~ surv_cog$A1bin, data = surv_cog)

F1 <- coxph(Surv(conTime, con) ~ surv_cog$F1bin + BLAgeS + sex + EDUCS + apoe, data= surv_cog)
F2 <- coxph(Surv(conTime, con) ~ surv_cog$F2bin + BLAgeS + sex + EDUCS + apoe, data= surv_cog)
F3 <- coxph(Surv(conTime, con) ~ surv_cog$F3bin + BLAgeS + sex + EDUCS + apoe, data= surv_cog)
F4 <- coxph(Surv(conTime, con) ~ surv_cog$F4bin + BLAgeS + sex + EDUCS + apoe, data= surv_cog)
K1 <- coxph(Surv(conTime, con) ~ surv_cog$K1bin + BLAgeS + sex + EDUCS + apoe, data= surv_cog)
G1 <- coxph(Surv(conTime, con) ~ surv_cog$G1bin + BLAgeS + sex + EDUCS + apoe, data= surv_cog)
A1 <- coxph(Surv(conTime, con) ~ surv_cog$A1bin + BLAgeS + sex + EDUCS + apoe, data= surv_cog)

write.csv(surv_cog, 'Manuscript/demo_survival_analysis.csv')

summary(F1) 
summary(F2) 
summary(F3) 
summary(F4) 
summary(K1) 
summary(G1) 
summary(A1)

splots <- list()
F1surv <- data.frame(matrix(ncol = 4, nrow = sum(data.frame(cog_Survfit_F1$strata))))
F1surv[,1] <- data.frame(cog_Survfit_F1$time)
F1surv[,2] <- data.frame(cog_Survfit_F1$surv)
F1surv[,3] <- data.frame(cog_Survfit_F1$lower)
F1surv[,4] <- data.frame(cog_Survfit_F1$upper)
colnames(F1surv) <- c("Time", "Surv", "Low", "High")
F1surv <- F1surv %>%
  mutate(Group = c(rep("High", times = cog_Survfit_F1$strata[1]), rep("Low", times = cog_Survfit_F1$strata[2])), Factor = "F1")

F1Hint <- F1surv %>% filter(Group == "High")
F1Hint <- data.frame(F1Hint[which.min(abs(0.5-F1Hint$Surv)),])$Time
F1Lint <- F1surv %>% filter(Group == "Low")
F1Lint <- data.frame(F1Lint[which.min(abs(0.5-F1Lint$Surv)),])$Time

F2surv <- data.frame(matrix(ncol = 4, nrow = sum(data.frame(cog_Survfit_F2$strata))))
F2surv[,1] <- data.frame(cog_Survfit_F2$time)
F2surv[,2] <- data.frame(cog_Survfit_F2$surv)
F2surv[,3] <- data.frame(cog_Survfit_F2$lower)
F2surv[,4] <- data.frame(cog_Survfit_F2$upper)
colnames(F2surv) <- c("Time", "Surv", "Low", "High")
F2surv <- F2surv %>%
  mutate(Group = c(rep("High", times = cog_Survfit_F2$strata[1]), rep("Low", times = cog_Survfit_F2$strata[2])), Factor = "F2")

F2Hint <- F2surv %>% filter(Group == "High")
F2Hint <- data.frame(F2Hint[which.min(abs(0.5-F2Hint$Surv)),])$Time
F2Lint <- F2surv %>% filter(Group == "Low")
F2Lint <- data.frame(F2Lint[which.min(abs(0.5-F2Lint$Surv)),])$Time

F3surv <- data.frame(matrix(ncol = 4, nrow = sum(data.frame(cog_Survfit_F3$strata))))
F3surv[,1] <- data.frame(cog_Survfit_F3$time)
F3surv[,2] <- data.frame(cog_Survfit_F3$surv)
F3surv[,3] <- data.frame(cog_Survfit_F3$lower)
F3surv[,4] <- data.frame(cog_Survfit_F3$upper)
colnames(F3surv) <- c("Time", "Surv", "Low", "High")
F3surv <- F3surv %>%
  mutate(Group = c(rep("High", times = cog_Survfit_F3$strata[1]), rep("Low", times = cog_Survfit_F3$strata[2])), Factor = "F3")

F3Hint <- F3surv %>% filter(Group == "High")
F3Hint <- data.frame(F3Hint[which.min(abs(0.5-F3Hint$Surv)),])$Time
F3Lint <- F3surv %>% filter(Group == "Low")
F3Lint <- data.frame(F3Lint[which.min(abs(0.5-F3Lint$Surv)),])$Time

F4surv <- data.frame(matrix(ncol = 4, nrow = sum(data.frame(cog_Survfit_F4$strata))))
F4surv[,1] <- data.frame(cog_Survfit_F4$time)
F4surv[,2] <- data.frame(cog_Survfit_F4$surv)
F4surv[,3] <- data.frame(cog_Survfit_F4$lower)
F4surv[,4] <- data.frame(cog_Survfit_F4$upper)
colnames(F4surv) <- c("Time", "Surv", "Low", "High")
F4surv <- F4surv %>%
  mutate(Group = c(rep("High", times = cog_Survfit_F4$strata[1]), rep("Low", times = cog_Survfit_F4$strata[2])), Factor = "F4")

F4Hint <- F4surv %>% filter(Group == "High")
F4Hint <- data.frame(F4Hint[which.min(abs(0.5-F4Hint$Surv)),])$Time
F4Lint <- F4surv %>% filter(Group == "Low")
F4Lint <- data.frame(F4Lint[which.min(abs(0.5-F4Lint$Surv)),])$Time

K1surv <- data.frame(matrix(ncol = 4, nrow = sum(data.frame(cog_Survfit_K1$strata))))
K1surv[,1] <- data.frame(cog_Survfit_K1$time)
K1surv[,2] <- data.frame(cog_Survfit_K1$surv)
K1surv[,3] <- data.frame(cog_Survfit_K1$lower)
K1surv[,4] <- data.frame(cog_Survfit_K1$upper)
colnames(K1surv) <- c("Time", "Surv", "Low", "High")
K1surv <- K1surv %>%
  mutate(Group = c(rep("High", times = cog_Survfit_K1$strata[1]), rep("Low", times = cog_Survfit_K1$strata[2])), Factor = "K1")

K1Hint <- K1surv %>% filter(Group == "High")
K1Hint <- data.frame(K1Hint[which.min(abs(0.5-K1Hint$Surv)),])$Time
K1Lint <- K1surv %>% filter(Group == "Low")
K1Lint <- data.frame(K1Lint[which.min(abs(0.5-K1Lint$Surv)),])$Time

G1surv <- data.frame(matrix(ncol = 4, nrow = sum(data.frame(cog_Survfit_G1$strata))))
G1surv[,1] <- data.frame(cog_Survfit_G1$time)
G1surv[,2] <- data.frame(cog_Survfit_G1$surv)
G1surv[,3] <- data.frame(cog_Survfit_G1$lower)
G1surv[,4] <- data.frame(cog_Survfit_G1$upper)
colnames(G1surv) <- c("Time", "Surv", "Low", "High")
G1surv <- G1surv %>%
  mutate(Group = c(rep("High", times = cog_Survfit_G1$strata[1]), rep("Low", times = cog_Survfit_G1$strata[2])), Factor = "G1")

G1Hint <- G1surv %>% filter(Group == "High")
G1Hint <- data.frame(G1Hint[which.min(abs(0.5-G1Hint$Surv)),])$Time
G1Lint <- G1surv %>% filter(Group == "Low")
G1Lint <- data.frame(G1Lint[which.min(abs(0.5-G1Lint$Surv)),])$Time

A1surv <- data.frame(matrix(ncol = 4, nrow = sum(data.frame(cog_Survfit_A1$strata))))
A1surv[,1] <- data.frame(cog_Survfit_A1$time)
A1surv[,2] <- data.frame(cog_Survfit_A1$surv)
A1surv[,3] <- data.frame(cog_Survfit_A1$lower)
A1surv[,4] <- data.frame(cog_Survfit_A1$upper)
colnames(A1surv) <- c("Time", "Surv", "Low", "High")
A1surv <- A1surv %>%
  mutate(Group = c(rep("High", times = cog_Survfit_A1$strata[1]), rep("Low", times = cog_Survfit_A1$strata[2])), Factor = "A1")

A1Hint <- A1surv %>% filter(Group == "High")
A1Hint <- data.frame(A1Hint[which.min(abs(0.5-A1Hint$Surv)),])$Time
A1Lint <- A1surv %>% filter(Group == "Low")
A1Lint <- data.frame(A1Lint[which.min(abs(0.5-A1Lint$Surv)),])$Time

Intercepts <- data.frame(rbind(F1Hint, F1Lint, F2Hint, F2Lint, F3Hint, F3Lint, F4Hint, F4Lint, K1Hint, K1Lint, A1Hint, A1Lint, G1Hint, G1Lint))%>%
  mutate(Factor = c("F1","F1", "F2", "F2", "F3", "F3", "F4", "F4", "K1", "K1", "A1", "A1", "G1", "G1"), Group = c("High","Low", "High","Low", "High","Low", "High","Low", "High","Low", "High","Low", "High","Low"))
colnames(Intercepts) <- c("Int", "Factor", "Group")

# Merge model data 
survival <- rbind(F1surv, F2surv, F3surv, F4surv, G1surv, K1surv, A1surv)
survival$Factor <- factor(survival$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))
Intercepts$Factor <- factor(Intercepts$Factor, levels=c("F1","F2","F3","F4","K1","A1","G1"))

# Plot using ggplot
survplot <- survival %>%
  ggplot(aes(x = Time, y = Surv,  color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin = Low, ymax = High, color = Group, fill = Group), linetype = "dashed", size = 0.2,alpha = 0.1)+
  scale_color_manual(values = c("black", "purple"), labels=c("Low", "High"))+
  scale_fill_manual(values = c("black", "purple"),  labels=c("Low", "High"))+
  labs (y="Survival Probability", x = "Time From Baseline Visit (days)", color = "Rate of Cognitive\nDecline", fill = "Rate of Cognitive\nDecline") +
  scale_y_continuous(labels=scaleFUN)+
  geom_vline(data= Intercepts, aes(xintercept = Int, color = Group), linetype="dotted", size=0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), axis.text = element_text(size = 10), axis.title = element_text (size = 12))+
  facet_wrap(Factor ~ ., ncol = 4, labeller = as_labeller(c(F1 = "Episodic Memory", F2 = "Semantic Memory", F3 = "Attention and Processing Speed", F4 = "Working Memory", K1 = "Knight ADRC PACC", G1 = "Global", A1 = "ADCS PACC")))+
  theme(strip.placement = "outside")+
  geom_hline(yintercept=0.5, linetype="dotted", color = "darkgrey", size=0.5)+
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
survplot
ggsave("Manuscript/survival.png", width=10,height =4)  

#keep(full_data, scaleFUN, ab, tau, mri, sure = TRUE)
#### FIGURE 9: POWER ANALYSES ####
pacman::p_load(lmerTest, ggbreak)
lmerCont <- lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun=10000))
power <- full_data%>% filter(Year > 2015)
power$F_na_count <- apply(power [, c("F1", "F2", "F3", "F4")], 1, function(x) sum(is.na(x)))

power <- power%>% 
  filter(F_na_count < 1)%>% 
  group_by(ID)%>% 
  arrange(psy_Date)%>% 
  mutate(count = "1", visit = cumsum(count), BLabpos = abpos[1], BLage = Age[1], minDate = psy_Date[1], Time = (as.numeric(psy_Date-minDate))/365.25, maxvisit = max(visit), BLcdr_bin = cdr_bin[1])%>% 
  ungroup()%>% mutate(EDUCS = as.numeric(scale(EDUC, scale = FALSE)), BLAgeS = as.numeric(scale(BLage)))%>% 
  filter(maxvisit > 1)%>% filter(BLcdr_bin == 0)%>% 
  select(ID, F1:A1, sex, apoe, race, BLcdr_bin, BLabpos, BLAgeS, EDUCS, Time, visit)

F1 <- lmer(F1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = power, na.action = na.omit, control = lmerCont)
summary(F1)
F1b <- lmmpower(F1, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

F2 <- lmer(F2 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = power, na.action = na.omit, control = lmerCont)
summary(F2)
F2b <- lmmpower(F2, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

F3 <- lmer(F3 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = power, na.action = na.omit, control = lmerCont)
summary(F3)
F3b <- lmmpower(F3, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

F4 <- lmer(F4 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = power, na.action = na.omit, control = lmerCont)
summary(F4)
F4b <- lmmpower(F4, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

K1 <- lmer(K1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = power, na.action = na.omit, control = lmerCont)
summary(K1)
K1b <- lmmpower(K1, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

Gb <- lmer(G1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = power, na.action = na.omit, control = lmerCont)
summary(Gb)
Gbb <- lmmpower(Gb, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

Ad <- lmer(A1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = power, na.action = na.omit, control = lmerCont)
summary(Ad)
Adb <- lmmpower(Ad, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

all1 <- rbind(F1b$n, F2b$n, F3b$n, F4b$n, K1b$n, Gbb$n, Adb$n)
all2 <- rbind(F1b$n.CI, F2b$n.CI, F3b$n.CI, F4b$n.CI, K1b$n.CI, Gbb$n.CI, Adb$n.CI)
all <- cbind(all1, all2)
all <- data.frame(all)%>% mutate(n = X1, Lower = X5, Upper = X3)%>% dplyr::select(n, Lower, Upper)%>%
  mutate(Label = c("F1", "F2", "F3", "F4", "K1", "G1", "A1"), Group = "All")

powerx <- power %>%
  filter(BLabpos == 1)

F1 <- lmer(F1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = powerx, na.action = na.omit, control = lmerCont)
F1b <- lmmpower(F1, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

F2 <- lmer(F2 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = powerx, na.action = na.omit, control = lmerCont)
F2b <- lmmpower(F2, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

F3 <- lmer(F3 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = powerx, na.action = na.omit, control = lmerCont)
F3b <- lmmpower(F3, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

F4 <- lmer(F4 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = powerx, na.action = na.omit, control = lmerCont)
F4b <- lmmpower(F4, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

K1 <- lmer(K1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = powerx, na.action = na.omit, control = lmerCont)
K1b <- lmmpower(K1, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

Gb <- lmer(G1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = powerx, na.action = na.omit, control = lmerCont)
Gbb <- lmmpower(Gb, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

Ad <- lmer(A1 ~ 1 + Time + (Time | ID) + BLAgeS + sex + apoe  + EDUCS, data = powerx, na.action = na.omit, control = lmerCont)
Adb <- lmmpower(Ad, pct.change = 0.50, t = seq (0, 4, 1), power = 0.80)

abpos <- rbind(F1b$n, F2b$n, F3b$n, F4b$n, K1b$n, Gbb$n, Adb$n)

abpos1 <- rbind(F1b$n, F2b$n, F3b$n, F4b$n, K1b$n, Gbb$n, Adb$n)
abpos2 <- rbind(F1b$n.CI, F2b$n.CI, F3b$n.CI, F4b$n.CI, K1b$n.CI, Gbb$n.CI, Adb$n.CI)
abpos <- cbind(abpos1, abpos2)
abpos <- data.frame(abpos)%>% mutate(n = X1, Lower = X5, Upper = X3)%>% dplyr::select(n, Lower, Upper)%>% mutate(Label = c("F1", "F2", "F3", "F4", "K1", "G1", "A1"), Group = "AB+")

all
abpos

#merge two dataframes 
power_forest <- data.frame(rbind(all, abpos))%>%
  mutate(Label = factor(Label, levels = c("G1", "A1", "K1", "F4", "F3", "F2", "F1"))) %>%
  mutate(Group = factor (Group, levels = c("All", "AB+")))%>%
  mutate(Label = str_replace_all(Label, c("G1" = "Global", "A1" = "ADCS PACC", "K1" = "Knight ADRC \nPACC", "F4" = "Working \nMemory", "F3" = "Attention and \nProcessing Speed", "F2" = "Semantic \nMemory", "F1" = "Episodic \nMemory")))%>%
  mutate(Group = str_replace_all(Group, c("All" = "Full Sample", "AB+" = "Amyloid \nPositive")))

power_forest%>% 
  mutate(across(Label, ~factor(., levels=c("Global","ADCS PACC","Knight ADRC \nPACC","Working \nMemory","Attention and \nProcessing Speed","Semantic \nMemory","Episodic \nMemory"))))%>%
  ggplot(., aes(x=Label, y=n, ymin=Lower, ymax=Upper, col = Group, fill = Group)) +
  geom_pointrange(position = position_dodge(width = 0.4))+ 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Estimated Sample Size (95% CI)") +
  scale_y_break(c(1200, 2500), scale = 0.4, space = 0.5)+
  scale_color_manual(values = c("black", "purple"))+
  scale_fill_manual(values = c("black", "purple"))+
  theme_apa()+  # use a white background 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave("Manuscript/Power_CI.png", width=10,height =4)  

#keep (full_data, sure = TRUE)
#### DEMOGRAPHICS: TABLE  ####
#remotes::install_github("ddsjoberg/gtsummary")
pacman::p_load(gtsummary,lsr)

# summary table
table1 <- full_data%>% group_by(ID)%>% 
  arrange(psy_Date)%>% slice_min(psy_Date)%>% mutate(apoe = as.factor(apoe), race = as.factor(race), sex = as.factor(sex), abpos = as.factor(abpos), taupos = as.factor(taupos))

table_one <- table1%>% tbl_summary(
  include = c(Age, apoe, EDUC, race, sex, abpos, taupos, eqMMSE, F1, F2, F3, F4, K1, A1, G1), by = cdr_bin, missing_text = "Missing", type = all_continuous() ~ "continuous2", statistic = all_continuous() ~ c("{mean} ({sd})", "{median} ({p0}-{p100})"))%>% add_n()%>% add_difference(test = list(all_continuous() ~ "cohens_d"))%>% add_p(test = list(all_continuous() ~ c("t.test") , all_categorical() ~ "kruskal.test"),pvalue_fun = ~style_pvalue(.x,digits =2))%>% modify_fmt_fun(statistic ~ style_sigfig)%>% bold_labels()
table_one%>% as_gt()%>%
  gt::gtsave(filename = "Manuscript/Table_One_Demographics.png")

# individuals statistical tests to get the effect sizes 
describe(table1)
describeBy(table1$Age, table1$cdr_bin)
mod <- t.test(Age ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(Age ~ cdr_bin, data = table1)

describeBy(table1$EDUC, table1$cdr_bin)
mod <- t.test(EDUC ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(EDUC ~ cdr_bin, data = table1)

test <- kruskal.test(sex ~ cdr_bin, data = table1); test
tbl = table(table1$sex, table1$cdr_bin); tbl
sqrt(test$statistic /nrow(table1)); test$p.value

test <- kruskal.test(apoe ~ cdr_bin, data = table1); test
tbl = table(table1$apoe, table1$cdr_bin); tbl
sqrt(test$statistic /nrow(table1)); test$p.value

test <- kruskal.test(race ~ cdr_bin, data = table1); test
tbl = table(table1$race, table1$cdr_bin); tbl
sqrt(test$statistic /nrow(table1)); test$p.value

describeBy(table1$eqMMSE, table1$cdr_bin)
mod <- t.test(eqMMSE ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(eqMMSE ~ cdr_bin, data = table1)

describeBy(table1$F1, table1$cdr_bin)
mod <- t.test(F1 ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(F1 ~ cdr_bin, data = table1)

describeBy(table1$F2, table1$cdr_bin)
mod <- t.test(F2 ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(F2 ~ cdr_bin, data = table1)

describeBy(table1$F3, table1$cdr_bin)
mod <- t.test(F3 ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(F3 ~ cdr_bin, data = table1)

describeBy(table1$F4, table1$cdr_bin)
mod <- t.test(F4 ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(F4 ~ cdr_bin, data = table1)

describeBy(table1$K1, table1$cdr_bin)
mod <- t.test(K1 ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(K1 ~ cdr_bin, data = table1)

describeBy(table1$A1, table1$cdr_bin)
mod <- t.test(A1 ~ cdr_bin, data = table1)
mod; mod$p.value  geom_density(alpha = 0.6)
cohensD(A1 ~ cdr_bin, data = table1)

describeBy(table1$G1, table1$cdr_bin)
mod <- t.test(G1 ~ cdr_bin, data = table1)
mod; mod$p.value
cohensD(G1 ~ cdr_bin, data = table1)

table1%>% mutate(cdr_bin = as.factor(cdr_bin))%>%
  ggplot(., aes(x=cdr, color = cdr_bin, fill = cdr_bin, group = cdr_bin))+ 
  geom_histogram(bins = 5, alpha = 0.5)+
  scale_color_manual(values = c("#EBB261","#5A4A6F"), labels=c('Cognitively Normal', 'Cognitively Impaired'))+
  scale_fill_manual(values = c("#EBB261","#5A4A6F"), labels=c( 'Cognitively Normal', 'Cognitively Impaired'))+
  scale_y_continuous(labels = scaleFUN)+
  labs (y="Frequency at Baseline", x = "Cognitive Dementia Rating Score", color = "Cognitive Status", fill = "Cognitive Status") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(linewidth = 0.5, colour = "black"), 
        axis.text = element_text(size = 10), axis.title = element_text (size = 12))
ggsave("Manuscript/CDR_Histograms.png", width=6,height =4) 


