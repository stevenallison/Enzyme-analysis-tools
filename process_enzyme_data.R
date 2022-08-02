library(openxlsx)
library(tidyverse)
library(lubridate)

# To run this script, download the extracellular_enzyme_data folder from the WTFProject Drive

# Make a list of the file names
enzyme_files <- list.files(path = "extracellular_enzyme_data/raw_data",pattern = ".txt",full.names = TRUE, recursive = TRUE)

# Read in plate data
enzymes <- data.frame()
for(x in 1:length(enzyme_files)) {
  file <- as.data.frame(read_table(enzyme_files[x],col_names=FALSE))
  names(file) <- c("Plate","Well","Read")
  enzymes <- rbind(enzymes,file)
}

# Parse out plate information
enzymes2 <- enzymes %>%
  separate(Plate, c("assay_date","PlateType","SampleID","extraID"), convert=T, sep="_") %>%
  unite(SampleID,extraID, col="SampleID",sep="_",na.rm=TRUE) %>%
  mutate(assay_date = ymd(assay_date),
         Read = as.numeric(Read))

# Correct naming issues, select correct samples, remove samples that need to be redone
# Some samples have multiple reads, the correct read is marked by "_reread"
enzyme_redo <- enzymes2 %>%
  separate(SampleID,c("SampleID","extraID"),sep="_") %>%
  filter(extraID == "reread") %>%
  mutate(is_redo = TRUE) %>%
  select(assay_date,PlateType,SampleID,is_redo) %>%
  unique()

enzymes3 <- enzymes2 %>%
  left_join(enzyme_redo,by=c("assay_date","PlateType","SampleID")) %>%
  mutate(is_redo = replace(is_redo,is.na(is_redo),FALSE)) %>%
  filter(is_redo == FALSE) %>%
  mutate(SampleID = gsub("_reread","",SampleID)) %>%
  select(-is_redo)

# Remove other samples with errors
enzymes4 <- enzymes3 %>%
  separate(SampleID,c("SampleID","extraID"),sep="_") %>%
  filter(is.na(extraID)) %>%
  filter(!SampleID%in%c("test","a")) %>%
  mutate(SampleID = as.numeric(SampleID)) %>% # fixes blank IDs for buffer plates caused by unite
  select(-extraID)

# Merge plate layout info
Sample.Info <-read.xlsx("extracellular_enzyme_data/EnzymePlateInfo.xlsx",1)
enzymes5 <- left_join(enzymes4,Sample.Info, by=c("PlateType","Well"))

# Clean clear plate readings by deleting samples with reads that could interfere with the MM fit
enzymes6 <- enzymes5 %>%
  separate(Well,into = c("Row", "Col"),sep ="(?<=[A-Za-z])(?=[0-9])",remove = FALSE) %>%
  mutate(Remove = "No",
         Col = as.numeric(Col))

enzymes_P <- enzymes6 %>%
  filter(PlateType == "P") %>%
  arrange(SampleID,Col)

# Input desired threshold here
# 1.65 is currently the lowest threshold we can use for clear plates without removing valid points
threshold <- 1.65
# Counter helps keep track of when we finish reading every plate
well_counter <- 0
for (row in 1:nrow(enzymes_P)) {
  well_counter <- well_counter + 1
  # Corrects for row A, the top row of the plate
  if (enzymes_P[row,]$Row == "A") {
    if (enzymes_P[row, 'Read'] > threshold*enzymes_P[row+1, 'Read']) {
      enzymes_P[row, 'Remove'] = 'Yes'
    }
  }
  # Corrects for row H, the bottom row of the plate
  else if (enzymes_P[row,]$Row == "H") {
    if (enzymes_P[row, 'Read'] > threshold*enzymes_P[row-1, 'Read']) {
      enzymes_P[row, 'Remove'] = 'Yes'
    }
  }
  # Corrects for rows in middle
  else {
    if (enzymes_P[row, 'Read'] > threshold*(enzymes_P[row+1, 'Read']+threshold*enzymes_P[row-1, 'Read'])/2) {
      enzymes_P[row, 'Remove'] = 'Yes'
    }
  }
  # Indicates when we have finished reading a well
  if (well_counter == 96) {
    # Loops through letters A-H
    for (platerow in LETTERS[1:8]) {
      # Gives subset of data that has only specific Sample ID and plate row
      enzymes_P_platerow <- enzymes_P[enzymes_P$SampleID == enzymes_P[row,'SampleID'] & enzymes_P$Row == platerow,]
      # Counts the number of Yes's from the specific Sample ID and plate row
      # If more than 4 in the row are Yes, we change entire row to Yes
      if (length(which(enzymes_P_platerow$Remove=="Yes")) >= 4) {
        enzymes_P[enzymes_P$SampleID == enzymes_P[row,'SampleID'] & enzymes_P$Row == platerow, 'Remove'] = 'Wipe'
      }
    }
    # Resets well_counter to 0 for the new plate
    well_counter = 0
  }
}

enzymes_clean <- enzymes_P %>%
  # Remove transfer errors that affected entire rows
  filter(Remove!="Wipe") %>%
  # Indicate buffer errors and change read to 0.05 when error is present
  mutate(Remove = ifelse(Col%in%c(1:4) & Remove=="Yes","buffer_error",Remove),
         Read = ifelse(Remove=="buffer_error",0.05,Read)) %>%
  # Filter out sample points that were marked as errors
  filter(Remove!="Yes") %>%
  rbind(filter(enzymes6,PlateType!="P"))

#################################################
# Calculate activity from fluorescence/absorbance readings

# Group together similar well types and calculate averages
Enz <- enzymes_clean %>% 
  group_by(assay_date,SampleID,AssayGroup,WellType,Variate,Concen) %>%
  summarize(Read.Mean = mean(Read,na.rm=TRUE)) %>%
  ungroup()

# Organize the well types for LAP assay
AMC <- filter(Enz,AssayGroup %in% c("AMC","AMCMUB")) %>% select(-AssayGroup)
AMC.Stan <- filter(AMC,WellType %in% c("Stan","Buffer")) %>% select(-SampleID,-Variate,-Concen) %>% spread(WellType,Read.Mean)
AMC.Quench <- filter(AMC,WellType %in% c("Quench","HomBlank")) %>% select(-Variate,-Concen) %>% spread(WellType,Read.Mean)
AMC.Assay <- filter(AMC,WellType=="Assay") %>% spread(WellType,Read.Mean)
AMC.SubCon <- filter(AMC,WellType=="SubCon") %>% select(-SampleID) %>% spread(WellType,Read.Mean)

# Organize the well types for other hydrolase assays
MUB <- filter(Enz,AssayGroup %in% c("MUB","AMCMUB")) %>% select(-AssayGroup)
MUB.Stan <- filter(MUB,WellType %in% c("Stan","Buffer")) %>% select(-SampleID,-Variate,-Concen) %>% spread(WellType,Read.Mean)
MUB.Quench <- filter(MUB,WellType %in% c("Quench","HomBlank")) %>% select(-Variate,-Concen) %>% spread(WellType,Read.Mean)
MUB.Assay <- filter(MUB,WellType=="Assay") %>% spread(WellType,Read.Mean)
MUB.SubCon <- filter(MUB,WellType=="SubCon") %>% select(-SampleID) %>% spread(WellType,Read.Mean)

# Read in and process dry weight data
Weights <-rbind(read.xlsx("extracellular_enzyme_data/enzyme_weights.xlsx",sheet="Natives"),
                read.xlsx("extracellular_enzyme_data/enzyme_weights.xlsx",sheet="Pines")) %>%
  select(sample_id,envelope_weight,envelope_fresh,envelope_dry,assay_weight)
names(Weights) <- c("SampleID","Env","EnvFresh","EnvDry","AssayWt")

# Calculate dry weight of assay sample
Weights <- filter(select(mutate(Weights,DryWt = AssayWt*(EnvDry - Env)/(EnvFresh - Env)),SampleID,DryWt),!is.na(DryWt))

# Merge hydrolase data with dry weights
AMC.MUB <- left_join(rbind(left_join(left_join(left_join(AMC.Assay,AMC.SubCon),AMC.Stan),AMC.Quench),
                           left_join(left_join(left_join(MUB.Assay,MUB.SubCon),MUB.Stan),MUB.Quench)),Weights)

# Calculate hydrolase activity in umol/h/g units based on 3.125 nmol standard and 0.125 ml homogenate, 4 hours, DryWt per 100 ml
Hydro <- mutate(AMC.MUB, Quench.coef = Quench/Stan) %>%
  mutate(Activity = ((Assay-HomBlank-SubCon*Quench.coef)/(Stan*Quench.coef/3.125))/(4*DryWt/100*0.125)/1000) %>%
  select(assay_date,SampleID,Concen,Variate,Activity)

# Organize the well types for oxidase assays
OX <- filter(Enz,AssayGroup=="Oxidase") %>% select(-AssayGroup)
OX.Blank <- filter(OX,WellType %in% c("Buffer","HomBlank")) %>% select(-Concen) %>% spread(WellType,Read.Mean)
OX.Assay <- filter(OX,WellType %in% c("Assay","SubCon")) %>% spread(WellType,Read.Mean)

# Calculate oxidase activity in umol/h/g units based on LDOPA EC of 7.9*125/300 (125 ml supernatant transferred), 0.125 ml homogenate, 24 hours, DryWt per 100 ml buffer
Oxi <- mutate(left_join(left_join(OX.Assay,OX.Blank),Weights), Activity = ((Assay-HomBlank-(SubCon-Buffer))*100)/(3.29*0.125*24*DryWt))
PER <- rename(select(filter(Oxi,Variate=="PER"),assay_date,SampleID,Concen,Activity),PER = Activity)
PPO <- rename(select(filter(Oxi,Variate=="PPO"),assay_date,SampleID,Concen,Activity),PPO = Activity)
Oxidase <- gather(mutate(full_join(PPO,PER),PER = PER - PPO),key=Variate,value=Activity,-assay_date,-SampleID,-Concen)

# Combine hydrolase and oxidase data
Enz_all <- rbind(Hydro,Oxidase) %>%
  filter(is.finite(Activity))

#################################################

# Fit data to Michaelis Menten curve
library(nlme)
library(nlstools)

# MM fit function
attach(Enz_all)
mmfit <- function(i){
  nls.m <- try(nls(y~Vmax*x/(Km+x),data=list(y=Activity[i],x=Concen[i]),
                   start=list(Vmax=max(Activity[i],na.rm=T),Km=max(Concen[i]/4,na.rm=T))))
  Vmax <- try(coef(nls.m)[[1]],silent=T)
  Km <- try(coef(nls.m)[[2]],silent=T)
  lowerV <- try(confint2(nls.m)[[1]],silent=T)
  upperV <- try(confint2(nls.m)[[3]],silent=T)
  lowerK <- try(confint2(nls.m)[[2]],silent=T)
  upperK <- try(confint2(nls.m)[[4]],silent=T)
  fit <- try(100*(upperV - lowerV)/Vmax,silent=T)
  model.pts <- try(Vmax*Concen[i]/(Km+Concen[i]),silent=T)
  plot(Activity[i]~Concen[i],main=paste("Sample ID:",SampleID[i][1],"Enzyme:",Variate[i][1],sep=" "))
  try(lines(model.pts~Concen[i],col="red"),silent=T)
  text(Concen[i],Activity[i],Concen[i],pos=4)
  mtext(paste("Vmax=",ifelse(is.numeric(Vmax),round(Vmax,3),NA)," Km=",ifelse(is.numeric(Km),round(Km,3),NA),
              " Fit=",ifelse(is.numeric(fit),round(fit,3),NA),sep=""))
  c(Vmax,Km,lowerV,upperV,lowerK,upperK)
}

# Apply function to values
pdf("fits.pdf")
x.seq <- seq(len=dim(Enz_all)[1])
fit <- tapply(x.seq, list(Enz_all$assay_date, Enz_all$SampleID, Enz_all$Variate), 
              mmfit)
dev.off()

Enz_fit <- as.data.frame.table(tapply(x.seq, list(Enz_all$assay_date, Enz_all$SampleID, Enz_all$Variate),function(i)1))
names(Enz_fit) <- c("assay_date","SampleID","Variate","Vmax")
detach(Enz_all)

# Extract successful Vmax and Km values
Enz_fit$Vmax <- as.numeric(unlist(sapply(sapply(fit,'[',1),'[[',1) %>% replace_na()))
Enz_fit$Km <- as.numeric(unlist(sapply(sapply(fit,'[',2),'[[',1) %>% replace_na()))
Enz_fit$lowerV <- as.numeric(unlist(sapply(sapply(fit,'[',3),'[[',1) %>% replace_na()))
Enz_fit$upperV <- as.numeric(unlist(sapply(sapply(fit,'[',4),'[[',1) %>% replace_na()))
Enz_fit$lowerK <- as.numeric(unlist(sapply(sapply(fit,'[',5),'[[',1) %>% replace_na()))
Enz_fit$upperK <- as.numeric(unlist(sapply(sapply(fit,'[',6),'[[',1) %>% replace_na()))

# Reorganize dataset and add back samples that failed the fit
# Full list samples
Enz_all2 <- Enz_all %>%
  select(assay_date,SampleID, Variate) %>%
  unique()

Enz_fit2 <- Enz_fit %>%
  filter(!is.na(Vmax)) %>%
  mutate(SampleID = as.numeric(as.character(SampleID)),
         Variate = as.character(Variate),
         assay_date = ymd(as.character(assay_date))) %>%
  right_join(Enz_all2,by=c("assay_date","SampleID","Variate"))

# Set Vmaxs that are negative to 0 and their Kms to NA
Enz_fit3 <- Enz_fit2 %>%
  mutate(Vmax = ifelse(Vmax<0,0,Vmax),
         lowerV = ifelse(Vmax==0,0,lowerV),
         upperV = ifelse(Vmax==0,0,upperV),
         Km = ifelse(Vmax==0,NA,Km),
         lowerK = ifelse(Vmax==0,NA,lowerK),
         upperK = ifelse(Vmax==0,NA,upperK))

# Export file
write_csv(Enz_fit3, "fitted_enzyme_data.csv")

#################################################
# The next section of code visualizes the samples that are still having fit issues

# Check samples that failed fit
Enz_fit_fail <- Enz_fit2 %>% 
  filter(is.na(Vmax)) %>%
  select(assay_date,SampleID,Variate) %>%
  unique() %>%
  # Add points that were used in MM fit
  left_join(Enz_all,by=c("assay_date","SampleID","Variate"))

# Determine which of the fit fails were just due to there being no detectable activity
Enz_fit_fail_0 <- Enz_fit_fail %>%
  # Take the average of the highest 3-4 concentration points
  filter(Concen > 60) %>%
  group_by(assay_date,SampleID,Variate) %>%
  summarize(Avg_activity = mean(Activity,na.rm=TRUE)) %>%
  # When averages are less than 0.05 for BX and BG and 0.2 for the remaining enzymes, the plots look like noise
  filter(Avg_activity < 0.2 & Variate%in%c("AG","AP","BG","CBH","LAP","PER","PPO") |
           Avg_activity < 0.05 & Variate%in%c("BX","NAG")) %>%
  ungroup() %>%
  select(-Avg_activity) %>%
  unique()
# Vmax can be set to 0 for these samples and enzymes

# Use the Enz0 if you want to visualize the plots
Enz_0 <- Enz_fit_fail_0 %>%
  left_join(Enz_all,by=c("assay_date","SampleID","Variate"))

# Determine which of the fit fails had detectable activity and need to be cleaned again for proper fits
Enz_fit_fail_1 <- Enz_fit_fail %>%
  filter(Concen > 60) %>%
  group_by(assay_date,SampleID,Variate) %>%
  summarize(Avg_activity = mean(Activity,na.rm=TRUE)) %>%
  filter(Avg_activity >= 0.2 & Variate%in%c("AG","AP","BG","CBH","LAP","PER","PPO") |
           Avg_activity >= 0.05 & Variate%in%c("BX","NAG")) %>%
  ungroup() %>%
  select(-Avg_activity) %>%
  unique() 
Enz_reclean <- Enz_fit_fail_1 %>%
  left_join(Enz_all,by=c("assay_date","SampleID","Variate"))

# Visualize
ggplot(filter(Enz_reclean,Variate=="AG"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="AP"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="BG"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="BX"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="CBH"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="LAP"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="NAG"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="PER"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)
ggplot(filter(Enz_reclean,Variate=="PPO"),aes(Concen,Activity)) + 
  geom_point() + facet_wrap(SampleID~.)


