#Cleaning the data table All Samples----

library(mediation)

set.seed(1000)

m.otu <- Oct23_16S.056.a.micro2.c$otu_table
m.meta <- Oct23_16S.056.a.micro2.c$sample_table
m.taxa <- Oct23_16S.056.a.micro2.c$tax_table
m.taxa.abund <- Oct23_16S.056.a.micro2.c$taxa_abund
m.taxa.abund.genus <- m.taxa.abund$Genus
m.taxa.abund.species <- m.taxa.abund$Species
m.taxa.abund.genus <- as.data.frame(t(m.taxa.abund.genus))
m.taxa.abund.species <- as.data.frame(t(m.taxa.abund.species))

#merging metadata and taxa abundance table (not changing "label" genus)
m.genus <- cbind(m.meta, m.taxa.abund.species)
#replace - with space
m.genus$Age <- gsub("-", "", m.genus$Age, fixed = TRUE) 
#replace - with space
m.genus$Semen_Volume <- gsub("-", "", m.genus$Semen_Volume, fixed = TRUE)

#Check Class for each column
m.class.check <- as.data.frame(sapply(m.genus, class))

#make numeric for Age
m.genus$Age <- as.numeric(m.genus$Age)

#make numeric for Normal Fertilization Rate
m.genus$Normal_Fertilization_Rate <- as.numeric(m.genus$Normal_Fertilization_Rate)
m.genus$Normal_Fertilization_Rate <- as.numeric(m.genus$Normal_Fertilization_Rate)

#make numeric for Semen Volume
m.genus$Semen_Volume <- as.numeric(m.genus$Semen_Volume)

#count number of NA
sum(is.na(m.genus))

#NA did not affect the Mediation lm and Outcome lm models, but affect 
#the mediation analysis package

#omit patient with NA data
m.genus.o <- na.omit(m.genus)

sum(is.na(m.genus.o)) #omitted 9 patients' data

#When no. of cells that has a non-zero value is less than or equal to 1,
#the mediation analysis did not work
#so, need to remove the column that has
# " no. of cells that has a non-zero value is less than or equal to 1

sum(m.genus.o[,225] > 0) <= 1

#Check no. of cells that has a non-zero value is less than or equal to 1 
#for @ column, similar as above
s <- as.data.frame(colSums(m.genus.o !=1))

#Microbiome Analysis - Control Sample----
m.genus.o2 <- m.genus.o[which(m.genus.o$Fertility=='Control'),]

#removal of column which has 1 cells or less containing zero value
m.genus.o2 <- m.genus.o2[, colSums(m.genus.o2 !=0) >= 2]

##DEFB119 > Microbiome > SemenParameter-Forward----
m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[13]
    Mediator <- colnames(m.genus.o2)[i]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.con_def_Micro_sem_f <- as.data.frame(m.vec)
write_xlsx(m.vec.con_def_Micro_sem_f,
           "DEFB119_Micro_Sem_Con_For.xlsx")


##DEFB119 > SemenParameter > Microbiome-Reverse----
m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[13]
    Mediator <- colnames(m.genus.o2)[i]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.con_def_Micro_sem_r <- as.data.frame(m.vec)
write_xlsx(m.vec.con_def_Micro_sem_r,
           "DEFB119_Micro_Sem_Con_Rev.xlsx")


###############################
###############################
#Microbiome Analysis - Infertile Sample----
m.genus.o2 <- m.genus.o[which(m.genus.o$Fertility=='Infertile'),]

#removal of column which has 1 cells or less containing zero value
m.genus.o2 <- m.genus.o2[, colSums(m.genus.o2 !=0) >= 2]


##DEFB119 > Microbiome > SemenParameter-Forward----
m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[13]
    Mediator <- colnames(m.genus.o2)[i]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.In_def_Micro_sem_f <- as.data.frame(m.vec)
write_xlsx(m.vec.In_def_Micro_sem_f,
           "DEFB119_Micro_Sem_In_For.xlsx")


##DEFB119 > SemenParameter > Microbiome-Reverse----
m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[13]
    Mediator <- colnames(m.genus.o2)[i]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.In_def_Micro_sem_r <- as.data.frame(m.vec)
write_xlsx(m.vec.In_def_Micro_sem_r,
           "DEFB119_Micro_Sem_In_Rev.xlsx")



##_________________________#####

#Microbiome as Treatment###########
##Control Sample----
m.genus.o2 <- m.genus.o[which(m.genus.o$Fertility=='Control'),]

#removal of column which has 1 cells or less containing zero value
m.genus.o2 <- m.genus.o2[, colSums(m.genus.o2 !=0) >= 2]

###Microbiome > Defb119 > SemenParameter-Forward----

m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[13]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.con_micro_def_sem_f <- as.data.frame(m.vec)
write_xlsx(m.vec.con_micro_def_sem_f,
           "Micro_DEFB119_Semen_Con_For.xlsx")


###Microbiome > SemenParameter > DEFB119 -Rev----

m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[j]
    Outcome <- colnames(m.genus.o2)[13]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.con_micro_def_sem_r <- as.data.frame(m.vec)
write_xlsx(m.vec.con_micro_def_sem_r,
           "Micro_DEFB119_Semen_Con_Rev.xlsx")

##Infertile Sample----
m.genus.o2 <- m.genus.o[which(m.genus.o$Fertility=='Infertile'),]

#removal of column which has 1 cells or less containing zero value
m.genus.o2 <- m.genus.o2[, colSums(m.genus.o2 !=0) >= 2]

###Microbiome > Defb119 > SemenParameter-Forward----

m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[13]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.In_micro_def_sem_f <- as.data.frame(m.vec)
write_xlsx(m.vec.In_micro_def_sem_f,
           "Micro_DEFB119_Semen_In_For.xlsx")

###Microbiome > SemenParameter > DEFB119 -Rev----

m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[j]
    Outcome <- colnames(m.genus.o2)[13]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.In_micro_def_sem_r <- as.data.frame(m.vec)
write_xlsx(m.vec.In_micro_def_sem_r,
           "Micro_DEFB119_Semen_In_Rev.xlsx")


###________________________########


#SemenParameter as Treatment###########
##Control Sample----
m.genus.o2 <- m.genus.o[which(m.genus.o$Fertility=='Control'),]

#removal of column which has 1 cells or less containing zero value
m.genus.o2 <- m.genus.o2[, colSums(m.genus.o2 !=0) >= 2]

###Semen Parameter > DEFB119 > Microbiome -Forward----

m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[13]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: SemenParameter > Microbiome
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: SemenParameter -> DEFB119
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}



library(writexl)

m.vec.con_sem_def_micro_f <- as.data.frame(m.vec)
write_xlsx(m.vec.con_sem_def_micro_f,
           "Semen_DEFB119_Micro_Con_For.xlsx")


###SemenParameter > Microbiome > DEFB119-Rev----

m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[j]
    Outcome <- colnames(m.genus.o2)[13]
    
    #Mediation lm: SemenParameter > Microbiome
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: SemenParameter -> DEFB119
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.con_sem_def_micro_r <- as.data.frame(m.vec)
write_xlsx(m.vec.con_sem_def_micro_r ,
           "Semen_DEFB119_Micro_Con_Rev.xlsx")


##Infertile Sample----
m.genus.o2 <- m.genus.o[which(m.genus.o$Fertility=='Infertile'),]

#removal of column which has 1 cells or less containing zero value
m.genus.o2 <- m.genus.o2[, colSums(m.genus.o2 !=0) >= 2]
###Semen Parameter > DEFB119 > Microbiome -Forward----

m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[13]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: SemenParameter > Microbiome
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: SemenParameter -> DEFB119
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}



library(writexl)

m.vec.In_sem_def_micro_f <- as.data.frame(m.vec)
write_xlsx(m.vec.In_sem_def_micro_f,
           "Semen_DEFB119_Micro_In_For.xlsx")


###SemenParameter > Microbiome > DEFB119-Rev----

m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[j]
    Outcome <- colnames(m.genus.o2)[13]
    
    #Mediation lm: SemenParameter > Microbiome
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: SemenParameter -> DEFB119
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.In_sem_def_micro_r <- as.data.frame(m.vec)
write_xlsx(m.vec.In_sem_def_micro_r ,
           "Semen_DEFB119_Micro_In_Rev.xlsx")







#Overal Mediation Analysis without in groups ----
library(mediation)

set.seed(1000)

m.otu <- Oct23_16S.056.a.micro2.c$otu_table
m.meta <- Oct23_16S.056.a.micro2.c$sample_table
m.taxa <- Oct23_16S.056.a.micro2.c$tax_table
m.taxa.abund <- Oct23_16S.056.a.micro2.c$taxa_abund
m.taxa.abund.genus <- m.taxa.abund$Genus
m.taxa.abund.species <- m.taxa.abund$Species
m.taxa.abund.genus <- as.data.frame(t(m.taxa.abund.genus))
m.taxa.abund.species <- as.data.frame(t(m.taxa.abund.species))

#merging metadata and taxa abundance table (not changing "label" genus)
m.genus <- cbind(m.meta, m.taxa.abund.species)
#replace - with space
m.genus$Age <- gsub("-", "", m.genus$Age, fixed = TRUE) 
#replace - with space
m.genus$Semen_Volume <- gsub("-", "", m.genus$Semen_Volume, fixed = TRUE)

#Check Class for each column
m.class.check <- as.data.frame(sapply(m.genus, class))

#make numeric for Age
m.genus$Age <- as.numeric(m.genus$Age)

#make numeric for Normal Fertilization Rate
m.genus$Normal_Fertilization_Rate <- as.numeric(m.genus$Normal_Fertilization_Rate)
m.genus$Normal_Fertilization_Rate <- as.numeric(m.genus$Normal_Fertilization_Rate)

#make numeric for Semen Volume
m.genus$Semen_Volume <- as.numeric(m.genus$Semen_Volume)

#count number of NA
sum(is.na(m.genus))

#NA did not affect the Mediation lm and Outcome lm models, but affect 
#the mediation analysis package

#omit patient with NA data
m.genus.o <- na.omit(m.genus)

sum(is.na(m.genus.o)) #omitted 9 patients' data

#When no. of cells that has a non-zero value is less than or equal to 1,
#the mediation analysis did not work
#so, need to remove the column that has
# " no. of cells that has a non-zero value is less than or equal to 1

##Prepare the data----
m.genus.o2 <- m.genus.o

#removal of column which has 1 cells or less containing zero value
m.genus.o2 <- m.genus.o2[, colSums(m.genus.o2 !=0) >= 2]

#DEFB119 as Treatment----
##DEFB119 > Microbiome > SemenParameter-Forward----
m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[13]
    Mediator <- colnames(m.genus.o2)[i]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.def_Micro_sem_f <- as.data.frame(m.vec)
write_xlsx(m.vec.def_Micro_sem_f,
           "DEFB119_Micro_Sem_For.xlsx")


##DEFB119 > SemenParameter > Microbiome-Reverse----
m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[13]
    Mediator <- colnames(m.genus.o2)[i]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.def_Micro_sem_r <- as.data.frame(m.vec)
write_xlsx(m.vec.def_Micro_sem_r,
           "DEFB119_Micro_Sem_Rev.xlsx")




#Microbiome as Treatment----
###Microbiome > Defb119 > SemenParameter-Forward----

m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[13]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.micro_def_sem_f <- as.data.frame(m.vec)
write_xlsx(m.vec.micro_def_sem_f,
           "Micro_DEFB119_Semen_For.xlsx")


###Microbiome > SemenParameter > DEFB119 -Rev----

m.vec <- NULL

for ( i in 17:38){
  for ( j in 8:12){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[j]
    Outcome <- colnames(m.genus.o2)[13]
    
    #Mediation lm: Microbiome > Defb119 >
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: Microbiome -> Semen parameters
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.micro_def_sem_r <- as.data.frame(m.vec)
write_xlsx(m.vec.micro_def_sem_r,
           "Micro_DEFB119_Semen_Rev.xlsx")





#SemenParameter as Treatment----
###Semen Parameter > DEFB119 > Microbiome -Forward----

m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[13]
    Outcome <- colnames(m.genus.o2)[j]
    
    #Mediation lm: SemenParameter > Microbiome
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: SemenParameter -> DEFB119
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}



library(writexl)

m.vec.sem_def_micro_f <- as.data.frame(m.vec)
write_xlsx(m.vec.sem_def_micro_f,
           "Semen_DEFB119_Micro_For.xlsx")


###SemenParameter > Microbiome > DEFB119-Rev----

m.vec <- NULL

for ( i in 8:12){
  for ( j in 17:38){
    
    Treat <- colnames(m.genus.o2)[i]
    Mediator <- colnames(m.genus.o2)[j]
    Outcome <- colnames(m.genus.o2)[13]
    
    #Mediation lm: SemenParameter > Microbiome
    med_formula = as.formula(paste('`',Mediator,'`', "~", '`',
                                   Treat,'`', "+", "Age", sep = ""))
    
    mediate_micro_def <- lm(med_formula, data=m.genus.o2)
    
    
    ###Outcome lm: SemenParameter -> DEFB119
    out_formula = as.formula(paste('`',Outcome,'`', "~",
                                   '`',Treat,'`', "+",'`', Mediator,'`', "+", 
                                   "Age", sep = ""))
    
    out_micro_def <- lm(out_formula, data=m.genus.o2)
    
    
    
    ###Mediation- Med vs Out lm
    results=mediate(mediate_micro_def, out_micro_def, treat=Treat,
                    mediator=Mediator, 
                    robustSE = TRUE, boot=F, sims = 1000)
    
    m.summary = (summary(results))
    
    id <- paste("T:", Treat,">", "M:",Mediator,">",
                "O:", Outcome, sep="") 
    
    m.res <- capture.output(m.summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7]) grab ACME p value
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME <- tmp[(2)]
    
    tmp <-  base::strsplit(m.res[grep("ACME",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ACME.p <- tmp[length(tmp)]
    
    tmp <-  base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    m.ADE <- tmp[(2)]
    
    tmp <- base::strsplit(m.res[grep("ADE",m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    i_str = which(grepl("Effect", tmp))
    m.totale <- tmp[(i_str + 1)]
    
    tmp <- base::strsplit(m.res[grep("Total Effect", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]  # remove stars in case the last element is star
    m.totale.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(m.res[grep("Prop. Mediated", m.res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    m.prop.mediated <- tmp[(i_str + 1)]
    
    m.vec.a = c(id, m.ACME, m.ACME.p, m.ADE, m.ADE.p, m.totale, m.totale.p, m.prop.mediated)
    names(m.vec.a) <- c("Treat_Mediator_Outcome", "ACME", "ACME.p", "ADE", "ADE.p", "Total.Effect",
                        "Total.Effect.p",  "prop.mediated")
    
    #Glue data together
    m.vec = rbind(m.vec, m.vec.a)
    
  }}

library(writexl)

m.vec.sem_def_micro_r <- as.data.frame(m.vec)
write_xlsx(m.vec.sem_def_micro_r ,
           "Semen_DEFB119_Micro_Rev.xlsx")






