#Example code for the main analyses described in Duffy et al. (Development of a human genetics-guided priority score for 19,365 genes and 347 drug indications). 2023.
							      
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(caret)
library(parallel)
library(logistf)
							      
#mi - binary variable used in dataset for drug indication presence.  
#genetic features used which comprise the GPS
geneticpredictors=c('clinvar', 'hgmd', 'omim','geneburden', 'singlevariant','eqtl_phenotype','locus2gene','pqtl_phenotype')
geneticpredictors_doe=c('clinvar_doe', 'hgmd_doe', 'omim_doe','geneburden_doe', 'singlevariant_doe','eqtl_phenotype_doe','locus2gene_doe','pqtl_phenotype_doe')

#phecode categories used as covariates in regression models
covariates=c('categorycongenital_anomalies','categorydermatologic','categorydigestive','categoryendocrine_metabolic','categorygenitourinary','categoryhematopoietic','categorymental_disorders','categorymusculoskeletal','categoryneurological','categoryrespiratory','categorysense_organs','categorysymptoms')
covariates=c('number_gene_targets','oe_dichotomized', covariates)

# Analysis 1a-1e: Code required to run five-fold cross-validation (CV) model to construct the genetic priority score (GPS) using the Open Targets dataset and then applied to Sider and the all genes dataset (19,365 genes and 348 drug indications).

# Analysis 1a: Get weights for each of the 5-CV Open target datasets 

Firthreg_weights<-mclapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  OT_dataset=fread(paste0('OT_drugdataset_80_CV', CVsample, '.txt'),data.table=F) #80% training Open target dataset 
  #genetic features + phecode category covariates
  Predictors=paste(paste(geneticpredictors,collapse='+'),'+', paste(covariates,collapse='+'),collapse='+')
   ##Run firth regression across each CV training dataset to get weights
  firth_mod <- try(logistf(as.formula(paste0('mi ~ ',Predictors)), data=OT_dataset, firth = TRUE))
  results <- cbind(beta=coef(firth_mod)[-1],lowerCI= firth_mod$ci.lower[-1], upperCI=firth_mod$ci.upper[-1],  P.val =firth_mod$prob[-1])
  write.table(results, paste0('Firth_weights_Opentargets_',CVsample,'.txt'), sep='\t',quote=F)
}, mc.cores=10)

#Analysis 1b: Use weights from Firth and create score using remaining 20% of data for each CV

Genescore_sum<-lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){

  Firth<-fread(paste0('Firth_weights_Opentargets_',CVsample,'.txt'),data.table=F) #weights from CV sample
  Firth=Firth %>% select(V1, beta) %>% rename(predictor=V1)
  Firth_weights=as.data.frame(t(Firth$beta))
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F) # OT remaining 20% test set 
  #extract gene, parentterm and genetic predictor columns
  OT_dataset_20_gene_phenotype=OT_dataset_20[, grepl(paste0('\\bgene\\b|parentterm|genescores_|phenotype_'), colnames(OT_dataset_20))] %>% distinct() 
  ## multiply betas from firth with values from 20% dataset
  Genescores_beta_genetic_values=mapply("*", OT_dataset_20_gene_phenotype[intersect(names(OT_dataset_20_gene_phenotype), names(Firth_weights))],
    Firth_weights[intersect(names(OT_dataset_20_gene_phenotype), names(Firth_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'), sep='\t', row.names=F, quote=F )

})

#Analysis 1c: Choose the CV sample with the max OR for validation set

max_filetype<-do.call(rbind,lapply(c(paste0('CVsample',seq(1:5))), function(samplecv){   
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F)
  genescorefile=fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'), data.table=F)
  # combine genescore sum file with drug and mi data for each test set 
  Dataset_genescores=merge(OT_dataset_20[c('drugname','gene','parentterm','category','mi')] ,genescorefile, by=c('gene', 'parentterm') )
  #run logistic model 
  model1 =glm(as.formula(paste0('mi ~ genescoresum + category + number_gene_targets')), data=Dataset_genescores,family = 'binomial')
  mod_output<-rbind(cbind.data.frame(CV=samplecv,OR=exp(summary(model2)$coefficient[2,1]),lowerCI=exp(summary(model2)$coefficient[2,1]-(1.96* summary(model2)$coefficient[2,2])),upperCI=exp(summary(model2)$coefficient[2,1]+(1.96* summary(model2)$coefficient[2,2])),P.val=summary(model2)$coefficient[2,4]))
  return(mod_output)
  }))

## take CV with the max OR
max_filetype<-max_filetype[order(max_filetype$OR, decreasing=T),]
OT_CV=as.character(max_filetype$CV[1])

#Analysis 1d: Combine 5-CV 20% datasets to create one dataset

#scores with predictors 
Combine_genescores<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'),data.table=F)
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','number_gene_targets','parentterm','category','mi')], Genescore_sumfile_20 ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores1<-do.call(rbind.fill,Combine_genescores)
write.table(Combine_genescores1, gzfile(paste0('All_genescoresum_across_all_predictors_opentargets.txt.gz')), sep='\t',quote=F,row.names=F)

#scores without predictors 
Combine_genescores_nopredictors<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'),data.table=F)
  Genescore_sumfile_20_nopred<-Genescore_sumfile_20 %>% distinct(gene,parentterm,genescoresum )
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','number_gene_targets','parentterm','category','mi')], Genescore_sumfile_20_nopred ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores_nopredictors1<-do.call(rbind.fill,Combine_genescores_nopredictors)
write.table(Combine_genescores_nopredictors1, gzfile(paste0('All_genescoresum_opentargets.txt.gz')), sep='\t',quote=F,row.names=F)

# Analysis 1e: Calculate GPS in Sider dataset and all genes dataset using the OT firth weights from Analysis 3 which gave the max OR 

samplecv=OT_CV
Firth<-fread(paste0('Firth_weights_Opentargets_',samplecv,'.txt'),data.table=F)
Firth=Firth %>% select(V1, beta) %>% rename(predictor=V1)

Genescore_sum<-lapply(c('Sider','Allgenes'), function(valdataset){
  Validation_dataset=fread(paste0(valdataset,'_drugdataset.txt'),data.table=F) #sider/allgenes dataset 
  Firth_weights=as.data.frame(t(Firth$beta))
  Validation_dataset_gene_phenotype=Validation_dataset[, grepl(paste0('\\bgene\\b|parentterm|genescores_|phenotype_'), colnames(Validation_dataset))] %>% distinct()
  ## multiply betas from firth with values from Sider/Allgenes dataset
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))],
    Firth_weights[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0('All_genescoresum_across_all_predictors_',valdataset,'.txt'), sep='\t', row.names=F, quote=F )
    #for sider add GPS sum to the sider-drug dataset with mi
  if(valdataset=='Sider'){
  Genescores_beta_genetic_genept=Validation_dataset %>% distinct(gene, parentterm, genescoresum) 
  genescorefile_drugs=inner_join(Validation_dataset[c('drugname','gene','number_gene_targets','parentterm','category','mi')],Genescores_beta_genetic_genept, by=c('gene', 'parentterm'))
  write.table(genescorefile_drugs, paste0('All_genescoresum_across_drugs_',valdataset,'.txt'), sep='\t', row.names=F, quote=F )
  }
})

# Analysis 2a-2e: Create GPS-D 
#GOF annotated as -1, LOF annotated as +1 and missing/no evidence/neutral annotated as 0

# Analysis 2a: Get weights for each of the 5-CV Open target datasets 
Firthreg_weights_doe<-mclapply(c(paste0('CVsample',rep(1:5))), function(CVsample){

  OT_dataset_doe=fread(paste0('OT_drugdataset_80_CV', CVsample, '_doe.txt'),data.table=F) #80% training Open target dataset 
  OT_dataset_doe[8:15][OT_dataset_doe[8:15]=='GOF'] <- 1
  OT_dataset_doe[8:15][OT_dataset_doe[8:15]=='LOF'] <-1
  OT_dataset_doe[8:15][OT_dataset_doe[8:15]=='Neutral'] <-0
 
  #genetic features + phecode category covariates
  Predictors_doe=paste(paste(geneticpredictors_doe,collapse='+'), '+', paste(covariates,collapse='+'),collapse='+' )
   ##Run firth regression across each CV training dataset to get weights - don't take DOE into account here
  firth_mod <- try(logistf(as.formula(paste0('mi ~ ',Predictors_doe)), data=OT_dataset_doe, firth = TRUE))
  results <- cbind(beta=coef(firth_mod)[-1],lowerCI= firth_mod$ci.lower[-1], upperCI=firth_mod$ci.upper[-1],  P.val =firth_mod$prob[-1])
  write.table(results, paste0('Firth_weights_Opentargets_',CVsample,'_doe.txt'), sep='\t',quote=F)
}, mc.cores=10)

#Analysis 2b: Use weights from Firth and create score using remaining 20% of data for each CV
Genescore_sum<-lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  Firth<-fread(paste0('Firth_weights_Opentargets_',CVsample,'_doe.txt'),data.table=F) #weights from CV sample
  Firth=Firth %>% select(V1, beta) %>% rename(predictor=V1)
  Firth_weights=as.data.frame(t(Firth$beta))
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F) # OT remaining 20% test set 
  OT_dataset_20[8:15][OT_dataset_20[8:15]=='GOF'] <- -1 
  OT_dataset_20[8:15][OT_dataset_20[8:15]=='LOF'] <-1
  OT_dataset_20[8:15][OT_dataset_20[8:15]=='Neutral'] <-0
  #extract gene, parentterm and genetic predictor columns
  OT_dataset_20_gene_phenotype=OT_dataset_20[, grepl(paste0('\\bgene\\b|parentterm|genescores_|phenotype_'), colnames(OT_dataset_20))] %>% distinct() 
  ## multiply betas from firth with values from 20% dataset
  Genescores_beta_genetic_values=mapply("*", OT_dataset_20_gene_phenotype[intersect(names(OT_dataset_20_gene_phenotype), names(Firth_weights))],
    Firth_weights[intersect(names(OT_dataset_20_gene_phenotype), names(Firth_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'), sep='\t', row.names=F, quote=F )

})

#Analysis 2c: Choose the CV sample with the max OR for validation set

max_filetype_doe<-do.call(rbind,lapply(c(paste0('CVsample',seq(1:5))), function(samplecv){   
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '_doe.txt'),data.table=F)
  genescorefile=fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'), data.table=F)
  # combine genescore sum file with drug and mi data for each test set 
  Dataset_genescores=merge(OT_dataset_20[c('drugname','gene','parentterm','category','mi')] ,genescorefile, by=c('gene', 'parentterm') )
  #run logistic model 
  model1 =glm(as.formula(paste0('mi ~ abs(genescoresum) + category + number_gene_targets')), data=Dataset_genescores,family = 'binomial')
  mod_output<-rbind(cbind.data.frame(CV=samplecv,OR=exp(summary(model2)$coefficient[2,1]),lowerCI=exp(summary(model2)$coefficient[2,1]-(1.96* summary(model2)$coefficient[2,2])),upperCI=exp(summary(model2)$coefficient[2,1]+(1.96* summary(model2)$coefficient[2,2])),P.val=summary(model2)$coefficient[2,4]))
  return(mod_output)
  }))

## take CV with the max OR
max_filetype_doe<-max_filetype_doe[order(max_filetype_doe$OR, decreasing=T),]
OT_CV_doe=as.character(max_filetype_doe$CV[1])

#Analysis 2d: Combine 5-CV 20% datasets to create one dataset

#scores with predictors 
Combine_genescores_doe<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '_doe.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'),data.table=F)
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','number_gene_targets','parentterm','category','mi')], Genescore_sumfile_20 ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores_doe1<-do.call(rbind.fill,Combine_genescores_doe)
write.table(Combine_genescores_doe1, gzfile(paste0('All_genescoresum_across_all_predictors_opentargets_doe.txt.gz')), sep='\t',quote=F,row.names=F)

#scores without predictors 
Combine_genescores_nopredictors_doe<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '_doe.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'),data.table=F)
  Genescore_sumfile_20_nopred<-Genescore_sumfile_20 %>% distinct(gene,parentterm,genescoresum )
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','number_gene_targets','parentterm','category','mi')], Genescore_sumfile_20_nopred ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores_nopredictors_doe1<-do.call(rbind.fill,Combine_genescores_nopredictors_doe)
write.table(Combine_genescores_nopredictors_doe1, gzfile(paste0('All_genescoresum_opentargets_doe.txt.gz')), sep='\t',quote=F,row.names=F)

# Analysis 2e - Calculate GPS in Sider dataset and all genes dataset using the OT firth weights from Analysis 3 which gave the max OR 

samplecv=OT_CV_doe
Firth_doe<-fread(paste0('Firth_weights_Opentargets_',samplecv,'_doe.txt'),data.table=F)
Firth_doe=Firth_doe %>% select(V1, beta) %>% rename(predictor=V1)

Genescore_sum<-lapply(c('Sider','Allgenes'), function(valdataset){
  Validation_dataset=fread(paste0(valdataset,'_drugdataset_doe.txt'),data.table=F) #sider/allgenes dataset 
  Validation_dataset[8:15][Validation_dataset[8:15]=='GOF'] <- -1 
  Validation_dataset[8:15][Validation_dataset[8:15]=='LOF'] <-1
  Validation_dataset[8:15][Validation_dataset[8:15]=='Neutral'] <-0
  Firth_weights=as.data.frame(t(Firth_doe$beta))
  Validation_dataset_gene_phenotype=Validation_dataset[, grepl(paste0('\\bgene\\b|parentterm|genescores_|phenotype_'), colnames(Validation_dataset))] %>% distinct()
  ## multiply betas from firth with values from Sider/Allgenes dataset
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))],
    Firth_weights[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0('All_genescoresum_across_all_predictors_',valdataset,'_doe.txt'), sep='\t', row.names=F, quote=F )
    #for sider add GPS sum to the sider-drug dataset with mi
  if(valdataset=='Sider'){
  Genescores_beta_genetic_genept=Validation_dataset %>% distinct(gene, parentterm, genescoresum) 
  genescorefile_drugs=inner_join(Validation_dataset[c('drugname','gene','number_gene_targets', 'parentterm','category','mi')],Genescores_beta_genetic_genept, by=c('gene', 'parentterm'))
  write.table(genescorefile_drugs, paste0('All_genescoresum_across_drugs_',valdataset,'_doe.txt'), sep='\t', row.names=F, quote=F )
  }
})


### Analysis 3 - (data for fig 4 and supplementary fig 5) Binned association analysis - GPS at 0.3 increments, logistic regression for Opentarget and Sider datasets

lapply(c('Opentargets','Sider'), function(datafile){

  if(datafile=='Opentargets'){
    genescoredataset=fread('All_genescoresum_opentargets.txt', data.table=F) #GPS OT data 
  } else {
    genescoredataset=fread('All_genescoresum_across_drugs_Sider.txt', data.table=F)#GPS Sider data 
  }
  #order by increasing GPS and add percentile column
  genescoredataset = genescoredataset %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
  # run logistic model for each 0.3 increment of score against data with GPS =0
  
  binned_gps<-do.call(rbind,lapply(c(seq(0.3,2.1,0.3)), function(bin) {
    df_filtered=rbind(genescoredataset %>% filter(genescoresum>=bin ) %>% mutate(bin=1),
        genescoredataset %>% filter(genescoresum==0 ) %>% mutate(bin=0))
    model1 =glm(as.formula(paste0('mi ~ bin + category + number_gene_targets')), data=df_filtered,family = 'binomial')
    mod_output<-rbind(cbind.data.frame( genescoresum =bin, Percentile=round(df_filtered$percent[1],2) , genes=length(unique(df_filtered$gene[df_filtered$bin==1])), parentterms=length(unique(df_filtered$parentterm[df_filtered$bin==1])), OR=exp(summary(model1)$coefficient[2,1]),lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4])) 
    },mc.cores=10)) 
    write.table(binned_gps, paste0('Binned_by_sum_binsize0.3_',datafile,'.txt'),sep='\t', quote=F, row.names=F)
})

### Analysis 4 - (data for fig 4 and supplementary fig 5) Binned association analysis Doe- GPS at 0.3 increments, logistic regression for Opentarget and Sider datasets
  # absolute value for all gps-d, +1 values for lof scores and -1 for gof 

lapply(c('Opentargets','Sider'), function(datafile){

  if(datafile=='Opentargets'){
    genescoredataset=fread('All_genescoresum_opentargets_doe.txt', data.table=F) #GPS OT data 
  } else {
    genescoredataset=fread('All_genescoresum_across_drugs_Sider_doe.txt', data.table=F)#GPS Sider data 
  }

  binned_gpsscore<-do.call(rbind,lapply(c('all','gof', 'lof'), function(score){
  #order by increasing GPS and add percentile column
  if(score=='all'){
      Dataset_genescores1= genescoredataset %>%mutate(genescoresum=abs(genescoresum)) %>% arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
      } else if(score=='gof') {
        genescoredataset$genescoresum_gof= ifelse(genescoredataset$genescoresum >0,  genescoredataset$genescoresum , 0)
          Dataset_genescores1= genescoredataset%>% select(-genescoresum) %>% rename(genescoresum=genescoresum_gof) %>% arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
      } else {
          genescoredataset$genescoresum_lof= ifelse(genescoredataset$genescoresum < 0,  genescoredataset$genescoresum , 0)
          Dataset_genescores1= genescoredataset%>% select(-genescoresum) %>% rename(genescoresum=genescoresum_lof) %>% mutate(genescoresum=genescoresum*-1) %>% arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
      }
  # run logistic model for each 0.3 increment of score against data with GPS =0
  
  binned_gps<-do.call(rbind,lapply(c(seq(0.3,2.1,0.3)), function(bin) {
    df_filtered=rbind(Dataset_genescores1 %>% filter(genescoresum>=bin ) %>% mutate(bin=1),
        Dataset_genescores1 %>% filter(genescoresum==0 ) %>% mutate(bin=0))
    model1 =glm(as.formula(paste0('mi ~ bin + category + number_gene_targets')), data=df_filtered,family = 'binomial')
    mod_output<-rbind(cbind.data.frame( genescoresum =bin, Doe=score, Percentile=round(df_filtered$percent[1],2) , genes=length(unique(df_filtered$gene[df_filtered$bin==1])), parentterms=length(unique(df_filtered$parentterm[df_filtered$bin==1])), OR=exp(summary(model1)$coefficient[2,1]),lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4])) 
    },mc.cores=10))
  }))
    write.table(binned_gpsscore, paste0('Binned_by_sum_binsize0.3_',datafile,'_doe.txt'),sep='\t', quote=F, row.names=F)
})

#Analysis 5 - (data for fig 5) - Simulation of prioritization experiment, carry out in validation datasetset
set.seed(155)
genescoredataset=fread('All_genescoresum_across_drugs_Sider.txt', data.table=F)#GPS Sider data 
Simulate_MI_prioritization<-lapply(c(seq(0,2.1,0.3)), function(threshold){

 #step 1. Randomly select 1000 gene-phenotypes with a high GPS and get percentage of sampled gene-phenotype pairs with MI
     Dataset_genescores1_abovepercentilebin=genescoredataset %>% filter(genescoresum>threshold)  # subset of dataset with an GPS > threshold 
  gene_pt=Dataset_genescores1_abovepercentilebin %>% distinct(gene,parentterm)  # gene-parentterms with a high GPS score
  Randomgenept=gene_pt[sample(nrow(gene_pt), 1000, replace=TRUE), ] # randomly sample 1000 gene-phenotypes with replacement
  Randomgenept$ID=paste0(Randomgenept$gene,'_',Randomgenept$parentterm)
  Dataset_genescores1_abovepercentilebin$ID=paste0(Dataset_genescores1_abovepercentilebin$gene,'_',Dataset_genescores1_abovepercentilebin$parentterm)
  Dataset_genescores1_abovepercentilebin_random=subset(Dataset_genescores1_abovepercentilebin, (ID %in%Randomgenept$ID )) # get the data for the gene-phenotype sampling 
  RemoveDups <- function(df, column) {
    inds = sample(1:nrow(df))  
    df   = df[inds, ]
    dups = duplicated(df[, column])
    df   = df[!dups, ]
    inds = inds[!dups]
    df[sort(inds, index=T)$ix, ]
  }
  random_gene_pt_above_gps=RemoveDups(Dataset_genescores1_abovepercentilebin_random, 'ID') 
  random_gene_pt_above_gps_iteration=cbind.data.frame(Iteration=1, MI=table(random_gene_pt_above_gps$mi)[[2]], MI_percent=round(table(random_gene_pt_above_gps$mi)[[2]]/nrow(random_SRS_above2) *100,2), data='Above', threshold=threshold) # Report how many of these gene-phenotypes have an actual MI.

  #step2.Go through each gene-phenotype from the random selection, match on phenotype and then randomly select gene target from the entire dataset with out restricting using GPS.
  #Carry out 1000 iterations 
  Null_distribution<- mclapply(seq(1,nrow(Randomgenept),1), function(it){
    Randomgenept_match=do.call(rbind,lapply(c(1:nrow(Randomgenept)), function(i){
    generow<-Randomgenept[['parentterm']][i]
    DF<-Dataset_genescores1[which(Dataset_genescores1[['parentterm']]==generow),]
    inds = sample(nrow(DF),1)  
    DF   = DF[inds, ]  
    return(DF)
    }))
  Randomgenept_match_output=cbind(Iteration=it, MI=table(Randomgenept_match$mi)[2], MI_percent=round(table(Randomgenept_match$mi)[[2]]/nrow(Randomgenept_match) *100,2), data='Null', threshold=threshold)  #Report how many gene-phenotypes have an actual MI.
  return(Randomgenept_match_output)}, mc.cores=10)
  Null_distribution<-as.data.frame(do.call(rbind,Null_distribution))
  Allsimulation_data<-rbind(random_gene_pt_above_gps_iteration,Null_distribution)
})
Simulate_MI_prioritization1<-do.call(rbind,Simulate_MI_prioritization)
write.table(Simulate_MI_prioritization1, paste0('Simulation_highgps_matchgenes_mi_sider.txt'),sep='\t', quote=F, row.names=F)

#Analysis 6 - (data for supplementary fig 6a) - Stratify by clinical phase and run a logistic regression model with the GPS score as the predictor for each phase
#maxPhaseForIndication - clinical phase variable
genescoredataset=fread('All_genescoresum_opentargets.txt.gz',data.table=F) 
OT_dataset_full=fread(paste0('OT_drugdataset.txt'),data.table=F)
opentargets_allgenescore_clinicalphase=inner_join(OT_dataset_full,genescoredataset )

model1 =glm(as.formula(paste0('mi ~ genescoresum + number_gene_targets + category')), data=opentargets_allgenescore_clinicalphase,family = 'binomial')
mod_output<-rbind(cbind.data.frame(clinicalphase='All',No.MI=length(unique(opentargets_allgenescore_clinicalphase$parentterm[opentargets_allgenescore_clinicalphase$mi==1])), No.drugs=length(unique(opentargets_allgenescore_clinicalphase$drugname[opentargets_allgenescore_clinicalphase$mi==1])), OR=exp(summary(model2)$coefficient[2,1]),lowerCI=exp(summary(model2)$coefficient[2,1]-(1.96* summary(model2)$coefficient[2,2])),upperCI=exp(summary(model2)$coefficient[2,1]+(1.96* summary(model2)$coefficient[2,2])),P.val=summary(model2)$coefficient[2,4]))

strat_clinicalphase<-lapply(c(1,2,3,4), function(phase){
    opentargets_allgenescore_clinicalphase_strat=opentargets_allgenescore_clinicalphase[which(opentargets_allgenescore_clinicalphase$maxPhaseForIndication==phase|opentargets_allgenescore_clinicalphase$mi==0),]
    miphase4=opentargets_allgenescore_clinicalphase_strat %>% distinct(parentterm, mi,maxPhaseForIndication) %>% filter(mi==1 & maxPhaseForIndication==phase)  
    drugsphase4=opentargets_allgenescore_clinicalphase_strat %>% distinct(drugname,parentterm, mi,maxPhaseForIndication) %>% filter(mi==1 & maxPhaseForIndication==phase) %>% distinct(drugname)
    model1_strat =glm(as.formula(paste0('mi ~ genescoresum +number_gene_targets + category')), data=opentargets_allgenescore_clinicalphase_strat,family = 'binomial')
    mod_output_stratified<-rbind(cbind.data.frame(clinicalphase=phase,No.MI=length(unique(miphase4$parentterm)), No.drugs=nrow(drugsphase4), OR=exp(summary(model1_strat)$coefficient[2,1]),lowerCI=exp(summary(model1_strat)$coefficient[2,1]-(1.96* summary(model1_strat)$coefficient[2,2])),upperCI=exp(summary(model1_strat)$coefficient[2,1]+(1.96* summary(model1_strat)$coefficient[2,2])),P.val=summary(model1_strat)$coefficient[2,4]))
  return(mod_output_stratified)
  })

strat_clinicalphase1<-do.call(rbind,strat_clinicalphase)
strat_clinicalphase_and_all<-rbind(mod_output,strat_clinicalphase1)
write.table(strat_clinicalphase_and_all, paste0('Stratified_clinicalphase_genescoresumregression.txt'), sep='\t', quote=F, row.names=F)

#Analysis 7 - (data for supplementary fig 6b) - fold enrichment of drug indications with support from a high genetic priority score (0.9,1.5, 2.1) compared to those without in each clinical phase divided by the total sum observed in phase I.

genescoredataset=fread('All_genescoresum_opentargets.txt', data.table=F) #GPS OT data 
OT_dataset_full=fread(paste0('OT_drugdataset.txt'),data.table=F)

phasedatacounts_cutoff<-do.call(cbind,lapply(c(0.9,1.5,2.1), function(cutoff){

	#Subset drug dataset with scores with GPS ==0 . Count number of MI in each phase
Dataset_genescore_OR_below=Dataset_genescores1[which(Dataset_genescores1$genescoresum == 0),]
Dataset_genescore_OR_below_phase=right_join(OT_dataset_full,Dataset_genescore_OR_below) %>% distinct()
tablecounts_below=data.frame(GPS_binb_cutoff=cutoff,cbind(No.mi=nrow(Dataset_genescore_OR_below_phase),rbind(table(opentargets_mi1_score_allscores_below$OT_maxPhaseForIndication_max))))
#Subset drug dataset with scores above GPS that gives an OR 4. Count number of MI in each phase
Dataset_genescore_OR_above=Dataset_genescores1[which( Dataset_genescores1$genescoresum>=cutoff),]
Dataset_genescore_OR_above_phase=right_join(OT_dataset_full,Dataset_genescore_OR_above) %>% distinct()
tablecounts_above=data.frame(GPS_bin_cutoff=cutoff,cbind(No.mi=nrow(Dataset_genescore_OR_above_phase),rbind(table(Dataset_genescore_OR_above_phase$OT_maxPhaseForIndication_max))))
#MI counts for above and below threshold
tablecounts<-rbind(tablecounts_below,tablecounts_above)

#divide counts for phase II/III/IV by phase I counts for above and below GPS cutoff
phasecounts<-apply(tablecounts,1,function(cutoff){
      cutoff[2:7]=lapply(cutoff[2:7], function(X) as.numeric(X))
        phase2=cutoff[[5]]/(cutoff[[4]])
        phase3=cutoff[[6]]/(cutoff[[4]])
        phase4=cutoff[[7]]/(cutoff[[4]])
        counts_phase=cbind(cutoff$Predictors,phase2,phase3,phase4)
        counts_phase[2:4]=lapply(counts_phase[2:4], function(X) round(as.numeric(X),3))
        counts_phase1=as.data.frame(counts_phase)
        colnames(counts_phase1)=c('GPS_bin_cutoff','Phase II','Phase III', 'Phase IV')
        return(counts_phase1)
  })  
  phasecounts1<-do.call(rbind,phasecounts)
  phasecounts1<-data.frame(t(phasecounts1))
  colnames(phasecounts1)= phasecounts1[1,]
  phasecounts1=phasecounts1[-1,]
  phasecounts1[1:ncol(phasecounts1)]=lapply(phasecounts1[1:ncol(phasecounts1)], function(X) round(as.numeric(X),3))
  #fold difference
  phasecounts1$Folddif=phasecounts1[[2]]/phasecounts1[[1]] 
  phasecounts1$phase=rownames(phasecounts1)
}))

write.table(phasedatacounts_cutoff,paste0('Fold_enrichment_phasescomparedtophase0.txt'),sep='\t', quote=F, row.names=F)

#Analysis 8 - (data for supplementary table 4) 

OT_dataset_full=fread(paste0('OT_drugdataset.txt'),data.table=F)#full open targets dataset with clinical trial phase
### regression model
getModelMI <- function(predictor,data){
    model <-glm(as.formula(paste0('mi ~', paste(predictor, number_gene_targets, categories,sep = " + "))), data=data,family = 'binomial')
	  return(rbind(summary(model)$coefficient[predictor, c(1, 2, 4)]))
  }
#loop through each predictor, shuffle MI, apply model
phenotypes_model<- do.call(rbind,lapply(geneticpredictors, function(predictor) {
  OT_dataset_predictor=OT_dataset_full[,grepl(paste0('drugname|\\bgene\\b|category|\\mi\\b|parentterm|',predictor), colnames(OT_dataset_full))]
  OT_dataset_predictor$category=as.character(OT_dataset_predictor$category)

  resUnivariateRandomMI <- mclapply(1:10000, function(i) {
    cat(paste(i, "\n"))
    dataR <- OT_dataset_predictor
    dataR[("mi")] <- apply(dataR[("mi")], 2, sample)
    return(cbind.data.frame(Replicate = i, getModelMI(predictor = predictor,data = dataR)))
    }, mc.cores = 10)

  resUnivariateRandomMI1<-do.call(rbind,resUnivariateRandomMI)
  resUnivariateRandomResults_MI <- cbind.data.frame(Predictor = predictor,
  `Average estimate` = Reduce("+", lapply(resUnivariateRandomMI, function(res) res[, 2])) / length(resUnivariateRandomMI),
  `Average OR` = Reduce("+", lapply(resUnivariateRandomMI, function(res) exp(res[, 2]))) / length(resUnivariateRandomMI),
  `Average SD` = Reduce("+", lapply(resUnivariateRandomMI, function(res) res[, 3])) / length(resUnivariateRandomMI),
  `Percentage P<0.05` = Reduce("+", lapply(resUnivariateRandomMI, function(res) res[, 4] < 0.05)) / length(resUnivariateRandomMI) * 100,
  `Min P-value` = min(sapply(resUnivariateRandomMI, function(res) res[, 4]))
  )
}))
write.table(phenotypes_model, paste0('Univar_regression_opentargets_all_predictors_random_shuffleoutcome_10000_permutations.txt'), sep='\t', row.names=F, quote=F )


### False negative
#Analysis 9 - (data for supplementary table X) 
## Permute a percentage of each predictor to 0 to evaluate impact impact of false negatives on the GPS
set.seed(125)

False_negative<-lapply(c(0,1,2,5,10,20,30), function(percent_sample){
 lapply(c(seq(1,100,1)), function(perm) {
 # leave one out analysis
 percent_sample=as.numeric(percent_sample)
 print(paste(percent_sample, perm))
genescoredataset=fread('All_genescoresum_opentargets.txt.gz',data.table=F) 
genescoredataset_misclassified_percent=genescoredataset %>% select(drugname, gene, parentterm,category,mi)
  percent_missclas<-do.call(cbind,lapply(c(geneticpredictors), function(predname){
    tmp=which(genescoredataset[[predname]]!=0 & genescoredataset$mi==1)
    genescoredataset[[predname]][sample(tmp,round(percent_sample*length(tmp)))]<-0
    genescoredataset_misclassified_percent1<-cbind(genescoredataset_misclassified_percent, genescoredataset[predname])
      return(genescoredataset_misclassified_percent1)
    }))
  percent_missclas <- percent_missclas[, !duplicated(colnames(percent_missclas))]
  #sum all predictor values for each gene - phenotype row
  percent_missclas$genescoresum=rowSums(percent_missclas[, !names(percent_missclas) %in% c("gene", "parentterm","drugname","category","mi")])
  Genescores_beta_genetic_genept=percent_missclas %>% distinct(gene, parentterm,drugname,mi, category, genescoresum) 
  write.table(Genescores_beta_genetic_genept, paste0('All_genescoresum_across_drugs_Opentargets_permutation_',perm,'_',percent_sample,'.txt'), sep='\t', row.names=F, quote=F )
}))

    getModelMI_allgenes <- function(data){
        model <-glm(as.formula(paste0('mi ~ genescoresum + category + number_gene_targets')), data=data,family = 'binomial')
        return(rbind(cbind.data.frame(Estimate=(summary(model)$coefficient[2,1]), lowerCI=(summary(model)$coefficient[2,1]-(1.96* summary(model)$coefficient[2,2])),upperCI=(summary(model)$coefficient[2,1]+(1.96* summary(model)$coefficient[2,2])), P.val=summary(model)$coefficient[2,4]))) 
    }
   
   sampled<-do.call(rbind,lapply(c(0.01,0.05, 0.02, 0.1, 0.20,0.3), function(percent_sample){
    resUnivariateRandomMI <- mclapply(1:100, function(perm) {
    dataset_file<-fread(paste0(Revision_folder, 'All_genescoresum_across_drugs_Opentargets_permutation_',perm,'_',percent_sample,'.txt'), data.table=F)          
		return(cbind.data.frame(Replicate = perm,percent_sample= percent_sample, getModelMI_allgenes(data = dataset_file)))
  }, mc.cores = 10)
		resUnivariateRandomMI1<-do.call(rbind,resUnivariateRandomMI)
  Sampled_OR <- cbind.data.frame(percent_sample=percent_sample, 
	`Average OR` = Reduce("+", lapply(resUnivariateRandomMI, function(res) exp(res[, 6]))) / length(resUnivariateRandomMI),
	`Average lower CI` = Reduce("+", lapply(resUnivariateRandomMI, function(res) exp(res[, 7]))) / length(resUnivariateRandomMI),
	`Average upper CI` = Reduce("+", lapply(resUnivariateRandomMI, function(res) exp(res[, 8]))) / length(resUnivariateRandomMI),
	`Percentage P<0.05` = Reduce("+", lapply(resUnivariateRandomMI, function(res) res[, 9] < 0.05)) / length(resUnivariateRandomMI) * 100,
	`Min P-value` = min(sapply(resUnivariateRandomMI, function(res) res[, 9]))
)

}))

write.table(sampled, paste0(Revision_folder, 'False_negative_missclassification_opentargets_bymi.txt'), sep='\t', quote=F, row.names=F)



