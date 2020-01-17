require(CCprofiler)
require(reshape2)
require(data.table)
require('R.utils')
require(reshape2)
require(tidyr)
require(dplyr)

# need to setwd() to correct file path

df<- fread(input = "feature_alignment.csv",data.table = T,showProgress = T)



# reformat OSW
df$filename <- gsub(pattern = "\\.mzXML\\.gz",replacement = "",x = df$filename)
df$ProteinName <- paste("1/",df$ProteinName,sep = "")
df$ProteinName  <- gsub(pattern = "1/DECOY_",replacement = "DECOY_1/",df$ProteinName)


# create fraction annotation
Annotation_txt<- data.table(filename = unique(df$filename), fraction_number = 1:61)

# import file
pepTraces <- importFromOpenSWATH(data = df, annotation_table = Annotation_txt, rm_requantified = TRUE,MS1Quant = F,verbose = T,rm_decoys = T)


# preprocessing
pepTraces_cons <- filterConsecutiveIdStretches(traces = pepTraces, 
                                               min_stretch_length = 3)

pepTraces_cons_sib <- filterBySibPepCorr(traces = pepTraces_cons,
                                         fdr_cutoff = NULL, 
                                         absolute_spcCutoff = 0.2, 
                                         plot = TRUE)

# quantify proteins
protTraces <- proteinQuantification(pepTraces_cons_sib, 
                                    topN = 2,
                                    keep_less = FALSE,
                                    rm_decoys = TRUE)


# complex centric analysis
complexHypotheses <- corumComplexHypotheses
binaryHypotheses <- generateBinaryNetwork(complexHypotheses)
pathLength <- calculatePathlength(binaryHypotheses)

corumTargetsPlusDecoys <- generateComplexDecoys(target_hypotheses=complexHypotheses,
                                                dist_info=pathLength,
                                                min_distance = 2,
                                                append=TRUE)



complexFeatures <- findComplexFeatures(traces=protTraces,
                                       complex_hypothesis = corumTargetsPlusDecoys)


complexFeaturesScored <- calculateCoelutionScore(complexFeatures)
qvalueComplexFeaturesScored <- calculateQvalue(complexFeaturesScored)
head(qvalueComplexFeaturesScored, n = 2)
qvalueComplexFeaturesScoredStats <- qvaluePositivesPlot(qvalueComplexFeaturesScored)
complexFeaturesFiltered <- subset(qvalueComplexFeaturesScored, qvalue <= 0.05)
summarizeFeatures(complexFeaturesFiltered)


# plot 26S proteasome
plotFeatures(feature_table = complexFeaturesFiltered,
             traces = protTraces,
             feature_id = "193",
             annotation_label="Entry_name",
             calibration = calibration,
             peak_area = TRUE)


