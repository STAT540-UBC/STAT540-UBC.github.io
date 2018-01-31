setwd("~/ws/STAT540-instructors-only/homework/homework_winter_2017/01_quality_control/data")

experimentDes <- read_tsv("data/GSE60019_expdesign.data.txt", skip = 9)
expressionDat <- read_tsv("data/GSE60019_expmat.data.txt", skip = 6)

experimentDes <- experimentDes %>% select(id = Bioassay,
                                          sample = ExternalID,
                                          cell_type = cell.type,
                                          organism_part = organism.part,
                                          time_point = timepoint)

experimentDes %>% arrange(cell_type, time_point, organism_part) -> experimentDes

expressionDat[c("GeneSymbol", experimentDes$id)] -> expressionDat
names(expressionDat) <- c("GeneSymbol", experimentDes$sample)

expressionDat %>% group_by(GeneSymbol) %>% 
  summarise_all(mean) %>% 
  select(gene = GeneSymbol, everything()) %>% 
  arrange(gene) -> expressionDat


experimentDes$cell_type <- factor(experimentDes$cell_type, levels = c("surrounding_cell", "sensory_hair_cell"))

experimentDes$organism_part <- factor(experimentDes$organism_part, levels = c("epithelium_of_utricle", "sensory_epithelium_of_spiral_organ"))

experimentDes$time_point <- factor(experimentDes$time_point, levels = c("E16", "P0",  "P4",  "P7"))

experimentDes %>% select(sample, organism_part, cell_type, time_point) -> experimentDes



batchInfo <- tibble(sample = c("GSM1463883",
                               "GSM1463882",
                               "GSM1463887",
                               "GSM1463881",
                               "GSM1463872",
                               "GSM1463880",
                               "GSM1463875",
                               "GSM1463879",
                               "GSM1463876",
                               "GSM1463885",
                               "GSM1463871",
                               "GSM1463886",
                               "GSM1463878",
                               "GSM1463874",
                               "GSM1463884",
                               "GSM1463877",
                               "GSM1463888",
                               "GSM1463873"),
                    fastq_header = c("@SRR1534788.1.1 HWI-EAS00214_0040_FC226b:3:1:1115:13251", 
                                     "@SRR1534787.1.1 HWI-EAS00214_0030_FC00202:5:1:1366:996",
                                     "@SRR1534793.1.1 HWI-EAS00214_0030_FC00202:6:1:1096:985",
                                     "@SRR1534785.1.1 HWI-EAS00184_0044_FC199x3:7:1:15556:976",
                                     "@SRR1534772.1.1 HWI-EAS00214_0037_FC192:2:1:1083:19639",
                                     "@SRR1534784.1.1 HWI-ST363_0168:8:1101:1208:2082", 
                                     "@SRR1534777.1.1 HWI-ST363_0143:5:1101:1230:1984",
                                     "@SRR1534783.1.1 HWI-ST363_0144:8:1101:1217:2042", 
                                     "@SRR1534779.1.1 HWI-ST363_0144:7:1101:1219:2087",
                                     "@SRR1534790.1.1 HWI-EAS00184_0044_FC199x3:6:1:1557:1044",
                                     "@SRR1534770.1.1 HWI-EAS00184_0044_FC199x3:5:1:1424:1035",
                                     "@SRR1534792.1.1 HWI-EAS00214_0037_FC192:1:1:1092:6300",
                                     "@SRR1534781.1.1 HWI-ST363_0143:6:1101:1430:1972",
                                     "@SRR1534775.1.1 HWI-EAS00214_0037_FC192:3:7:1087:18479",
                                     "@SRR1534789.1.1 HWI-EAS00214_0040_FC:8:1:7957:1059",
                                     "@SRR1534780.1.1 HWI-ST363_0168:5:1101:1080:2090",
                                     "@SRR1534794.1.1 HWI-EAS00184_0040_FC196:2:1:1965:1068",
                                     "@SRR1534773.1.1 HWI-EAS00184_0044_FC199x3:8:1:1277:1066"))


batchInfo$batch <- batchInfo$fastq_header %>% strsplit(":") %>% 
  lapply(function(currHeader) { currHeader[1] %>% str_extract("HWI-.*(?=_0[0-9]{3})") }) %>% unlist() %>% 
  factor()

batchInfo <- batchInfo %>% select(sample, batch)

experimentDes %>% left_join(batchInfo) %>% saveRDS("data/gse60019_experiment_design.RDS")

for (currSample in experimentDes$sample) {
  expressionDat[[currSample]] <- 2^expressionDat[[currSample]]
}



# expressMatrix <- readRDS("gse60019_expression_matrix.RDS")
expressionDat <- expressionDat %>% filter(!is.na(gene))



expressionDat %>% saveRDS("gse60019_expression_matrix.RDS")

