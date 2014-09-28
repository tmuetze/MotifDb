# MotifDb/inst/scripes/import/HOCOMOCO/test.R
#------------------------------------------------------------------------------------------------------------------------
library (RUnit)
#------------------------------------------------------------------------------------------------------------------------
source("import.R")
#------------------------------------------------------------------------------------------------------------------------
run.tests = function (dataDir)
{
  test.parsePwm()
  matrices.raw <- test.readRawMatrices (dataDir="./")
  matrices <- test.extractMatrices(matrices.raw)
  
  tbl.md <- test.createMetadataTable(dataDir="./", matrices, "md-raw.tsv")
  matrices <- test.normalizeMatrices(matrices)
  matrices.renamed <- test.renameMatrices(matrices, tbl.md)
  
} # run.tests
#------------------------------------------------------------------------------------------------------------------------
# create representative text here, a list of 5 lines, and make sure parsePWM can translate it into
# a one-element list, a 4 x 6 matrix (4 nucleotides by 6 positions) appropriately named.
test.parsePwm <- function()
{
  print("--- test.parsePwm")
  # a sample matrix, obtained as if from
  # all.lines <- scan (filename, what=character(0), sep='\n', quiet=TRUE)
  # text <- all.lines[1:5]
  text <-  c("> AHR_si",
             "40.51343240527031 10.877470982533044 21.7165707818416 2.5465132509466635 0.0 3.441039751299748 0.0 0.0 43.07922333291745", 
             "18.259112547756697 11.870876719950774 43.883079837598544 1.3171620263517245 150.35847450464382 0.7902972158110341 3.441039751299748 0.0 66.87558226865211", 
             "56.41253757072521 34.66312982331297 20.706746561638717 145.8637051322628 1.4927836298652875 149.37613720253387 0.7024864140542533 153.95871737667187 16.159862546986584", 
             "38.77363485291994 96.54723985087516 67.6523201955933 4.231336967110781 2.1074592421627525 0.3512432070271259 149.81519121131782 0.0 27.844049228115868") 
  pwm <- parsePwm(text)
  checkEquals(names(pwm), c("title", "matrix"))
  checkEquals(pwm$title, "> AHR_si")
  m <- pwm$matrix
  checkEquals(dim(m),c(4, 9))
  
  #tests would not turn true for an egality check for the last "checkTrue" case,
  #thus I used the closest bounds I could get using Rs colSum function
  checkTrue(all(as.numeric(colSums(m))>153.958717376671 & as.numeric(colSums(m)) < 153.958717376673))
  
} # test.parsePwm
#------------------------------------------------------------------------------------------------------------------------
test.readRawMatrices = function (dataDir)
{
  print ('--- test.readRawMatrices')
  list.pwms = readRawMatrices (dataDir)
  checkEquals (length (list.pwms), 426) #426 matrices in HOCOMOCO
  checkEquals (names (list.pwms [[1]]), c ("title", "matrix"))
  checkEquals (rownames (list.pwms[[1]]$matrix),  c ("A", "C", "G", "T"))
  invisible (list.pwms)
  
} # test.readRawMatrices
#------------------------------------------------------------------------------------------------------------------------
test.extractMatrices = function (matrices.raw)
{
  print ('--- test.extractMatrices')
  
  matrices.fixed <- extractMatrices(matrices.raw)
  checkEquals(length(matrices.fixed), 426)
  # checks the first 6 matrix names
  checkEquals(head(names(matrices.fixed)),
              c("AHR_si", "AIRE_f2", "ALX1_si", "ANDR_do", "AP2A_f2", "AP2B_f1"))
  
  invisible (matrices.fixed)
  
} # test.extractMatrices
#------------------------------------------------------------------------------------------------------------------------
test.normalizeMatrices = function (matrices)
{
  print ('--- test.normalizeMatrices')
  
  matrices.fixed <- normalizeMatrices(matrices)
  checkEquals(length(matrices.fixed), 426)
  
  checkEquals(head(names(matrices.fixed)),
              c("AHR_si", "AIRE_f2", "ALX1_si", "ANDR_do", "AP2A_f2", "AP2B_f1"))
  
  # make sure all columns in all matrices sum to 1.0
  checkTrue (all(sapply(matrices.fixed,
                        function(m)all(abs(colSums(m)-1.0) < 1e-10))))
  
  invisible (matrices.fixed)
  
} # test.extracNormalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
test.createMetadataTable = function (dataDir, matrices, raw.metadata.filename)
{
  print ('--- test.createMetadataTable')
  
  # try it with just the first matrix
  tbl.md = createMetadataTable (dataDir, matrices, raw.metadata.filename)
  
  checkEquals (dim (tbl.md), c (length(matrices), 15))
  checkEquals (colnames (tbl.md), c ("providerName", "providerId", "dataSource", "geneSymbol", "geneId", "geneIdType", 
                                     "proteinId", "proteinIdType", "organism", "sequenceCount", "bindingSequence",
                                     "bindingDomain", "tfFamily", "experimentType", "pubmedID"))
  
  with(tbl.md[1,], checkEquals(providerName,   "AHR_si"),
       checkEquals(providerId,     "AHR_si"),
       checkEquals(dataSource,     "AHR_si HOCOMOCOv9_AD_PLAINTEXT_H_PWM_hg19"),
       checkEquals(geneSymbol,     "AHR"),
       checkEquals(geneId,         "9606"),
       checkEquals(geneIdType,     "ENTREZ"),
       checkEquals(proteinId,      "P35869"),
       checkEquals(proteinIdType,  "uniprot"),
       checkEquals(organism,       "Hsapiens"),
       checkEquals(sequenceCount,  153.95872),
       checkEquals(experimentType, "low- and high-throughput methods"),
       checkEquals(pubmedID,       "23175603"),
       checkTrue(is.na(tfFamily)),
       checkTrue(is.na(bindingSequence)),
       checkTrue(is.na(bindingDomain)))
  
  invisible (tbl.md)
  
} # test.createMetadataTable
#------------------------------------------------------------------------------------------------------------------------
test.renameMatrices = function (matrices, tbl.md)
{
  print("--- test.renameMatrices")
  
  # try it with just the first two matrices
  matrix.pair <- matrices[1:2]
  tbl.pair <- tbl.md[1:2,]
  matrix.pair.renamed <- renameMatrices (matrix.pair, tbl.pair)
  checkEquals (names (matrix.pair.renamed), 
               c("Hsapiens-HOCOMOCOv9_AD_PLAINTEXT_H_PWM_hg19-AHR_si",
                 "Hsapiens-HOCOMOCOv9_AD_PLAINTEXT_H_PWM_hg19-AIRE_f2"))
  
} # test.renameMatrices
#------------------------------------------------------------------------------------------------------------------------
test.normalizeMatrices = function (matrices)
{
  print ('--- test.normalizeMatrices')
  
  colsums = as.integer (sapply (matrices, function (mtx) as.integer (mean (round (colSums (mtx))))))
  checkTrue (all (colsums > 1))
  
  matrices.norm = normalizeMatrices (matrices)
  
  colsums = as.integer (sapply (matrices.norm, function (mtx) as.integer (mean (round (colSums (mtx))))))
  checkTrue (all (colsums == 1))
  
  invisible (matrices.norm)
  
} # test.normalizeMatrices
#------------------------------------------------------------------------------------------------------------------------
