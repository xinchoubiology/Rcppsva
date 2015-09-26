#' get TCGA-Assembler path
#' @title getTCGAssemble
#' @description TCGA Assembler is a two-R-script module for TCGA data indexing. \code{getTCGAssemble}
#'              is used for TCGA Assembler environment preparation and directionary file retrieval
#' @param ...
#' @return string TCGA traverse result path
#' @export
getTCGAssemble <- function(...) {
  source(system.file("extdata", "Module_A.r", package = "Rcppsva"))
  source(system.file("extdata", "Module_B.r", package = "Rcppsva"))
  DirectorytraverseTCGA <- system.file("extdata", "DirectoryTraverseResult.rda", package = "Rcppsva")
  DirectorytraverseTCGA
}

#' get 450K url
#' @title get450Kurl
#' @description \code{get450Kurl} function is used for specific TCGA archive meta-database file(.sdrf)
#'              downloading
#' @param type cancer type of TCGA 450K database
#' @param save save directionary
#' @param verbose FALSE(default)
#' @return absolute path of tcga data archive for specific cancer type
#' @importFrom stringr str_locate_all
#' @export
get450Kurl <- function(type = NULL, save = ".", verbose = FALSE, ...) {
  options(warn = -1)
  if(is.null(type)){
    stop("cancer type is necessary by TCGA search ...")
  }
  if(verbose){
    cat(sprintf("  downloading 450K methylation data's meta for cancer %s", type))
  }
  load(getTCGAssemble())
  specif_id  <- grep(pattern = toupper(paste("/", type, "/cgcc/jhu-usc\\.edu/", "HumanMethylation450", "/", sep = "")), x = upper_file_url, ignore.case = FALSE)
  cancer_idx <- specif_id[grepEnd(pattern = toupper("\\.sdrf\\.txt"), x = upper_file_url[specif_id], ignore.case = FALSE)]
  URL        <- GetNewestURL(AllURL = file_url[cancer_idx])
  sep_id     <- str_locate_all(string = URL, pattern = "/")[[1]][,1]
  sdrfname   <- paste(save, substr(URL, tail(sep_id,1)+1, str_length(URL)), sep="/")
  downloadFile(URL, sdrfname)
  sdrfname
}

#' parse sample sheet and minfi compatible
#' @title get450Ksheet
#' @description Rcppsva is compatible with minif package, as the result \code{get450Ksheet} is a function
#'              for samplesheet.csv preparation for user's methylation 450K microarray raw data
#' @param save sample sheet file save directory
#' @param sdrf sdrf directionary
#' @return sample table
#' @export
get450Ksheet <- function(save = ".", sdrf = NULL, verbose = FALSE) {
  
  write.sheet <- function(file, ctx) {
    header      <- data.frame(matrix(rep('', 28), 4, 7), stringsAsFactors=F)
    header[,1]  <- c("Investigator Name", "Project Name", "Experiment Name", "Date")
    header[4,2] <- as.character(Sys.Date())
    suppressWarnings(
      write.table(c("[Header]"), file = file, sep=",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    )
    suppressWarnings(
      write.table(header, file = file, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    )
    suppressWarnings(
      write.table(c("\n[Data]"), file = file, append = TRUE, sep=",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    )
    suppressWarnings(
      write.table(ctx, file = file, append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    )
  }
  
  as.age <- function(str) {
    str <- as.character(str)
    str[str == '[Not Available]'] <- NaN
    abs(as.numeric(str) / 365.25)
  }
  
  read.clinic <- function(cancer = NULL, clinical = c("patient", "drug", "follow_up")) {
    traversal <- getTCGAssemble()
    load(traversal)
    
    clinic.path <- paste0(save, "/nationwidechildrens.org_clinical_patient_", tolower(cancer), ".txt")
    DownloadClinicalData(traverseResultFile = traversal, saveFolderName = save, cancerType = cancer, clinicalDataType = clinical)
    
    clinics    <- read.csv(clinic.path, header  = TRUE, sep ="\t", skip = 1)
    barcode    <- as.character(clinics[,"bcr_patient_barcode"])
    clinics    <- data.frame(age  = as.age(clinics[,"days_to_birth"]),
                             sex  = as.numeric(clinics[,"gender"] == "MALE"))
    rownames(clinics) <- barcode
    na.omit(clinics)
  }
  
  if (is.null(sdrf)) {
    cat("sdrf file is not found, function will rebuild it")
    sdrf <- get450Kurl(type = "BRCA", save = save)
  }
  
  sample_sheet <- paste0(save, "/SampleSheet.csv")
  if (file.exists(sample_sheet)) {
    file.remove(sample_sheet)
  }
  file.create(sample_sheet)
  if (verbose) {
    message("Loading data from ", save, "......")
  }
  map   <- read.csv(paste0(save, "/",list.files(save)[grep("MAP", list.files(save), ignore.case = TRUE)]),
                    header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, row.names = 1)
  
  idat  <- list.files(paste0(save, "/DNA_Methylation/JHU_USC__HumanMethylation450/Level_1"))[grep("\\Grn.idat", list.files(paste0(save, "/DNA_Methylation/JHU_USC__HumanMethylation450/Level_1")), ignore.case = TRUE)]
  sdat  <- read.csv(sdrf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  slide <- sdat[[15]]  # 'Comment..TCGA.Archive.Name.'
  name  <- map[idat,]
  
  vars  <- read.clinic(cancer = "BRCA", clinical = "patient")
  keep.idx <- which(sapply(str_split(name, "-"), function(x){paste0(x[1:3], collapse = "-")}) %in% rownames(vars))
  name  <- name[keep.idx]
  idat  <- idat[keep.idx]
  vars  <- vars[sapply(name, function(x){paste0(str_split(x, "-")[[1]][1:3], collapse = "-")}, USE.NAMES = FALSE),]
  
  group <- sapply(str_split(name, pattern = "-"), function(xx) {
    if (as.numeric(str_sub(xx[4], start = 0, end = 2)) <= 9) {
      "C" ## cancer
    } else {
      "N"
    }
  })
  
  pData <- data.frame(
    Sample_Name  = sapply(str_split(name, pattern = ","), `[`, 1),
    Sample_Well  = rep('', length(name)),
    Sample_Plate = sapply(str_split(name, pattern = "-"), `[`, 6),
    Sample_Group = group,
    Pool_ID      = rep('', length(name)),
    Sentrix_ID   = sapply(str_split(idat, pattern = "_"), `[`, 1),
    Sentrix_Position = sapply(str_split(idat, pattern = "_"), `[`, 2),
    person       = sapply(str_split(name, pattern = "-"), `[`, 3),
    age          = vars$age,
    sex          = vars$sex,
    status       = group,
    stringsAsFactors = FALSE
  )
  rownames(pData) <- NULL
  # remove single sample batch
  batches <- table(pData$Sentrix_ID)
  batches <- names(batches)[batches > 1]
  pData   <- pData[pData$Sentrix_ID %in% batches, ]
  
  # paired batch filter
  pList   <- get450Kzip(pData)
  pData   <- pData[pData$Sentrix_ID %in% pList$Sentrix_ID, ]
  write.sheet(sample_sheet, pData)
  
  pData
}

#' when the ComBat analysis has been completed, we only needed to get the paired
#' methylation data for future research
#' @title get450Kzip
#' @param pdata Level 1 methylation 450K data matrix for pair extraction
#' @return matrix for paired methylation data
#' @export
get450Kzip <- function(pdata = NULL) {
  # get C vs N pair
  pList <- split(pdata, pdata$person)
  for (i in 1:length(pList)) {
    if (all(c("C", "N") %in% pList[[i]][,4])) {
      if (nrow(pList[[i]]) > 2) {
        normal <- pList[[i]][pList[[i]][,4] == "N",]
        cancer <- pList[[i]][pList[[i]][,4] == "C",]
        pList[[i]] <- list()
        for (x in 1 : nrow(normal)) {
          for (y in 1 : nrow(cancer)) {
            if (normal[x, 3] == cancer[y, 3] || normal[x, 6] == cancer[y, 6] || normal[x, 7] == cancer[y, 7]) {
              pList[[i]] <- rbind(normal[x,], cancer[y,])
            }
          }
        }
      }
    }
  }
  
  pData <- do.call(rbind,
                   pList[sapply(1:length(pList),
                                function(ix) {
                                  nrow(pList[[ix]]) == 2 && all(c("C","N") %in% pList[[ix]][,4])
                                })
                         ]
  )
  rownames(pData) <- NULL
  pData
}


