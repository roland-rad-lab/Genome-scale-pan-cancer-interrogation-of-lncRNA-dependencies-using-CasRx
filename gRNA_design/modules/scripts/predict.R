prediction <- function(arg1, arg2, arg3, arg4, arg5, directory, prediction_dir){
  #!/usr/bin/env Rscript
  #arg1 is sequence, arg2 is './data/Cas13designGuidePredictorInput.csv', arg3 is for output plot, arg4 is for the name, arg5 is the guide length
  ####################load_dependencies################################################## 
  ##### load dependencies
  
  # set this path relative to your working directory
  #dir='./'
  #dir.create(arg4)
  #setwd(arg4)
  dir <- directory
  #dir <- paste(directory, arg4, sep = "/")
  
  
  
  # Check to see if packages are installed. Install them if they are not, then load them into the R session.
  # Installation of BiocManager requires R v3.6 or higher, else please use BiocLite to install Biostrings
  # All prediction have been done using R v3.5.1
  libpath = paste0(dir,"scripts/R/")
  if (!dir.exists(paths = libpath)) {
    dir.create(path = libpath, recursive = TRUE)
  }
  .libPaths(new=libpath)
  check.packages <- function(pkg){
    for (p in pkg) {
      if (!requireNamespace(p, quietly = TRUE)) {
        if(p == "Biostrings"){
          suppressMessages(BiocManager::install(p , dependencies = TRUE))
        }else{
          suppressMessages(install.packages(p, repos = "http://cran.us.r-project.org"))
        }
      }
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
  
  packages = c("BiocManager","Biostrings","stringr","randomForest","ggplot2","dplyr")
  
  check.packages(packages)
  
  
  #########################input_files############################################# 
  ##### input files
  #args = commandArgs(trailingOnly=TRUE)
  args = c("1", "2", "3", "4", "5")
  if (length(args)!=5) {
    stop("Exiting! Please provide a single file.fasta input: <file.fasta>, Random Forest Model input <Cas13designGuidePredictorInput.csv> and if you would like the results plotted <true> or <false>", call.=FALSE)
  } else if (length(args)==5) {
    print(arg1)
    print(paste("there are", length(args), "args"))
    FA <- DNAStringSet(arg1) # FA <- Biostrings::readDNAStringSet(filepath = "./data/test.fa", format = "fasta", use.names = T)
    names(FA) <- arg4
    print(length(FA))
    if (length(FA) == 1){
      name = gsub(">","",strsplit(names(FA), split = "\\|")[[1]][1])
      print(name)
      if(nchar(FA) < 80){
        warning("RNAplfold requires a minimum of 80nt length for the target site accessibility calculation")
      }
      
      if(nchar(FA) < 30){
        stop("Exiting! Input sequence to short. Please provide an input sequence >30nt to be able to assess the target site nucleotide context.")
      }
    }
    else{
      stop("Exiting. Please provide single-entry fasta files. It is faster to run multiple jobs on individual sequences than providing multiple sequences in one job.")
    }
    fullset = read.delim( arg2, header = T, sep = ",", stringsAsFactors = F)  # fullset = read.delim( "./data/Cas13designGuidePredictorInput.csv", header = T, sep = ",", stringsAsFactors = F) 
    fields = c('normCS', 'MFE','DR','Gquad','Log10_Unpaired','hybMFE_3.12','hybMFE_15.9' ,
               'NTdens_max_A','NTdens_max_C','NTdens_max_G','NTdens_max_T','NTdens_max_AT' , 'NTdens_max_GC',
               'NTdens_min_A' ,  'NTdens_min_C' ,  'NTdens_min_G' ,  'NTdens_min_T' , 'NTdens_min_AT' ,
               'pA','pC','pG',
               'pAA','pAC','pAG','pAT','pCA','pCC','pCG','pCT','pGA','pGC','pGG','pGT','pTA','pTC','pTG')
    
    if (arg3 %in% c("true","True","TRUE","T","treu","Treu","TREU","ture","Ture","TURE")){
      PLOT = TRUE
    }else{PLOT = FALSE}
  }
  
  
  
  
  
  ############################set_parapemeteres########################################## 
  # Set Parameters
  guideLength = arg5#23
  
  WindowOffset = 0 # Offset to start scoring guides. Some of the features will depend of the target sequence context upstream. Without the offset, NAs will be assigned for the respective feature
  GCmin = 0.1 # probability value between 0 and 1
  GCmax = 0.9 # probability value between 0 and 1
  # MinimumFreeEnergyCutoff = tbd
  # Tag = DNAString("aaac") # antiTag sequence as defined by Marraffini et al. 2018. As not being present for CasRx, it is disregarded here
  # antiTag <- reverseComplement(Tag) # As above. "GTTT" # reverse complement of Tag, this can be found, targetRNA sequences should not be 3' followed by GUUUCA with decreasing nt weight from G... -> ...A
  HomopolymerLength.T = 4
  HomopolymerLength.nonT = 5
  DirectRepeat = "aacccctaccaactggtcggggtttgaaac" # Used for crRNA MFE calculation using RNAfold Vienna package 2.4.10
  RNAfold = "/usr/local/bin/RNAfold" # Change if different
  RNAplfold = "/usr/local/bin/RNAplfold" # Change if different
  RNAhyb = paste0(dir,"scripts/RNAhyb.sh")
  
  # RNAhybrid 2.1.2 has to be installed and added to your PATH. Here I am assuming /usr/local/bin/RNAhybrid
  Sys.setenv("PATH" = paste(Sys.getenv('PATH'), '/usr/local/bin', sep = ':'))  
  
  
  
  
  
  
  ######################################################################
  # change nothing below here
  #######################subroutines############################################### 
  # subroutines
  
  # get position of match start
  extractPos <- function(y, END = "end"){
    coord <- strsplit( strsplit(y, split = ":")[[1]][2] ,split = "_")[[1]][1]
    if (END == "end" ){
      return(as.integer(strsplit(coord, split = "-")[[1]][2]))
    }
    else{
      return(as.integer(strsplit(coord, split = "-")[[1]][1]))
    }
  }
  
  #Gets the GC content of a guide given an XStringset
  GetGC <- function(y){
    return(sum(alphabetFrequency(y, baseOnly=TRUE, as.prob=TRUE )[c("C", "G")]))
  }
  GetConsecutiveBases <- function(y){
    consecutiveBases <- vector()
    consecutiveBases[1] <- longestConsecutive(as.character(y), "A")
    consecutiveBases[2] <- longestConsecutive(as.character(y), "C")
    consecutiveBases[3] <- longestConsecutive(as.character(y), "G")
    consecutiveBases[4] <- longestConsecutive(as.character(y), "T")
    names(consecutiveBases) <- c("A","C","G","T")
    return(consecutiveBases)
  }
  
  GetRawGuides <- function(x, GC.min = GCmin, GC.max = GCmax , HomopolymerLengthT=4, HomopolymerLengthnonT=5 ){
    # get all possible guide start coordinates (transcriptomic starts are guide ends)
    starts <- seq(from = 1, to = nchar(x) - guideLength + 1 , by = 1)
    # get all possible guide end coordinates (transcriptomic ends are guide starts)
    ends <- seq(from = guideLength, to = nchar(x) , by = 1)
    # get all possible guide sequences
    SubStrings <- Views( x , start=starts, end=ends)
    # get the guide GC content
    GC <-  sapply(SubStrings,GetGC)
    # Get length of Homopolymers (longest consecutive base stretch)
    ConsecutiveBases <-  do.call( rbind,lapply(SubStrings,GetConsecutiveBases))
    
    # filter
    GC.index <- which(GC > GC.min & GC < GC.max)
    ConsecutiveBases.index <- which(apply(ConsecutiveBases[,c("T","C","G")] , MARGIN = 1 , max) <= HomopolymerLengthnonT & ConsecutiveBases[,c("A")] <= HomopolymerLengthT)
    
    
    
    if (length(GC.index) > 0 & length(ConsecutiveBases.index) > 0){
      
      idx <- intersect(GC.index, ConsecutiveBases.index)
      
      if (length(idx) > 0){
        
        return(SubStrings[idx])
        
      }
      else{
        return("no intersect: GC and Homopolymers out of range")
      }
      
    }
    else{
      if ( length(GC.index) == 0 & length(ConsecutiveBases.index) > 0){
        return("GC out of range")
      }
      else if (  length(GC.index) > 0 & length(ConsecutiveBases.index) == 0 ){
        return("Homopolymers out of range")
      }
      else{
        return("GC and Homopolymers out of range")
      }
    }
    
  }
  
  
  
  Get_NT_density_Vector = function(fa = FA,NT = "G", WINDOW = 30){
    
    D = WINDOW
    
    ma <- rep(NA,ncol = width(fa))  
    
    
    if((D %% 2) == 0) {
      
      d = D/2
      
      for (p in 1:width(fa)){
        
        if ( (p-(d-1)) < 1 | (p+d) > width(fa) ){
          ma[p] <- NA
        }
        else{
          ma[p] <-  letterFrequency( subseq(fa, start=p-(d-1), end=p+d)  , letters = NT ,as.prob = T) 
        }
      }
    } 
    
    else {
      
      d = (D-1)/2
      for (p in 1:width(fa)){
        if ( (p-d) < 1 | (p+d) > width(fa) ){
          ma[p] <- NA
        }
        else{
          ma[p] <-  letterFrequency( subseq(fa, start=p-d, end=p+d)  , letters = NT ,as.prob = T) 
        }
      }
    }
    
    return(ma)
  }
  
  GetNTdensitities = function(fa = FA){
    NTs = c("A","C","G","T","AT","GC")
    ma.A = Get_NT_density_Vector(fa = fa, NT = NTs[1] , WINDOW = max[which(max$NT  == NTs[1]),"W"] )
    ma.C = Get_NT_density_Vector(fa = fa, NT = NTs[2] , WINDOW = max[which(max$NT  == NTs[2]),"W"] )
    ma.G = Get_NT_density_Vector(fa = fa, NT = NTs[3] , WINDOW = max[which(max$NT  == NTs[3]),"W"] )
    ma.T = Get_NT_density_Vector(fa = fa, NT = NTs[4] , WINDOW = max[which(max$NT  == NTs[4]),"W"] )
    ma.AT = Get_NT_density_Vector(fa = fa, NT = NTs[5] , WINDOW = max[which(max$NT == NTs[5]),"W"] )
    ma.GC = Get_NT_density_Vector(fa = fa, NT = NTs[6] , WINDOW = max[which(max$NT == NTs[6]),"W"] )
    
    mi.A = Get_NT_density_Vector(fa = fa, NT = NTs[1] , WINDOW = min[which(min$NT  == NTs[1]),"W"] )
    mi.C = Get_NT_density_Vector(fa = fa, NT = NTs[2] , WINDOW = min[which(min$NT  == NTs[2]),"W"] )
    mi.G = Get_NT_density_Vector(fa = fa, NT = NTs[3] , WINDOW = min[which(min$NT  == NTs[3]),"W"] )
    mi.T = Get_NT_density_Vector(fa = fa, NT = NTs[4] , WINDOW = min[which(min$NT  == NTs[4]),"W"] )
    mi.AT = Get_NT_density_Vector(fa = fa, NT = NTs[5] , WINDOW = min[which(min$NT == NTs[5]),"W"] )
    mi.GC = Get_NT_density_Vector(fa = fa, NT = NTs[6] , WINDOW = min[which(min$NT == NTs[6]),"W"] )
    
    L = list(ma.A,ma.C,ma.G,ma.T,ma.AT,ma.GC  ,  mi.A,mi.C,mi.G,mi.T,mi.AT,mi.GC)
    names(L) = c(paste0('max_',NTs),paste0('min_',NTs)) 
    return(L)
  }
  
  GetVal = function(j,vec=VEC, POINT = -11){
    
    # if j is NA, return vector of NAs
    if (is.na(j) == T){
      return(NA)
    }
    
    else{
      if ((j + POINT) < 1){
        return(NA)
      }
      else{
        return(vec[j + POINT])
      }
    }
  }
  
  GetNTpointdensities = function(DAT = CD71 , Vec.list = NTdensitities.cd71){
    
    ma.A = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["max_A"]] , POINT = max[which(max$NT  == "A"),"P"])
    ma.C = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["max_C"]] , POINT = max[which(max$NT  == "C"),"P"])
    ma.G = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["max_G"]] , POINT = max[which(max$NT  == "G"),"P"])
    ma.T = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["max_T"]] , POINT = max[which(max$NT  == "T"),"P"])
    ma.AT = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["max_AT"]] , POINT = max[which(max$NT  == "AT"),"P"])
    ma.GC = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["max_GC"]] , POINT = max[which(max$NT  == "GC"),"P"])
    
    mi.A = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["min_A"]] , POINT = min[which(min$NT  == "A"),"P"])
    mi.C = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["min_C"]] , POINT = min[which(min$NT  == "C"),"P"])
    mi.G = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["min_G"]] , POINT = min[which(min$NT  == "G"),"P"])
    mi.T = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["min_T"]] , POINT = min[which(min$NT  == "T"),"P"])
    mi.AT = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["min_AT"]] , POINT = min[which(min$NT  == "AT"),"P"])
    mi.GC = sapply( as.list(DAT$MatchPos) , FUN = GetVal ,vec=Vec.list[["min_GC"]] , POINT = min[which(min$NT  == "GC"),"P"])
    L = list(ma.A,ma.C,ma.G,ma.T,ma.AT,ma.GC    ,   mi.A,mi.C,mi.G,mi.T,mi.AT,mi.GC)
    names(L) = paste0( "NTdens_" , names(Vec.list))
    dens = do.call( cbind , L)
    out = cbind.data.frame( DAT , dens)
    return(out)
  }
  
  
  
  PrePareModelInput = function(x, FIELDS=fields){
    
    
    # The effect size value is normCS
    
    # remove incomplete entries
    x = x[complete.cases(x),FIELDS]
    
    # scale numeric values for training
    numeric = sapply(x, is.numeric) & colnames(x) != "normCS"  & colnames(x) != "Gquad"  & colnames(x) != "DR" & grepl("[A,C,G,T]_",colnames(x)) == FALSE #DO NOT SCALE response value to [0,1] interval
    
    tmpmean = apply(x[,numeric], 2, function(x) quantile(x, 0.05))
    tmpsd = apply(x[,numeric], 2, function(x) quantile(x, 0.95)-quantile(x, 0.05)) 
    x[,numeric] = scale(x[,numeric], center=tmpmean, scale=tmpsd) 
    x[,numeric][x[,numeric] > 1] = 1
    x[,numeric][x[,numeric] < 0] = 0
    
    # set seed for reproducibility
    set.seed(1234)
    #generate model  
    model = randomForest(x[,-1], x[,1], importance = FALSE, ntree = 2000)
    
    L = list(model, tmpmean, tmpsd)
    names(L) = c("model", "tmpmean", "tmpsd")
    return(L)
    
  }
  
  
  PredictGuideScores = function(x = Guide.df, FIELDS=fields[-1], Fit=Model, Mean = ModelInputMeans, SD = ModelInputSDs){
    
    # remove incomplete entries
    x = x[complete.cases(x),FIELDS]
    
    
    
    # scale Features
    numeric = sapply(x, is.numeric) & colnames(x) != "normCS"  & colnames(x) != "Gquad"  & colnames(x) != "DR" & grepl("[A,C,G,T]_",colnames(x)) == FALSE #DO NOT SCALE response value to [0,1] interval
    x[,numeric] = scale(x[,numeric], center=Mean, scale=SD)
    
    x[,numeric][x[,numeric] > 1] = 1
    x[,numeric][x[,numeric] < 0] = 0
    
    return(predict(Fit, newdata = x))
    
  }
  
  standardizeScores = function(x, MIN = min(fullset$normCS) , MAX = max(fullset$normCS)){
    if(is.na(x)){
      return(NA)
    }
    else{
      (x - MIN) / (MAX - MIN)
    }
  }
  
  
  AddScores = function( df = Guide.df , SCORES = GuideScores){
    
    #match by guide name
    idx = match(names(SCORES) , rownames(df))
    
    if (all(names(GuideScores) == rownames(Guide.df)[idx])){
      
      # initiate GuideScore column
      df$GuideScores = NA
      # add scores given the positional index 
      df$GuideScores[idx] <- GuideScores
      # initiate rank column
      df$Rank = NA
      # Assign Rank
      # For guides that reside to close to the target's 5' end it may be that not all features are assigned. Thus, all guides with NA features will not be ranked.
      # The rank ranges from 0 to 1, with 1 being the highest rank
      df$Rank[is.na(df$GuideScores) == F] <- signif(1- (rank(-df$GuideScores[is.na(df$GuideScores) == F], na.last = T, ties.method = "first") /length(df$GuideScores[is.na(df$GuideScores) == F])),4) 
      # order by rank
      df = df[order(df$Rank, decreasing = T),]
      # return object
      return(df)
      
    }
    else{
      stop("Exiting! Guide names do not correspond.")
    }
  }
  
  
  AssignQuartiles = function(x,q = quantile(x)){
    quartiles = vector()
    for ( i in 1:length(x)){
      if (is.na(x[i]) == T){
        quartiles[i] = NA
      }
      else if (x[i] <= q["25%"]){
        quartiles[i] = 1
      }
      else if (x[i] > q["25%"] & x[i] <= q["50%"]){
        quartiles[i] = 2
      }
      else if (x[i] > q["50%"] & x[i] <= q["75%"]){
        quartiles[i] = 3
      }
      else if (x[i] > q["75%"]){
        quartiles[i] = 4
      }
      else{
        stop("Exiting. Value outside quartiles")
      }
    }
    return(quartiles) 
  }
  
  GetCDSregionBounderies = function(x){
    tmp = strsplit(x, split = "\\|")[[1]]
    r = tmp[grep("CDS",tmp)]
    if (length(r) == 0){
      return(NULL)
    }
    else{
      return(as.numeric(strsplit( strsplit(r , split = "\\:")[[1]][2] , split = "-")[[1]]))
    }
  }
  
  evalFOLD <- function(x){
    if(substr(x, 1, 24) == "((((((.(((....))).))))))"){
      return(1)
    }
    else{
      return(0)
    }
  }
  
  GetMFE = function(g,DR = DirectRepeat){
    crRNA = paste0(DR,g)
    # You may need to change the path to your RNAfold executable 
    cmd = paste0( "echo ",crRNA, " | ", RNAfold , " --gquad --noPS")
    output = system(cmd , intern = TRUE)
    mfe = as.numeric(gsub("\\)","",gsub("\\(","",strsplit(output[2], split = " ")[[1]][2])))
    gq = ifelse( grepl("\\+",output[2]) == TRUE, 1, 0)
    dr = evalFOLD(output[2])
    return(c(mfe,gq,dr))
  }
  
  ReadUnpairedPorbabilities = function(x){
    UnpairedProbabilities <- read.delim(x, sep="\t", skip = 1, row.names = "X.i.")[,1:50]
    colnames(UnpairedProbabilities) <- seq(1,50,1)
    UnpairedProbabilities=t(UnpairedProbabilities)
    return(UnpairedProbabilities)
  }
  
  Transform_RNAplfold_predictions <- function(x){
    
    
    ma <- matrix(NA,ncol = ncol(x), nrow = nrow(x) )  
    colnames(ma) <- colnames(x)
    rownames(ma) <- rownames(x)
    for (i in 1:nrow(x)){
      
      if((i %% 2) == 0) {
        
        d = (i/2)-1
        
        ma[i, 1:(ncol(x)-d) ]  <- x[i, (1+d):ncol(x) ]
        
      } 
      else {
        
        d = (i-1)/2
        
        ma[i, 1:(ncol(x)-d) ]  <- x[i, (1+d):ncol(x) ]
        
      }
    }
    return(ma)
  }
  
  
  SliceMatrix <- function(j,mat=MA, w=50){
    
    if (is.na(j) == T){
      tmp = matrix(NA, ncol = 2*w+1, nrow = nrow(mat))
      colnames(tmp) <- seq(-w,w,1)
      rownames(tmp) <- rownames(mat)
    }
    else{
      tmp = matrix(NA, ncol = length((j-w):(j+w)), nrow = nrow(mat))
      rownames(tmp) <- rownames(mat)
      
      if ((j-w) < 1){
        
        tmp[,( ((w+1)-j+1):ncol(tmp))] <-  mat[,1:(j+w)]
        
      }
      else if ((j+w) > ncol(mat)){
        
        tmp[,( 1 : ( (w+1) + (ncol(mat) - j)) )] <-  mat[,(j-w):(ncol(mat))]
        
      }
      else{
        
        tmp <- mat[,(j-w):(j+w)]
        
      }
      
      colnames(tmp) <- seq(-w,w,1)
      return(tmp)
    }
  }   
  
  GetValue <- function(x,d,w){ 
    if (is.null(nrow(x)) == T){
      return(NA)
    }
    else{
      return(x[d,w])
    }
  }
  
  
  GetUnpairedProb <- function( x = Guide.df , MA = log10(UnpairedProbabilities.tranformed), W=50, X=40 , Y=23){
    
    # slice the nucleotide density matrix for each guide to obtain a window of +/- the window size (default 50nt) centered on guide match position 1
    Slices = lapply( as.list(x$MatchPos) , FUN = SliceMatrix ,mat=MA, w=W)
    names(Slices) <- rownames(x)
    
    dens <- sapply(Slices , FUN = GetValue, d = Y, w = X)
    
    return(dens)
    
  }  
  
  
  GetTargetSiteAccessibility = function( dat = Guide.df, fa = FA ){
    
    #generate random string for tmp file that will be written to the hard drive to avoid collisions
    RanStr = paste0(sample( letters , size = 6 , replace = F), collapse = "")
    names(fa) = RanStr
    
    # writing the fa file back to hard drive as a tmp file. This is done to be independent of any naming issues
    writeXStringSet(fa, filepath = paste0('./',RanStr,'.fa'), append=FALSE, format="fasta")
    
    # You may need to change the path to your RNAplfold executable 
    cmd = paste0( "cat ", paste0('./',RanStr,'.fa')  , " | ", RNAplfold , " -L 40 -W 80 -u 50 ")
    output = system(cmd , intern = TRUE )
    
    UnpairedProbabilities = ReadUnpairedPorbabilities(x = paste0('./',RanStr,'_lunp') )
    UnpairedProbabilities.tranformed <- Transform_RNAplfold_predictions(UnpairedProbabilities)
    
    # As there was no clear pattern, this will only record the unpaired probability covering the entire guide match
    Log10_Unpaired <- GetUnpairedProb(x=dat , MA = log10(UnpairedProbabilities.tranformed), W=50, X=40 , Y=23)
    
    # clean up
    tmp = file.remove( paste0('./',RanStr,'.fa') ,paste0('./',RanStr,'_lunp'),  paste0('./',RanStr,'_dp.ps'))
    
    
    return(Log10_Unpaired)
  }
  
  GetRNAhybMFE = function( g ){
    cmd = paste0( "RNAhybrid -c -s 3utr_human ",g," ",reverseComplement(g))
    output = system(cmd , intern = TRUE)
    return(as.numeric(strsplit(output, split = ":")[[1]][5])) # Gets MFE
  }
  
  
  GetRNAhybMFE_bulk = function( dat = Guide.df , POS = 3 , WIDTH = 12 ){
    
    # transform to DNAstringSet
    GuideSeq = DNAStringSet(dat$GuideSeq)
    # extract guide
    g = as.character(subseq( x = GuideSeq , start = POS , width = WIDTH))
    # extract target
    t = as.character(reverseComplement(subseq( x = GuideSeq , start = POS , width = WIDTH)))
    # write tmp file to hard disk 
    tmp = cbind( g , t )
    RanStr = paste0(sample( letters , size = 6 , replace = F), collapse = "")
    write.table( tmp , file =  paste0(arg4, RanStr,'_hybMFE_',POS,'.',WIDTH,'.txt') , sep = "," , quote = F, col.names = F, row.names = F)
    # calculate RNA hyb MFE in bulk 
    cmd = paste0( "bash ",RNAhyb," ", paste0(arg4, RanStr,'_hybMFE_',POS,'.',WIDTH,'.txt') )
    output = system(cmd , intern = TRUE )
    # extract the MFE
    hybMFE = as.numeric(sapply( output , FUN=function(j){as.numeric(strsplit(j, split = ":")[[1]][5])})) # Gets MFE)
    # clean up
    tmp = file.remove( paste0(arg4, RanStr,'_hybMFE_',POS,'.',WIDTH,'.txt'))
    # return value
    return(hybMFE)
  }
  
  
  
  GetLetterProbs = function( x = Guide.df, S=1, E=23){
    
    
    G = DNAStringSet(x$GuideSeq)
    G.sub = subseq( G,start = S,end = E)
    G = G.sub
    
    A.prob <- letterFrequency( G , letters = c("A") ,as.prob = T)
    C.prob <- letterFrequency( G , letters = c("C") ,as.prob = T)
    G.prob <- letterFrequency( G , letters = c("G") ,as.prob = T)
    U.prob <- letterFrequency( G , letters = c("T") ,as.prob = T)
    
    AU.prob <- letterFrequency( G , letters = c("AT") ,as.prob = T)
    GC.prob <- letterFrequency( G , letters = c("GC") ,as.prob = T)
    
    diNucleotide.prob <- dinucleotideFrequency(G, as.prob = T)
    
    
    out <- cbind.data.frame(A.prob,C.prob,G.prob,U.prob,GC.prob,AU.prob,diNucleotide.prob)
    rownames(out) <- names(G)
    colnames(out) <- c("pA","pC","pG","pT",
                       "pG|pC", "pA|pT", 
                       "pAA","pAC","pAG","pAT",
                       "pCA","pCC","pCG","pCT",
                       "pGA","pGC","pGG","pGT",
                       "pTA","pTC","pTG","pTT")
    
    return(out)
    
  }
  
  
  
  ####
  
  # Plotting subroutines 
  GetCDSregionBounderies = function(x){
    tmp = strsplit(x, split = "\\|")[[1]]
    r = tmp[grep("CDS",tmp)]
    if (length(r) == 0){
      return(NULL)
    }
    else{
      return(as.numeric(strsplit( strsplit(r , split = "\\:")[[1]][2] , split = "-")[[1]]))
    }
  }
  
  
  
  gg_plot <- function(df=Guide.df,LEN=width(FA), NAME = name , CDS = cds_coord, f=0.2){
    
    
    df = filter(df, !is.na(quartiles))
    if(LEN < 500){b=25}else if(LEN >= 500 & LEN < 1000 ){b=50}else if(LEN >= 1000 & LEN < 2000 ){b=100}else if(LEN >= 2000 & LEN < 5000 ){b=500}else if(LEN >= 5000 & LEN < 10000 ){b=1000}else{b=4000}
    
    # # Plot predictions
    if (is.null(CDS) == T){
      mx = max(df$standardizedGuideScores, na.rm = T)
      mn = min(df$standardizedGuideScores, na.rm = T)
      su = mx-mn
      df$quartiles = factor(df$quartiles, levels = (unique(df$quartiles)))
      Ymax=1
      Ymin=0
      
      g=ggplot(df, aes(x=MatchPos, y=standardizedGuideScores)) +
        geom_point(shape=20,  aes(color = quartiles))  +
        geom_smooth(span = 0.025 , method = "loess", formula = y~x, colour = "#525252" ) +
        theme_classic() +
        scale_color_manual(values = rev(c('#ca0020','#f4a582','#92c5de','#0571b0'))) +
        ggtitle(NAME) +
        ylim(c(Ymin,Ymax)) + ylab("standardized guide score") + xlab("guide match position [nt]") +
        coord_fixed(ratio=(LEN/(Ymax-Ymin))*f) +
        scale_x_continuous(breaks = seq(0,LEN,b)) +
        theme(text = element_text(size=12)) 
      
    }
    else{
      cds=data.frame(matrix(CDS, ncol = 2))
      colnames(cds) = c("start","end")
      mx = max(df$standardizedGuideScores, na.rm = T)
      mn = min(df$standardizedGuideScores, na.rm = T)
      su = mx-mn
      df$quartiles = factor(df$quartiles, levels = (unique(df$quartiles)))
      Ymax=1
      Ymin=0
      
      g=ggplot(df, aes(x=MatchPos, y=standardizedGuideScores)) +
        geom_rect(data=cds, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=Ymin, ymax=Ymax), color="transparent", fill="#969696", alpha=0.2) +
        geom_point(shape=20,  aes(color = quartiles))  +
        geom_smooth(span = 0.025 , method = "loess", formula = y~x, colour = "#525252" ) +
        theme_classic() +
        scale_color_manual(values = rev(c('#ca0020','#f4a582','#92c5de','#0571b0'))) +
        ggtitle(NAME) +
        ylim(c(Ymin,Ymax)) + ylab("standardized guide score") + xlab("guide match position [nt]") +
        coord_fixed(ratio=(LEN/(Ymax-Ymin))*f) +
        scale_x_continuous(breaks = seq(0,LEN,b)) +
        theme(text = element_text(size=12)) +
        annotate(geom = "text", x = cds$start , y = 0.01, label = "CDS" , hjust = 0, vjust=0)
      
    }
    print(g)
  }
  
  
  
  ##########################main############################################ 
  # Guide design
  
  
  # shorten the windows by a user defined window edge offset
  if (WindowOffset > 0){
    FA <- subseq(x = FA,start = 0+WindowOffset, end = (nchar(FA)-WindowOffset))
  }
  
  cat( paste0("Getting guides for: ", name ,"\n"))
  
  # Get all sequences of a given guide length and filter by GC content and homopolymers
  cat( paste0( "Get raw CasRx guides sequences started on " , date() ,"\n"))
  # get raw guides sense
  rawGuides <- lapply(FA,GetRawGuides, GC.min = GCmin, GC.max = GCmax, HomopolymerLengthT=HomopolymerLength.T, HomopolymerLengthnonT=HomopolymerLength.nonT) 
  # get reverse complement
  rawGuides.rc <- reverseComplement(rawGuides[[1]])
  # transform to DNAStringSet
  rawGuides.dnaStringSet <- DNAStringSet(rawGuides.rc) 
  # Assign names
  names(rawGuides.dnaStringSet) <- paste0("crRNA", str_pad(seq(1,length(rawGuides.rc),1), nchar(length(rawGuides.rc)), pad = "0"),":",nchar(FA)-end(rawGuides.rc)+1, "-" ,nchar(FA)-start(rawGuides.rc)+1) 
  # Transform to data frame
  Guide.df = as.data.frame(rawGuides.dnaStringSet)
  colnames(Guide.df) = "GuideSeq"
  Guide.df$MatchPos = sapply( rownames(Guide.df) , extractPos)
  
  
  
  # Calculate minimum free enegergy for crRNA fold nt including the direct repeat sequence
  cat( paste0( "Calculating crRNA MFE started on " , date() ,"\n"))
  mfe.info = do.call( rbind , lapply(  Guide.df$GuideSeq , FUN = GetMFE, DR = DirectRepeat) )
  Guide.df$MFE = mfe.info[,1] 
  Guide.df$DR = mfe.info[,2] 
  Guide.df$Gquad = mfe.info[,3] 
  
  
  # Calculate target site accessibility
  cat( paste0( "Calculating target site accessibility started on " , date() ,"\n"))
  Guide.df$Log10_Unpaired <- GetTargetSiteAccessibility(dat = Guide.df, fa = FA)
  
  
  
  
  # Calculate guide::target hybridization energies
  # Takes about 2.8 seconds per 100 RNA:RNA hybrids
  cat( paste0( "Calculating guide::target hybridization energies started on " , date() ,"\n"))  
  
  # kind of slow
  # GuideSeq = DNAStringSet(Guide.df$GuideSeq)
  # Guide.df$hybMFE_3.12 = sapply( subseq( x = GuideSeq , start = 3 , width = 12) , GetRNAhybMFE )
  # Guide.df$hybMFE_15.9 = sapply( subseq( x = GuideSeq , start = 15 , width = 9) , GetRNAhybMFE )
  # this is ~ 2-fold faster
  Guide.df$hybMFE_3.12 = GetRNAhybMFE_bulk( dat = Guide.df , POS = 3 , WIDTH = 12  )
  Guide.df$hybMFE_15.9 = GetRNAhybMFE_bulk( dat = Guide.df , POS = 15 , WIDTH = 9  )
  
  
  
  
  
  
  
  # Get nucleotide densities of target sequences
  cat( paste0( "Get nucleotide density of target sequences started on " , date() ,"\n"))
  
  cors = read.delim(paste0(dir,'data/LocalNTdensityCorrelations.txt'), sep = '\t', header = T, stringsAsFactors = F)
  comb = cors[grep("combined",cors$Screen),]
  max = comb[grep("max",comb$COR),]
  min = comb[grep("min",comb$COR),]
  
  # This caluculates only the specific vector of the nt-density matrix that I am interested in
  NTdensitities = GetNTdensitities(fa = FA)
  Guide.df = GetNTpointdensities( DAT = Guide.df , Vec.list = NTdensitities )
  
  
  
  
  
  
  
  
  # Get nucleotide proportions in guide
  cat( paste0( "Get nucleotide proportions in guide RNA started on " , date() ,"\n"))
  
  LetterProbs = GetLetterProbs(x = Guide.df,  S=1, E=guideLength)
  Guide.df = cbind.data.frame( Guide.df, LetterProbs)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Generate Model based on original on-target screens. 
  cat( paste0( "Generate model started on " , date() ,"\n"))
  out = PrePareModelInput(x=fullset, FIELDS=fields)
  Model = out[[1]]
  # Get feature Means and SD for scaling
  ModelInputMeans = out[[2]]
  ModelInputSDs = out[[3]]
  
  # Standardize original Screen data to bin predicted guide scores according to screen efficacy quartiles
  Min=quantile(fullset$normCS, na.rm = T , probs = 0.05)
  Max=quantile(fullset$normCS, na.rm = T , probs = 0.95)
  fullset$standardizedGuideScores = sapply( fullset$normCS , FUN = standardizeScores , MIN = Min, MAX = Max)
  
  # Predict guide scores
  cat( paste0( "Guide scoring started on " , date() ,"\n"))
  GuideScores = PredictGuideScores(x = Guide.df, FIELDS=fields[-1], Fit=Model, Mean = ModelInputMeans, SD = ModelInputSDs)
  
  # Add guide Scores
  Guide.df = AddScores(df = Guide.df , SCORES = GuideScores)
  
  
  # Standardize GuideScores relative to ModelInput
  Min=quantile(fullset$normCS, na.rm = T , probs = 0.05)
  Max=quantile(fullset$normCS, na.rm = T , probs = 0.95)
  Guide.df$standardizedGuideScores = sapply( Guide.df$GuideScores , FUN = standardizeScores , MIN = Min, MAX = Max)
  Guide.df$standardizedGuideScores[Guide.df$standardizedGuideScores > 1] = 1
  Guide.df$standardizedGuideScores[Guide.df$standardizedGuideScores < 0] = 0
  
  # add quartile info
  Guide.df$quartiles = AssignQuartiles(Guide.df$standardizedGuideScores, q = quantile(fullset$standardizedGuideScores, na.rm=T))
  
  # remove feature columns
  Guide.df = Guide.df[, -which( colnames(Guide.df) %in% c( fields ,c('NTdens_min_GC','pT','pG|pC','pA|pT','pTT')) )]
  Guide.df = cbind.data.frame(rownames(Guide.df) , Guide.df)
  rownames(Guide.df)=NULL
  colnames(Guide.df)[1]='GuideName'
  
  
  
  # write guides to file system
  Seqs = DNAStringSet(Guide.df$GuideSeq)
  names(Seqs) =  paste0(rownames(Guide.df),"_",as.character(round(Guide.df$standardizedGuideScores,4)),"_",as.character(Guide.df$Rank),"_Q",Guide.df$quartiles)
  writeXStringSet(Seqs, filepath = paste0(prediction_dir, "/", name,"_","CasRxguides.fa"), append = F, format = "fasta")
  write.table(Guide.df, file = paste0(prediction_dir, "/", name,"_","CasRxguides.csv"), quote = F, sep=",", col.names = T, row.names = F)
  
  cat( paste0( "done " , date() ,"\n"))
  
  
  ######### optional plotting ############################################ 
  
  
  if (PLOT == TRUE){
    cds_coord = GetCDSregionBounderies(names(FA))
    pdf(paste0(prediction_dir, "/", name,"_CasRxguides.pdf"), width = 10, height = 3, useDingbats = F )
    suppressWarnings( gg_plot(df=Guide.df,LEN=width(FA), NAME = name , CDS = cds_coord, f=0.2) ) # There may be a warning if the smoothing exceeds the 0 to 1 y range
    dev.off()
  }
  
  #setwd("../")
  
  
  
  
  
  
  
}