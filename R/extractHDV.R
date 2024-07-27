#' @title Function for obtaining the HDV (High Dimensional Vector) matrix
#' @name extractHDV
#' @description Function for obtaining the HDV matrix without projecting it. 
#' The HDV corresponds to the high-dimensional matrix that represents the biological sequences
#' in a vectorial and structured way.
#' 
#' @param input There are two input formats available:
#'              (a) `BStringSet' (variants: `AAStringSet', `RNAStringSet', `DNAStringSet'). Biological sequence format loaded in memory;
#'              (b) `character'. String containing a path to a folder with FASTA files.
#' @param mask		readging mask. Default for amino acids is `c(2,1,2)` and for nucleotides c(5,5,5)#' 
#' @param seqtype   type of data: AA for amino acid, NT for nucleotide. The default is `AA`
#' @param bin       binary mode (TRUE), or counting mode (FALSE) for HDV construction. Default is TRUE
#' @param extension         extension of files desired to concatenate (Optional).   Available only for input type path to folder with FASTA files.
#' @param concatenate   defines whether to treat each sequence individually or to concatenate them into a single sequence. 
#'                  Available only for inputs in biological sequence format. The default is FALSE.
#' @param verbose   verbose mode. The default is TRUE
#' @param ...       other arguments of the function itself
#' 
#' @return 
#' `extractHDV' returns a `list` containing:
#' \itemize{
#'   \item HDV: a `matrix' containing the High Dimensional Vectors of the given FASTAS
#'   \item headers: a `character' containing the list of samples 
#'   \item mask: a `integer' containing the mask used
#'   \item SequenceType: a `character' containing the type of the sequence (amino acid: AA, ou nucleotide: NT)
#'   \item extension: a `character' containing the list of extensions considered
#'   \item bin: a `character' containing if binary or counting
#'   \item saturation: a `vector' containing the filled (non-zero) percentage of the HDV for each sample
#'   \item timeElapsed: a `double' containing the elapsed time in seconds
#' } 
#' 
#' @examples
#' 
#' # get the path to the folder containing the FASTA files
#' path = paste (system.file("examples/aaMitochondrial/",package = "rSWeeP"),'/', sep = '')
#' 
#' # define the parameters
#' mask = c(2,1,2)
#' 
#' # get the vectors that represent the sequences in high dimension (without projection)
#' HDV = extractHDV(input=path,mask=mask,seqtype='AA',bin=FALSE,extension=c('.faa','.fas','.fasta'))
#' 
#' 
#' @rdname extractHDV
#' @import foreach doParallel Biostrings methods
#' @export
setGeneric(
  "extractHDV",
  function(input, mask=NULL,seqtype='AA',...) standardGeneric("extractHDV"),
    signature = "input"
)
# setMethod("extractHDV", "character", function(path, extension,mask,seqtype,bin,verbose) {

.extractHDVfromFolder <-function(input, mask=NULL,seqtype='AA',bin=FALSE,extension='',verbose=TRUE){

    start_time = proc.time()

    mask = Defmask(mask,seqtype)

    ## Create the lists of files and names
    fastalist = list.files(input,pattern=extension,full.names=TRUE,recursive=TRUE)
    namesfasta = list.files(input,pattern=extension,full.names=FALSE,recursive=TRUE)
    # remove extension from filename
    namesfasta = tools::file_path_sans_ext(namesfasta)   
    
    
    ## input checks ---------------------- vv
    SW.checks('extension',extension)
    SW.checks('mask',mask)
    SW.checks('seqtype',seqtype)
    SW.checks('bin',bin)
    SW.checks('fastalist',fastalist)
    ## input checks ---------------------- ^^

	## some parameters ---------------------------------
	# spacer between fastas of same individual
    spacer = paste(rep("*",length(mask)),collapse='')
    # define maximum length of HDV
    if (seqtype=='AA'){
		lenmax = 20^sum(mask)
	} else if(seqtype=='NT'){
		lenmax = 4^sum(mask)
	}

	mask = as.integer(mask)

    
   	# number of sequences
    N = length(fastalist)
    N_ = as.character(N)


    # empty output  matrix
    output = createSWobjHDV(N,lenmax,namesfasta,mask,bin,seqtype,extension,FALSE)
    

    # for each fasta file, sweep
    for (k in 1:N){
		if(verbose){
            cat(paste('starting sequence ',as.character(k),'of', N_))
        }
        actualfasta = concatenaEach(fastalist[k],spacer)
       	seq = seq2num(actualfasta,seqtype)

        output$HDV[k,] = readmaskVEC(seq, mask,seqtype,bin,'none',lenmax)

		# atual sequence
		if(verbose){
            cat(paste(' - complete\n')) 
        }
    }
  
    

    if(N>1){    output$info$saturation = rowSums(output$HDV!=0)/lenmax    } 
    else {      output$info$saturation = sum(output$HDV!=0)/lenmax    }

    end_time = proc.time()
    output$info$timeElapsed = (end_time - start_time)[[3]] 
  	return(output)

}



 .extractHDVinRAM<-function(input, mask=NULL,seqtype=NULL,bin=FALSE,concatenate=FALSE, verbose=TRUE){
    start_time = proc.time()


    if (inherits(input,"AAStringSet"))                                      { seqtype='AA' }
    else if (inherits(input,"RNAStringSet") | inherits(input,"DNAStringSet")) { seqtype='NT' }
    else if (is.null(seqtype)){ 
        stop("Please provide the type of data. Use seqtype='AA' for amino acids and seqtype='NT' for nucleotides. ")
    }
 
    mask = Defmask(mask,seqtype)


    ## input checks ---------------------- vv
    SW.checks('mask',mask)
    SW.checks('seqtype',seqtype)
    SW.checks('bin',bin)
    SW.checks('concatenate',concatenate)
    ## input checks ---------------------- ^^


    if(concatenate){ # concatenate with spacer
        spacer = paste(rep("*",length(mask)),collapse='')
        input = stringi::stri_join_list(list(paste(input)),sep = spacer,collapse = NULL) # require stringi
    }

     # PARAMETERS --------------
    N = length(input) # number of sequences
    N_ = as.character(N)
    lenmax <- if (seqtype == 'AA') 20^sum(mask) else if (seqtype == 'NT') 4^sum(mask)
    mask = as.integer(mask)
    

    # empty output  matrix
    output = createSWobjHDV(N,lenmax,names(input),mask,bin,seqtype,'',concatenate)
    
    for (k in 1:N){
        if(verbose){ cat(paste('starting sequence ',as.character(k),'of', N_)) }

        seq = seq2num(input[k],seqtype)

        output$HDV[k,] = readmaskVEC(seq, mask,seqtype,bin,'none',lenmax)

        # atual sequence
        if(verbose){ cat(paste(' - complete\n')) }

    } # end for k in N

  
    if(N>1){    output$saturation = rowSums(output$HDV!=0)/lenmax    } 
    else {      output$saturation = sum(output$HDV!=0)/lenmax    }

    end_time = proc.time()
    output$timeElapsed = (end_time - start_time)[[3]] 

    return(output)

}



#' @rdname extractHDV
setMethod("extractHDV", "AAStringSet",  .extractHDVinRAM) 
#' @rdname extractHDV
setMethod("extractHDV", "DNAStringSet", .extractHDVinRAM) 
#' @rdname extractHDV
setMethod("extractHDV", "RNAStringSet", .extractHDVinRAM) 
#' @rdname extractHDV
setMethod("extractHDV", "BStringSet",   .extractHDVinRAM) 
#' @rdname extractHDV
setMethod("extractHDV", "BString",      .extractHDVinRAM) 
#' @rdname extractHDV
setMethod("extractHDV", "character",    .extractHDVfromFolder) 
