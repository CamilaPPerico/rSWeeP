# internal function, not to export


convertMask = function(x){
	#input  format: x = c(2,1,2)
	#output format: y = c(1,1,0,1,1)
	N=length(x)
	y=1:N%%2
	return(rep(y,x))
}



# -----------------------------------------------------------------------------
createSWobj <- function(output,psz,namesfasta,mask,bin,seqtype,extension,concatenate,norm,timeElapsed){

    output$info=NULL
    output$info$ProjectionSize =  psz
    output$info$headers = namesfasta
    output$info$mask = mask
    output$info$SequenceType = seqtype
    output$info$extension = extension
    output$info$bin = ifelse(bin, 'binary (TRUE)', 'counting (FALSE)')
    output$info$concatenate = concatenate
    output$info$version = packageVersion('rSWeeP')
    output$info$norm =  norm
    output$info$timeElapsed = timeElapsed

    return(output)
}




# -----------------------------------------------------------------------------
createSWobjHDV <- function(N,lenmax,namesfasta,mask,bin,seqtype,extension,concatenate){
    output 		   = NULL # create an output object
    output$info=NULL
    output$HDV = matrix(0,nrow=N,ncol=lenmax)
    output$info$headers = namesfasta
    output$info$mask = mask
    output$info$SequenceType = seqtype
    output$info$extension = extension
    output$info$concatenate = concatenate
    output$info$bin = ifelse(bin, 'binary (TRUE)', 'counting (FALSE)')
    output$info$version = packageVersion('rSWeeP')
    return(output)
}





# -----------------------------------------------------------------------------
concatenaEach <- function(namefile,spacer) {
	# read the multifasta file
	fastaFile <- Biostrings::readBStringSet(namefile) 
	seq = stringi::stri_join_list(list(fastaFile),sep = spacer,collapse = NULL) 
	return(seq)
}  





# -----------------------------------------------------------------------------
concatenaAll<- function(namefile,spacer,N) {
	ALLseqs = list()

	for (k in 1:N){
        ALLseqs[[k]] = concatenaEach(namefile[k],spacer)
    }
	
	return(Biostrings::BStringSet(unlist(ALLseqs)))
}  





# -----------------------------------------------------------------------------
seq2num <- function(fastaFile,seqtype){
  seq = toupper(unlist(strsplit(as.character(fastaFile[1]),"")))
  vecseq = unlist(strsplit(seq, split = ""))

  if (seqtype=='AA'){

    aa2int = data.frame(
    aa= c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*','-','U','O','J','?'),  
    val=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,0))

    # aa rares:
    # U = Selenocysteine (rare)
    # O = Pyrrolysine (rare)
    # aa dubious
    # J = Leucine (L) or Isoleucine (I)
    # B = Aspartic acid (D) or Asparagine (N)

    for (i in 1:29) {
      vecseq[vecseq==aa2int$aa[i]] = aa2int$val[i]
    }
    
  } else if(seqtype=='NT'){

    nt2int = data.frame(
    nt= c('A','C','G','T','U','R','Y','K','M','S','W','B','D','H','V','N','-','*','E','F','I','J','L','O','P','Q','Z'),
    val=c(1,2,3,4,4,5,6,7,8,9,10,11,12,13,14,15,16,0,0,0,0,0,0,0,0,0,0))

    for (i in 1:27) {
      vecseq[vecseq==nt2int$nt[i]] = nt2int$val[i]
    }
    # 27
  }

  seq = as.integer(vecseq)
  return(seq)
}



# -----------------------------------------------------------------------------
returnPOT <- function(seq,mask,seqtype){

	# length of mask
    lmer = length(mask) 
    m = sum(mask)

    if (seqtype=='AA'){
        KK = 20
    } else if(seqtype=='NT'){
        KK = 4
    }

    # "dna2list"
    nwindows = length(seq)-lmer+1
    pot = vector(mode="integer",length=nwindows)
    maskwin = matrix(0,nrow=nwindows,ncol=lmer)
    lmer2 = lmer-1

    
    for (i in 1:nwindows){
        maskwin[i,] = seq[i:(i+lmer2)]
    }

    if (lmer>m){
    	remove = which(mask==0)
	    maskwin = maskwin[,-remove]
    }
	maskwin = maskwin-1	
    

    REPETECO = matrix(rep( rep(KK,m)**c(0:(m-1)) ,nwindows),nrow=nwindows,byrow=T)

    pot = rowSums(maskwin*REPETECO )+1

    pot[which(rowSums( maskwin>(KK-1) | maskwin<0 ) >0 )] = -1

    return(pot[pot>0])

}






# -----------------------------------------------------------------------------
readmask <- function(seq, mask,seqtype,bin,norm,lenmax){
    
    pot = returnPOT(seq,mask,seqtype)

    aux = tabulate(pot,nbins=lenmax)

    outseq = list()
   	outseq$idx = which(aux!=0)

    if (bin) {

    	switch(norm,
    		'none'   = {	outseq$count = NULL		},
    		'log'    = {	outseq$count = NULL		},
    		'logNeg' = {
    					outseq$count = rep(-1,lenmax)
    					outseq$count[outseq$idx] = 1
    					outseq$idx = 1:lenmax 
	    			}
			) # end switch

	
	} else { # count ----------

		switch(norm,
    		'none'   = {	outseq$count = aux[outseq$idx]			},
    		'log'    = {	outseq$count = log(aux[outseq$idx],base=10)		},
    		'logNeg' = {
    					outseq$count = rep(-1,lenmax)
    					outseq$count[outseq$idx] = log(aux[outseq$idx],base=10)
    					outseq$idx = 1:lenmax 
    				}
			) # end switch

	}


    return(outseq)
}





# -----------------------------------------------------------------------------
readmaskVEC <- function(seq, mask,seqtype,bin,norm,lenmax){ 
   	
    pot = returnPOT(seq,mask,seqtype)
   
	outseq = tabulate(pot,nbins=lenmax) 

	if (bin) {
		outseq[outseq!=0] = 1

		if (norm=='logNeg') { 	outseq[outseq==0] = -1  	}
	
	} else { # count ----------

		switch(norm,
    		'log'    = {	outseq[outseq!=0] = log(outseq[outseq!=0],base=10)	},
    		'logNeg' = {
    						idx = which(outseq==0)
    					 	outseq[outseq!=0] = log(outseq[outseq!=0],base=10)
    					 	outseq[idx] = -1  	
    				}
			) # end switch
	}

	return(outseq)
	
}






# -----------------------------------------------------------------------------
readmasklist <- function(seqlist, mask,seqtype,bin,norm){
   	# length of mask
    lmer = length(mask) 
    maskBoo = as.logical(mask)
	m = sum(mask)

	outseqlist = list()
	# "dna2list"
	for (k in 1:length(seqlist)){
		outseqlist[[k]] = readmask(seqlist[[k]], mask,seqtype,bin,norm)
	}
	return(outseqlist)
}







# -----------------------------------------------------------------------------
seq2numlist <- function(namefile,seqtype){

	# read the multifasta file
	fastaFile <- Biostrings::readBStringSet(namefile) # require Biostrings
	# extract the sequences in a list
	sequences = paste(fastaFile)
	# concatenate with spacer
	outseq = list()
	for (k in 1:length(sequences)){
		outseq[[k]] = seq2num(sequences[[k]],seqtype)
	}


	return(outseq)
}
# -----------------------------------------------------------------------------





COREloop <- function(xnorm,psz,Mproj,pslist,nps,hdv_vec,nk,bin,norm){
    ntimes = ceiling(length(hdv_vec$idx)/nk);
    lim = 1
    for (i in 1:ntimes){
        lim = c(lim,nk*i)
    }
    lim[ntimes+1] = length(hdv_vec$idx) 

    

    x <- foreach(
      i = 1:ntimes, 
      .combine = 'cbind'
    ) %dopar% {
        # sqrt(i)
        mret = (matrix(rep(hdv_vec$idx[lim[i]:lim[i+1]],nps),ncol=nps))%%(matrix(rep(pslist,length(hdv_vec$idx[lim[i]:lim[i+1]])),ncol=nps,byrow=T))
        dt = (1+(mret))%*%Mproj
        if(!bin | norm=='logNeg'){
	        bs = (matrix(rep(hdv_vec$count[lim[i]:lim[i+1]],psz),ncol=psz))*((dt - floor(dt))-0.5)/0.5
        } else {
	        bs = ((dt - floor(dt))-0.5)/0.5
        }
        as.double(colSums(bs)) /xnorm
      }

    if(ntimes>1){ return(rowSums(x))  } 
    else { return(x) }
}




# -----------------------------------------------------------------------------
liteParam <- function(mask,input,seqtype,N,psz,lenmax){

    a = NULL
       
    set.seed(647474747) # fixed
    a$pslist = c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229)
    a$nps = length(a$pslist)
    a$Mproj =  matrix(stats::runif(a$nps * psz), ncol=psz)
    a$xnorm = sqrt(lenmax/3)

    return(a)
} # end function liteParam





# -----------------------------------------------------------------------------
NCoresDef <- function(ncores) {
	if (is.null(ncores) || ncores<=0){
		return(2) # Bioconductor requirement
	} else {
		return(ncores)
	}
}






# -----------------------------------------------------------------------------
MakeLOG <- function(input,bin,norm){

	if(norm=='log' & !bin){
        input[input!=0] = log(input[input!=0],base=10) # os log(0) -> -inf
	}
	if(norm=='logNeg'){
		if(bin) {
	        input[input==0]=-1
	    } else {
	        input = log(input,base=10) # os log(0) -> -inf
	        input[is.infinite(input)]=-1
	    }	
	}
    
    return(input)
}





# -----------------------------------------------------------------------------
Defmask <- function(mask,seqtype){
	
	if ( is.null(mask) & seqtype == 'AA' ){ # [212]
		return(c(1,1,0,1,1))
		# return(c(2,1,2))
	} else if (is.null(mask) & seqtype == 'NT' ){ # [555]
		# return(c(5,5,5))
		return(c(1,1,1,1,1,0,0,0,0,0,1,1,1,1,1))
	} else {
		if (max(mask)>1){
			mask = convertMask(mask)
			return(mask)
		}
		return(mask)
	}

}






# -----------------------------------------------------------------------------
# hdv_vec = matrix(0,nrow=N,ncol=lenmax)

HDVparallel <- function(N,ncores,input,seqtype,mask,bin,norm,lenmax){

	# START PARALLEL 
	ncores = NCoresDef(ncores)
	sw.cluster <- parallel::makeCluster(ncores, type = "FORK") 
	doParallel::registerDoParallel(cl = sw.cluster)
	foreach::getDoParWorkers()
	i=NULL
	x <- foreach(
		i = 1:N, 
		.combine = 'rbind',
		.export='seq2num'
		) %dopar% {
			seq=NULL # necessary
			seq = seq2num(input[i],seqtype)

			# retorna:
			readmaskVEC(seq, mask,seqtype,bin,norm,lenmax) 
	}


	# STOP PARALLEL 
	parallel::stopCluster(cl = sw.cluster)

	return(x)
}
