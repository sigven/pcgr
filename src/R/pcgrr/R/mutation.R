
#' Function that annotates underlying mechanism (repeat-mediated, micro-homology etc) of indels
#'
#' @param vcf_data_df Data frame with somatic variants
#'
#' @return vcf_data_df
#'

get_deletion_mechanism <- function(vcf_data_df, genome_assembly = 'hg19'){

  rlogging::message("Determination of deletion mechanisms (repeat, microhomology etc) by investigation of junction sequence")
  genome_seq <- BSgenome.Hsapiens.UCSC.hg38
  if(genome_assembly == 'hg19'){
    genome_seq <- BSgenome.Hsapiens.UCSC.hg19
  }
  vcf_data_df_valid <- pcgrr::get_valid_chromosomes(vcf_data_df, chromosome_column = 'CHROM', bsg = genome_seq)
  deletions <- vcf_data_df_valid[vcf_data_df_valid$VARIANT_CLASS == 'deletion',]
  if(nrow(deletions) > 0){
    deletions <- dplyr::select(deletions, CHROM, POS, end, REF, ALT, VARIANT_CLASS, GENOMIC_CHANGE)
    deletions$DELETED_BASES <- substring(deletions$REF,2)
    deletions <- dplyr::filter(deletions, nchar(DELETED_BASES) > 1)
    if(nrow(deletions) > 0){
      seqinfo <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(genome_seq)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(genome_seq)), genome = genome_assembly)
      vcf_df_gr <- GenomicRanges::makeGRangesFromDataFrame(deletions, keep.extra.columns = T, seqinfo = seqinfo, seqnames.field = 'CHROM',start.field = 'POS', end.field = 'end', ignore.strand = T, starts.in.df.are.0based = F)

      downstream_gr <- GenomicRanges::flank(vcf_df_gr, width = nchar(mcols(vcf_df_gr)$DELETED_BASES), start=F)
      downstream_seq <- Biostrings::getSeq(genome_seq, downstream_gr)
      upstream_gr <- GenomicRanges::flank(vcf_df_gr, width = nchar(mcols(vcf_df_gr)$DELETED_BASES), start=T)
      upstream_seq <- Biostrings::getSeq(genome_seq, upstream_gr)
      df_up <- data.frame('FLANK_UPSTREAM'=toupper(unlist(strsplit(toString(upstream_seq),", "))),stringsAsFactors=F)
      df_down <- data.frame('FLANK_DOWNSTREAM'=toupper(unlist(strsplit(toString(downstream_seq),", "))),stringsAsFactors=F)
      deletions_flank_df <- cbind(deletions, df_up, df_down)
      deletions_flank_df$downstream_startswith_del <- FALSE
      if(nrow(deletions_flank_df[substring(deletions_flank_df$DELETED_BASES,1,2) == substring(deletions_flank_df$FLANK_DOWNSTREAM,1,2),]) > 0){
        deletions_flank_df[substring(deletions_flank_df$DELETED_BASES,1,2) == substring(deletions_flank_df$FLANK_DOWNSTREAM,1,2),]$downstream_startswith_del <- TRUE
      }
      deletions_flank_df$DELETION_MECHANISM <- NA
      if(nrow(deletions_flank_df[deletions_flank_df$DELETED_BASES == deletions_flank_df$FLANK_DOWNSTREAM,]) > 0){
        deletions_flank_df[deletions_flank_df$DELETED_BASES == deletions_flank_df$FLANK_DOWNSTREAM,]$DELETION_MECHANISM <- 'Repeat-mediated'
      }
      if(nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & deletions_flank_df$downstream_startswith_del == TRUE,]) > 0){
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & deletions_flank_df$downstream_startswith_del == TRUE,]$DELETION_MECHANISM <- 'Microhomology'
      }
      if(nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & stringr::str_detect(deletions_flank_df$DELETED_BASES,'^(TTT|AAA)') & stringr::str_detect(deletions_flank_df$FLANK_UPSTREAM,'(TTT|AAA)$'),]) > 0){
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & stringr::str_detect(deletions_flank_df$DELETED_BASES,'^(TTT|AAA)') & stringr::str_detect(deletions_flank_df$FLANK_UPSTREAM,'(TTT|AAA)$'),]$DELETION_MECHANISM <- 'Mononucleotide-repeat-mediated'
      }
      if(nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM),]) > 0){
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM),]$DELETION_MECHANISM <- 'Other'
      }

      deletions_flank_df <- dplyr::select(deletions_flank_df, GENOMIC_CHANGE, DELETION_MECHANISM)
      vcf_data_df <- dplyr::left_join(vcf_data_df, deletions_flank_df, by=c("GENOMIC_CHANGE"))

    }
  }

  return(vcf_data_df)

}


#' Function that assigns one of six mutation types to a list of mutations
#'
#' @param var_df data frame with variants
#'
#' @return var_df
#'
assign_mutation_type <- function(var_df){

  var_df$MUTATION_TYPE <- NA
  if('VARIANT_CLASS' %in% colnames(var_df) & 'REF' %in% colnames(var_df) & 'ALT' %in% colnames(var_df)){
    var_df <- var_df %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "G" & ALT == "A","C>T:G>A",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "C" & ALT == "T","C>T:G>A",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "G" & ALT == "C","C>G:G>C",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "C" & ALT == "G","C>G:G>C",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "C" & ALT == "A","C>A:G>T",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "G" & ALT == "T","C>A:G>T",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "A" & ALT == "G","A>G:T>C",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "T" & ALT == "C","A>G:T>C",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "A" & ALT == "T","A>T:T>A",MUTATION_TYPE)) %>%
      dplyr::mutate(MUTATION_TYPE = dplyr::if_else(VARIANT_CLASS == "SNV" & REF == "T" & ALT == "A","A>T:T>A",MUTATION_TYPE))
  }
  return(var_df)
}


get_proper_maf_alleles <- function(maf_df, genome_seq, seqinfo){

  maf_df_valid <- pcgrr::get_valid_chromosomes(maf_df, chromosome_column = 'Chromosome', bsg = genome_seq)
  if("end" %in% colnames(maf_df_valid)){
    maf_df_valid <- dplyr::select(maf_df_valid,-end)
  }

  maf_SNV <- dplyr::filter(maf_df_valid, Variant_Type == 'SNP')
  maf_SNV$REF <- maf_SNV$Reference_Allele
  maf_SNV$ALT <- maf_SNV$Tumor_Seq_Allele2
  maf_SNV$POS <- maf_SNV$Start_Position
  maf_ALL <- maf_SNV

  maf_INS <- dplyr::filter(maf_df_valid, Variant_Type == 'INS')
  maf_DEL <- dplyr::filter(maf_df_valid, Variant_Type == 'DEL')

  if(nrow(maf_DEL) > 0){
    ## get appropriate alleles (VCF-like) of reference and alternate (DELETIONS)
    maf_del_gr <- GenomicRanges::makeGRangesFromDataFrame(maf_DEL, keep.extra.columns = T, seqinfo = seqinfo, seqnames.field = 'Chromosome',start.field = 'Start_Position', end.field = 'End_Position', ignore.strand = T, starts.in.df.are.0based = F)

    maf_del_flank_gr <- GenomicRanges::flank(maf_del_gr, width = 1, start=T)
    maf_del_flank_seq <- Biostrings::getSeq(genome_seq, maf_del_flank_gr)
    maf_del_seq <- Biostrings::getSeq(genome_seq, maf_del_gr)
    vcf_alleles_alt <- data.frame('ALT'=toupper(unlist(strsplit(toString(maf_del_flank_seq),", "))),stringsAsFactors=F)
    vcf_alleles_ref <- data.frame('REF'=toupper(unlist(strsplit(toString(maf_del_seq),", "))),stringsAsFactors=F)
    vcf_alleles <- cbind(vcf_alleles_ref, vcf_alleles_alt)
    vcf_alleles$REF <- paste0(vcf_alleles$ALT,vcf_alleles$REF)
    maf_DEL <- cbind(maf_DEL,vcf_alleles)
    maf_DEL$POS <- maf_DEL$Start_Position - 1

    maf_ALL <- rbind(maf_ALL, maf_DEL)
  }

  if(nrow(maf_INS) > 0){
    ## get appropriate alleles (VCF-like) of reference and alternate (INSERTIONS)
    maf_ins_gr <- GenomicRanges::makeGRangesFromDataFrame(maf_INS, keep.extra.columns = T, seqinfo = seqinfo, seqnames.field = 'Chromosome',start.field = 'Start_Position', end.field = 'Start_Position', ignore.strand = T, starts.in.df.are.0based = F)
    maf_ins_seq <- Biostrings::getSeq(genome_seq, maf_ins_gr)
    vcf_alleles_alt <- data.frame('REF'=toupper(unlist(strsplit(toString(maf_ins_seq),", "))),stringsAsFactors=F)
    maf_INS <- cbind(maf_INS,vcf_alleles_alt)
    maf_INS$ALT <- paste0(maf_INS$REF,maf_INS$Tumor_Seq_Allele2)
    maf_INS$POS <- maf_INS$Start_Position

    maf_ALL <- rbind(maf_ALL, maf_INS)
  }


  maf_ALL$CHROM <- stringr::str_replace(maf_ALL$Chromosome,'chr','')
  maf_ALL$GENOMIC_CHANGE <- paste(paste(paste(paste0("g.chr",maf_ALL$CHROM),maf_ALL$POS,sep=":"),maf_ALL$REF,sep=":"),maf_ALL$ALT,sep=">")
  return(maf_ALL)

}

