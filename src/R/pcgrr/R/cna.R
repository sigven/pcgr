
#' Function that gets the chromosome bands of copy number segments
#'
#' @param cna_gr genomic ranges object with copy number aberrations
#' @param cytoband_gr genomic ranges object with chromosomal cytobands
#' @param normalization_method metod for normalization of context counts (deconstructSigs)
#' @param signatures_limit max number of contributing signatures
#'
#'
get_cna_cytoband <- function(cna_gr, cytoband_gr){

  cyto_hits <- GenomicRanges::findOverlaps(cna_gr, cytoband_gr, type="any", select="all")
  ranges <- cytoband_gr[subjectHits(cyto_hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(cna_gr[queryHits(cyto_hits)]))
  cyto_df <- as.data.frame(mcols(ranges))
  cyto_df$segment_start <- start(ranges(cna_gr[queryHits(cyto_hits)]))
  cyto_df$segment_end <- end(ranges(cna_gr[queryHits(cyto_hits)]))
  cyto_df$segment_length <- width(ranges(cna_gr[queryHits(cyto_hits)]))

  cyto_stats <- as.data.frame(dplyr::group_by(cyto_df,segmentID,segment_length) %>% dplyr::summarise(cytoband = paste(name, collapse=", "), chromosome_arm = paste(unique(arm), collapse=","), focalCNAthresholds = paste(unique(focalCNAthreshold), collapse=",")))
  if(nrow(cyto_stats[stringr::str_detect(cyto_stats$focalCNAthresholds,","),]) > 0){
    cyto_stats[stringr::str_detect(cyto_stats$focalCNAthresholds,","),]$focalCNAthresholds <- NA
  }
  cyto_stats$focalCNAthresholds <- as.numeric(cyto_stats$focalCNAthresholds)
  cyto_stats$event_type <- 'broad'
  if(nrow(cyto_stats[!is.na(cyto_stats$focalCNAthresholds) & cyto_stats$segment_length < cyto_stats$focalCNAthresholds,]) > 0){
    cyto_stats[!is.na(cyto_stats$focalCNAthresholds) & cyto_stats$segment_length < cyto_stats$focalCNAthresholds,]$event_type <- 'focal'
  }

  cyto_stats <- pcgrr::df_string_replace(cyto_stats, c("cytoband"), pattern = ", (\\S+, ){0,}", replacement = " - ")
  cyto_stats <- tidyr::separate(cyto_stats,segmentID,sep=":",into = c('chrom','start','stop'),remove=F)
  cyto_stats$cytoband <- paste0(cyto_stats$chrom,":",cyto_stats$cytoband)
  cyto_stats <- dplyr::select(cyto_stats, segmentID, cytoband, event_type)

  return(cyto_stats)
}



#' Function that annotates CNV segment files
#'
#' @param cna_file CNV file name
#' @param logR_gain logR treshold for annotating copy number amplifications
#' @param logR_homdel logR threshold for annotating homozygous deletions
#' @param pcgr_data List of data frames with PCGR data annotations
#'
#' @return cna_data
#'

cna_segment_annotation <- function(cna_file, logR_gain, logR_homdel, pcgr_data, transcript_overlap_frac = 100){

  rlogging::message(paste0("Annotation of copy number segment file ",cna_file))
  cna_df_raw <- read.table(file=cna_file,header = T,stringsAsFactors = F,comment.char="", quote="")
  cna_df_raw <- dplyr::rename(cna_df_raw, chromosome = Chromosome, LogR = Segment_Mean, segment_start = Start, segment_end = End) %>%
    dplyr::distinct()
  cna_df_raw$chromosome <- stringr::str_replace(cna_df_raw$chromosome,"^chr","")

  ## VALIDATE INPUT CHROMOSOMES
  cna_df <- pcgrr::get_valid_chromosomes(cna_df_raw, chromosome_column = 'chromosome', bsg = BSgenome.Hsapiens.UCSC.hg19)
  cna_df <- cna_df %>% dplyr::filter(!is.na(LogR))
  cna_df$LogR <- as.numeric(cna_df$LogR)
  cna_df$LogR <- round(cna_df$LogR,digits=2)
  cna_df$segmentID <- paste0(cna_df$chromosome,":",cna_df$segment_start,":",cna_df$segment_end)

  ## MAKE GRANGES OBJECT OF INPUT
  cna_gr <- GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T, seqinfo = pcgr_data$seqinfo_hg19, seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end', ignore.strand = T, starts.in.df.are.0based = T)
  cytoband_df <- pcgrr::get_cna_cytoband(cna_gr, pcgr_data$cytoband_gr)
  cna_df <- dplyr::left_join(cna_df, cytoband_df,by="segmentID")

  cna_segments <- cna_df
  cna_segments$segment_link <- paste0("<a href='",paste0('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',paste0(cna_segments$chromosome,':',cna_segments$segment_start,'-',cna_segments$segment_end)),"' target=\"_blank\">",paste0(cna_segments$chromosome,':',cna_segments$segment_start,'-',cna_segments$segment_end),"</a>")
  cna_segments$segment_length_Mb <- round((as.numeric((cna_segments$segment_end - cna_segments$segment_start)/1000000)),digits = 4)
  cna_segments <- dplyr::rename(cna_segments, SEGMENT_LENGTH_MB = segment_length_Mb, SEGMENT = segment_link)
  cna_segments <- dplyr::select(cna_segments, SEGMENT, SEGMENT_LENGTH_MB, cytoband, LogR, event_type) %>% dplyr::distinct()

  cna_gr <- GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T, seqinfo = pcgr_data$seqinfo_hg19, seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end', ignore.strand = T, starts.in.df.are.0based = T)

  hits <- GenomicRanges::findOverlaps(cna_gr, pcgr_data$gencode_genes_gr, type="any", select="all")
  ranges <- pcgr_data$gencode_genes_gr[subjectHits(hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(cna_gr[queryHits(hits)]))

  df <- as.data.frame(mcols(ranges))
  df$segment_start <- start(ranges(cna_gr[queryHits(hits)]))
  df$segment_end <- end(ranges(cna_gr[queryHits(hits)]))
  df$segment_length_Mb <- round((as.numeric((df$segment_end - df$segment_start)/1000000)),digits = 4)

  df$transcript_start <- start(ranges)
  df$transcript_end <- end(ranges)
  df$chrom <- as.character(seqnames(ranges))
  df <- as.data.frame(df %>% dplyr::rowwise() %>% dplyr::mutate(transcript_overlap_percent = round(as.numeric((min(transcript_end,segment_end) - max(segment_start,transcript_start)) / (transcript_end - transcript_start)) * 100, digits = 2)))

  df$segment_link <- paste0("<a href='",paste0('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',paste0(df$chrom,':',df$segment_start,'-',df$segment_end)),"' target=\"_blank\">",paste0(df$chrom,':',df$segment_start,'-',df$segment_end),"</a>")
  df_print <- df
  df_print <- dplyr::select(df_print,chrom,segment_start,segment_end,segment_length_Mb,event_type,cytoband,LogR,ensembl_gene_id,symbol,ensembl_transcript_id,transcript_start,transcript_end,transcript_overlap_percent,name,biotype,disgenet_cui,tsgene,ts_oncogene,intogen_drivers,chembl_compound_id,gencode_gene_biotype,gencode_tag,gencode_release)


  chrOrder <- c(as.character(paste0('chr',c(1:22))),"chrX","chrY")
  df_print$chrom <- factor(df_print$chrom, levels=chrOrder)
  df_print <- df_print[order(df_print$chrom),]
  df_print$segment_start <- as.integer(df_print$segment_start)
  df_print$segment_end <- as.integer(df_print$segment_end)

  df_print_sorted <- NULL
  for(chrom in chrOrder){
    if(nrow(df_print[df_print$chrom == chrom,]) > 0){
      chrom_regions <- df_print[df_print$chrom == chrom,]
      chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(segment_start, segment_end)),]
      df_print_sorted <- rbind(df_print_sorted, chrom_regions_sorted)
    }
  }

  df <- dplyr::select(df, -ensembl_transcript_id) %>% dplyr::filter(gencode_gene_biotype == 'protein_coding') %>% dplyr::distinct()

  df <- dplyr::rename(df, CHEMBL_COMPOUND_ID = chembl_compound_id)
  df$VAR_ID <- rep(1:nrow(df))
  df <- pcgrr::annotate_variant_link(df, vardb = 'DGIDB', pcgr_data = pcgr_data)
  df <- dplyr::rename(df, ONCOGENE = ts_oncogene, TUMOR_SUPPRESSOR = tsgene, ENTREZ_ID = entrezgene, GENE = symbol, CHROMOSOME = chrom, GENENAME = name, TARGETED_DRUGS = DGIDBLINK, SEGMENT_LENGTH_MB = segment_length_Mb, SEGMENT = segment_link, TRANSCRIPT_OVERLAP = transcript_overlap_percent)
  df$ENTREZ_ID <- as.character(df$ENTREZ_ID)
  df <- dplyr::left_join(df,pcgr_data$kegg_gene_pathway_links, by=c("ENTREZ_ID" = "gene_id"))
  df <- dplyr::rename(df, KEGG_PATHWAY = kegg_pathway_urls)
  df <- pcgrr::annotate_variant_link(df, vardb = 'NCBI_GENE')
  df <- dplyr::rename(df, GENE_NAME = NCBI_GENE_LINK)

  df <- dplyr::select(df, CHROMOSOME, GENE, GENE_NAME, KEGG_PATHWAY, TUMOR_SUPPRESSOR, ONCOGENE, TARGETED_DRUGS,SEGMENT_LENGTH_MB, SEGMENT, gencode_gene_biotype,LogR, TRANSCRIPT_OVERLAP) %>% dplyr::distinct()
  df <- df %>% dplyr::distinct()
  targeted_drugs <- dplyr::select(df, GENE, TARGETED_DRUGS) %>% dplyr::filter(!is.na(TARGETED_DRUGS)) %>% dplyr::distinct()

  cna_segments_filtered <- data.frame()
  cna_segments_filtered <- dplyr::filter(cna_segments, LogR >= logR_gain | LogR <= logR_homdel)
  cna_segments_filtered <- cna_segments_filtered %>% dplyr::arrange(desc(LogR))
  rlogging::message(paste0("Detected ",nrow(cna_segments_filtered)," segments subject to amplification/deletion"))

  oncogene_amplified <- data.frame()
  oncogene_amplified <- dplyr::filter(df, !is.na(ONCOGENE) & TRANSCRIPT_OVERLAP == transcript_overlap_frac & LogR >= logR_gain & gencode_gene_biotype == 'protein_coding')
  if(nrow(oncogene_amplified) > 0){
    oncogene_amplified <- dplyr::select(oncogene_amplified, -c(TUMOR_SUPPRESSOR, ONCOGENE,TRANSCRIPT_OVERLAP,gencode_gene_biotype,TARGETED_DRUGS)) %>% dplyr::distinct()
    oncogene_amplified <- dplyr::left_join(oncogene_amplified, targeted_drugs, by="GENE")
    oncogene_amplified <- dplyr::arrange(oncogene_amplified, TARGETED_DRUGS, LogR)
    oncogene_amplified <- dplyr::select(oncogene_amplified, CHROMOSOME, GENE, GENE_NAME, KEGG_PATHWAY, TARGETED_DRUGS, dplyr::everything())
    oncogene_amplified$CNA_TYPE <- 'gain'
    rlogging::message(paste0("Detected proto-oncogene(s) subject to amplification (log(2) ratio >= ",logR_gain,"): ",paste0(oncogene_amplified$GENE,collapse=", ")))
  }
  tsgene_homozygous_deletion <- data.frame()
  tsgene_homozygous_deletion <- dplyr::filter(df, !is.na(TUMOR_SUPPRESSOR) & TRANSCRIPT_OVERLAP == transcript_overlap_frac  & LogR <= logR_homdel & gencode_gene_biotype == 'protein_coding')
  if(nrow(tsgene_homozygous_deletion) > 0){
    tsgene_homozygous_deletion <- dplyr::select(tsgene_homozygous_deletion, -c(TUMOR_SUPPRESSOR, ONCOGENE,TRANSCRIPT_OVERLAP,gencode_gene_biotype,TARGETED_DRUGS)) %>% dplyr::distinct()
    tsgene_homozygous_deletion <- dplyr::left_join(tsgene_homozygous_deletion, targeted_drugs, by="GENE")
    tsgene_homozygous_deletion <- dplyr::arrange(tsgene_homozygous_deletion, TARGETED_DRUGS, LogR)
    tsgene_homozygous_deletion <- dplyr::select(tsgene_homozygous_deletion, CHROMOSOME, GENE, GENE_NAME, KEGG_PATHWAY, TARGETED_DRUGS, dplyr::everything())

    tsgene_homozygous_deletion$CNA_TYPE <- 'loss'
    rlogging::message(paste0("Detected tumor suppressor gene(s) subject to homozygous deletions (log(2) ratio <= ",logR_homdel,"): ",paste0(tsgene_homozygous_deletion$GENE,collapse=", ")))
  }
  civic_cna_biomarkers <- dplyr::filter(pcgr_data$civic_biomarkers, alteration_type == 'CNA' & !is.na(eitem_consequence)) %>% dplyr::select(genesymbol,evidence_type,evidence_level,evidence_description,cancer_type,evidence_direction,pubmed_html_link,therapeutic_context,rating,clinical_significance,eitem_consequence)
  names(civic_cna_biomarkers) <- toupper(names(civic_cna_biomarkers))
  civic_cna_biomarkers <- dplyr::rename(civic_cna_biomarkers, GENE = GENESYMBOL, CNA_TYPE = EITEM_CONSEQUENCE, DESCRIPTION = EVIDENCE_DESCRIPTION, CITATION = PUBMED_HTML_LINK)

  cbmdb_cna_biomarkers <- dplyr::filter(pcgr_data$cbmdb_biomarkers, alteration_type == 'CNA') %>% dplyr::select(genesymbol,evidence_type,evidence_level,evidence_description,cancer_type,evidence_direction,pubmed_html_link,therapeutic_context,rating,clinical_significance,variant_name)
  names(cbmdb_cna_biomarkers) <- toupper(names(cbmdb_cna_biomarkers))
  cbmdb_cna_biomarkers <- dplyr::rename(cbmdb_cna_biomarkers, GENE = GENESYMBOL, CNA_TYPE = VARIANT_NAME, DESCRIPTION = EVIDENCE_DESCRIPTION, CITATION = PUBMED_HTML_LINK)

  cna_biomarkers <- NULL
  cna_biomarker_segments <- data.frame()
  cna_evidence_items_diagnostic <- data.frame()
  cna_evidence_items_predisposing <- data.frame()
  cna_evidence_items_prognostic <- data.frame()
  cna_evidence_items_predictive <- data.frame()
  if(!is.null(tsgene_homozygous_deletion)){
    if(nrow(tsgene_homozygous_deletion) > 0){
      biomarker_hits <- NULL
      biomarker_hits_civic <- as.data.frame(dplyr::inner_join(tsgene_homozygous_deletion, civic_cna_biomarkers, by=c("GENE","CNA_TYPE")))
      if(nrow(biomarker_hits_civic) == 0){
        biomarker_hits_cbmdb <- as.data.frame(dplyr::inner_join(tsgene_homozygous_deletion, cbmdb_cna_biomarkers, by=c("GENE","CNA_TYPE")))
        biomarker_hits <- biomarker_hits_cbmdb
      }
      else{
        biomarker_hits <- biomarker_hits_civic
      }
      if(nrow(biomarker_hits) > 0){
        cna_biomarkers <- rbind(cna_biomarkers,biomarker_hits)
      }
    }
  }
  if(!is.null(oncogene_amplified)){
    if(nrow(oncogene_amplified) > 0){
      biomarker_hits <- NULL
      biomarker_hits_civic <- as.data.frame(dplyr::inner_join(oncogene_amplified, civic_cna_biomarkers, by=c("GENE","CNA_TYPE")))
      if(nrow(biomarker_hits_civic) == 0){
        biomarker_hits_cbmdb <- as.data.frame(dplyr::inner_join(oncogene_amplified, cbmdb_cna_biomarkers, by=c("GENE","CNA_TYPE")))
        biomarker_hits <- biomarker_hits_cbmdb
      }
      else{
        biomarker_hits <- biomarker_hits_civic
      }
      if(nrow(biomarker_hits) > 0){
        cna_biomarkers <- rbind(cna_biomarkers,biomarker_hits)
      }
    }
  }

  if(!is.null(cna_biomarkers)){
    if(nrow(cna_biomarkers) > 0){
      cna_biomarkers <- cna_biomarkers[c("GENE","CANCER_TYPE","CNA_TYPE","EVIDENCE_LEVEL","CLINICAL_SIGNIFICANCE","EVIDENCE_TYPE","DESCRIPTION","EVIDENCE_DIRECTION","THERAPEUTIC_CONTEXT","CITATION","RATING","GENE_NAME","KEGG_PATHWAY","TARGETED_DRUGS","SEGMENT_LENGTH_MB", "SEGMENT","LogR")]
      cna_biomarkers <- cna_biomarkers %>% dplyr::arrange(EVIDENCE_LEVEL,RATING)

      if(nrow(dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Predictive')) > 0){
        cna_evidence_items_predictive <- dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Predictive')
      }
      if(nrow(dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Prognostic')) > 0){
        cna_evidence_items_prognostic <- dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Prognostic')
      }
      if(nrow(dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Diagnostic')) > 0){
        cna_evidence_items_diagnostic <- dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Diagnostic')
      }
      if(nrow(dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Predisposing')) > 0){
        cna_evidence_items_predisposing <- dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == 'Predisposing')
      }
      cna_biomarker_segments <- dplyr::select(cna_biomarkers, SEGMENT, LogR) %>% dplyr::distinct()
    }
  }

  cna_data <- list(ranked_segments = cna_segments_filtered, oncogene_amplified = oncogene_amplified, tsgene_homozygous_deletion = tsgene_homozygous_deletion,cna_df_for_print = df_print_sorted, cna_biomarkers = cna_biomarkers, cna_biomarker_segments = cna_biomarker_segments, cna_evidence_items_diagnostic = cna_evidence_items_diagnostic, cna_evidence_items_predictive = cna_evidence_items_predictive, cna_evidence_items_prognostic = cna_evidence_items_prognostic, cna_evidence_items_predisposing = cna_evidence_items_predisposing)
  return(cna_data)
}

