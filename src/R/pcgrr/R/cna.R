
#' Function that removes copy number segments that go beyond chrosomal lengths for the given assembly
#'
#' @param cna_segments_df genomic segments with CNA
#' @param genome_assembly genome assembly (grch37/grch38)
#' @param bsg BSgenome object
#'
get_valid_chromosome_segments <- function(cna_segments_df, genome_assembly = 'grch37', bsg = BSgenome.Hsapiens.UCSC.hg19){
  chromosome_lengths <- data.frame('chromosome' = head(names(seqlengths(bsg)),24), 'chrom_length' = head(seqlengths(bsg),24), stringsAsFactors = F, row.names = NULL)
  cna_segments_df <- dplyr::left_join(cna_segments_df, chromosome_lengths,by=c("chromosome"))
  cna_segments_df <- as.data.frame(cna_segments_df %>% dplyr::rowwise() %>% dplyr::mutate(segment_error = segment_end > chrom_length))
  if(nrow(dplyr::filter(cna_segments_df, segment_error == T)) > 0){
    n_removed <- nrow(dplyr::filter(cna_segments_df, segment_error == T))
    rlogging::warning('Skipping ',n_removed,' copy number segments that span beyond chromosomal lengths for ', genome_assembly,' (make sure chromosomal segments are consistent with assembly)')
  }
  cna_segments_df <- dplyr::filter(cna_segments_df, segment_error == F) %>% dplyr::select(-c(segment_error,chrom_length))
  return(cna_segments_df)
}

#' Function that gets the chromosome bands of copy number segments
#'
#' @param cna_gr genomic ranges object with copy number aberrations
#' @param cytoband_gr genomic ranges object with chromosomal cytobands
#' @param normalization_method metod for normalization of context counts (deconstructSigs)
#' @param signatures_limit max number of contributing signatures
#'
get_cna_cytoband <- function(cna_gr, cytoband_gr){

  cyto_hits <- GenomicRanges::findOverlaps(cna_gr, cytoband_gr, type="any", select="all")
  ranges <- cytoband_gr[subjectHits(cyto_hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(cna_gr[queryHits(cyto_hits)]))
  cyto_df <- as.data.frame(mcols(ranges))
  cyto_df$segment_start <- start(ranges(cna_gr[queryHits(cyto_hits)]))
  cyto_df$segment_end <- end(ranges(cna_gr[queryHits(cyto_hits)]))
  cyto_df$segment_length <- width(ranges(cna_gr[queryHits(cyto_hits)]))

  cyto_stats <- as.data.frame(
      dplyr::group_by(cyto_df,segmentID,segment_length) %>%
      dplyr::summarise(cytoband = paste(name, collapse=", "), chromosome_arm = paste(unique(arm), collapse=","), focalCNAthresholds = paste(unique(focalCNAthreshold), collapse=","))
  )

  cyto_stats <- cyto_stats %>%
    dplyr::mutate(focalCNAthresholds = dplyr::if_else(!is.na(focalCNAthresholds) & stringr::str_detect(focalCNAthresholds,","),as.character(NA),as.character(focalCNAthresholds))) %>%
    dplyr::mutate(focalCNAthresholds = as.numeric(focalCNAthresholds))

  cyto_stats <- as.data.frame(cyto_stats %>% dplyr::rowwise() %>% dplyr::mutate(broad_cnv_event = segment_length > focalCNAthresholds))

  cyto_stats$event_type <- 'broad'
  cyto_stats <- cyto_stats %>%
    dplyr::mutate(event_type = dplyr::if_else(!is.na(broad_cnv_event) & broad_cnv_event == F,'focal',event_type)) %>%
    dplyr::mutate(event_type = dplyr::if_else(is.na(broad_cnv_event),as.character(NA),event_type))

  cyto_stats <- pcgrr::df_string_replace(cyto_stats, c("cytoband"), pattern = ", (\\S+, ){0,}", replacement = " - ")
  cyto_stats <- tidyr::separate(cyto_stats,segmentID,sep=":",into = c('chrom','start','stop'),remove=F)
  cyto_stats$cytoband <- paste0(cyto_stats$chrom,":",cyto_stats$cytoband)
  cyto_stats <- dplyr::select(cyto_stats, segmentID, cytoband, event_type)

  return(cyto_stats)
}

#' Function that annotates CNV segment files
#'
#' @param cna_file CNV file name with chromosomal log(2)-ratio segments
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param transcript_overlap_pct required aberration overlap fraction (percent) for reported transcripts (default 100 percent)
#'
generate_report_data_cna <- function(cna_file, pcgr_data, sample_name, pcgr_config, transcript_overlap_pct = 100){

  pcg_report_cna <- pcgrr::init_pcg_report(pcgr_config, sample_name, class = 'cna')
  logR_homdel <- pcgr_config[['cna']][['logR_homdel']]
  logR_gain <- pcgr_config[['cna']][['logR_gain']]

  rlogging::message('------')
  rlogging::message(paste0("Generating report data for copy number segment file ",cna_file))
  cna_df_raw <- read.table(file=cna_file,header = T,stringsAsFactors = F,sep="\t",comment.char="", quote="") %>%
    dplyr::rename(chromosome = Chromosome, LogR = Segment_Mean, segment_start = Start, segment_end = End) %>%
    dplyr::distinct() %>%
    dplyr::mutate(chromosome = stringr::str_replace(chromosome,"^chr",""))

  ## VALIDATE INPUT CHROMOSOMES AND SEGMENTS
  cna_df <- pcgrr::get_valid_chromosomes(cna_df_raw, chromosome_column = 'chromosome', bsg = pcgr_data[['assembly']][['bsg']])
  cna_df <- cna_df %>%
    pcgrr::get_valid_chromosome_segments(genome_assembly = pcgr_data[['assembly']][['grch_name']], bsg = pcgr_data[['assembly']][['bsg']]) %>%
    dplyr::filter(!is.na(LogR)) %>%
    dplyr::mutate(LogR = round(as.numeric(LogR),digits=3)) %>%
    dplyr::mutate(segmentID = paste0(chromosome,":",segment_start,":",segment_end))

  ## MAKE GRANGES OBJECT OF INPUT
  cna_gr <- GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T, seqinfo = pcgr_data[['assembly']][['seqinfo']],
                                                    seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end',
                                                    ignore.strand = T, starts.in.df.are.0based = T)
  cytoband_df <- pcgrr::get_cna_cytoband(cna_gr, pcgr_data[['genomic_ranges']][['cytoband']])
  cna_df <- dplyr::left_join(cna_df, cytoband_df,by="segmentID")

  cna_segments <- cna_df
  cna_segments <- cna_segments %>%
    pcgrr::add_ucsc_segment_link(hgname = pcgr_data[['assembly']][['hg_name']], chrom = "chromosome", start = "segment_start", end = "segment_end") %>%
    dplyr::mutate(segment_length_Mb = round((as.numeric((segment_end - segment_start)/1000000)),digits = 4)) %>%
    dplyr::rename(SEGMENT_LENGTH_MB = segment_length_Mb, SEGMENT = segment_link) %>%
    dplyr::select(SEGMENT, SEGMENT_LENGTH_MB, cytoband, LogR, event_type) %>%
    dplyr::distinct()

  cna_gr <- GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T, seqinfo = pcgr_data[['assembly']][['seqinfo']],
                                                    seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end',
                                                    ignore.strand = T, starts.in.df.are.0based = T)

  hits <- GenomicRanges::findOverlaps(cna_gr, pcgr_data[['genomic_ranges']][['gencode_genes']], type="any", select="all")
  ranges <- pcgr_data[['genomic_ranges']][['gencode_genes']][subjectHits(hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(cna_gr[queryHits(hits)]))

  local_df <- as.data.frame(as.data.frame(mcols(ranges)) %>%
    dplyr::mutate(segment_start = as.integer(start(ranges(cna_gr[queryHits(hits)])))) %>%
    dplyr::mutate(segment_end = as.integer(end(ranges(cna_gr[queryHits(hits)])))) %>%
    dplyr::mutate(segment_length_Mb = round((as.numeric((segment_end - segment_start)/1000000)),digits = 4)) %>%
    dplyr::mutate(transcript_start = start(ranges)) %>%
    dplyr::mutate(transcript_end = end(ranges)) %>%
    dplyr::mutate(sample_id = sample_name) %>%
    dplyr::mutate(chrom = as.character(seqnames(ranges))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(transcript_overlap_percent =
                    round(as.numeric((min(transcript_end,segment_end) - max(segment_start,transcript_start)) / (transcript_end - transcript_start)) * 100, digits = 2)) %>%
    pcgrr::add_ucsc_segment_link(hgname = pcgr_data[['assembly']][['hg_name']], chrom = "chrom", start = "segment_start", end = "segment_end")
  )

  chrOrder <- c(as.character(paste0('chr',c(1:22))),"chrX","chrY")
  local_df_print <- local_df
  local_df_print <- local_df_print %>%
    dplyr::select(chrom,segment_start,segment_end,segment_length_Mb,event_type,cytoband,
                  LogR,sample_id, ensembl_gene_id,symbol,ensembl_transcript_id,transcript_start,
                  transcript_end,transcript_overlap_percent,name,biotype,disgenet_cui,
                  tsgene,p_oncogene,intogen_drivers,chembl_compound_id,gencode_gene_biotype,
                  gencode_tag,gencode_release) %>%
    dplyr::mutate(chrom = factor(chrom, levels=chrOrder)) %>%
    dplyr::arrange(chrom)

  local_df_print_sorted <- NULL
  for(chrom in chrOrder){
    if(nrow(local_df_print[!is.na(local_df_print$chrom) & local_df_print$chrom == chrom,]) > 0){
      chrom_regions <- local_df_print[local_df_print$chrom == chrom,]
      chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(segment_start, segment_end)),]
      local_df_print_sorted <- rbind(local_df_print_sorted, chrom_regions_sorted)
    }
  }

  avg_transcript_overlap <- as.data.frame(
    local_df %>%
      dplyr::filter(biotype == 'protein_coding') %>%
      dplyr::group_by(segmentID, symbol) %>%
      dplyr::summarise(MEAN_TRANSCRIPT_CNA_OVERLAP = mean(transcript_overlap_percent),
                       TRANSCRIPTS = paste0(ensembl_transcript_id, collapse=", ")) %>%
      dplyr::rename(SYMBOL = symbol) %>%
      dplyr::mutate(MEAN_TRANSCRIPT_CNA_OVERLAP = round(MEAN_TRANSCRIPT_CNA_OVERLAP, digits=2))
  )

  local_df <- dplyr::select(local_df, -ensembl_transcript_id) %>%
    dplyr::filter(biotype == 'protein_coding') %>%
    dplyr::distinct() %>%
    dplyr::mutate(symbol = as.character(symbol)) %>%
    dplyr::rename(CHEMBL_COMPOUND_ID = chembl_compound_id, SYMBOL = symbol) %>%
    dplyr::mutate(VAR_ID = as.character(rep(1:nrow(.)))) %>%
    pcgrr::annotate_variant_link(vardb = 'ANTINEOPHARMA', pcgr_data = pcgr_data) %>%
    dplyr::rename(ONCOGENE = p_oncogene,
                  TUMOR_SUPPRESSOR = tsgene,
                  ENTREZ_ID = entrezgene,
                  CHROMOSOME = chrom,
                  GENENAME = name,
                  TARGETED_DRUGS = ANTINEOPHARMALINK,
                  SEGMENT_LENGTH_MB = segment_length_Mb,
                  SEGMENT = segment_link,
                  TRANSCRIPT_OVERLAP = transcript_overlap_percent) %>%
    dplyr::mutate(ENTREZ_ID = as.character(ENTREZ_ID)) %>%
    dplyr::left_join(pcgr_data[['kegg']][['pathway_links']], by=c("ENTREZ_ID" = "gene_id")) %>%
    dplyr::rename(KEGG_PATHWAY = kegg_pathway_urls)
  entrezgene_annotation_links <- pcgrr::generate_annotation_link(local_df,
                                                                 group_by_var = "VAR_ID",
                                                                 url_prefix = "http://www.ncbi.nlm.nih.gov/gene/",
                                                                 link_key_var = "ENTREZ_ID",
                                                                 link_display_var = "GENENAME")
  local_df <- dplyr::left_join(local_df, dplyr::rename(entrezgene_annotation_links,GENE_NAME = link),by=c("VAR_ID"))

  local_df <- local_df %>%
    dplyr::select(segmentID, CHROMOSOME, SYMBOL, GENE_NAME, KEGG_PATHWAY,
                  TUMOR_SUPPRESSOR, ONCOGENE, TARGETED_DRUGS,SEGMENT_LENGTH_MB,
                  SEGMENT, biotype,LogR) %>%
    dplyr::distinct() %>%
    dplyr::left_join(avg_transcript_overlap,by=c("segmentID","SYMBOL"))

  targeted_drugs <- dplyr::select(local_df, SYMBOL, TARGETED_DRUGS) %>%
    dplyr::filter(!is.na(TARGETED_DRUGS)) %>%
    dplyr::distinct()

  n_cna_loss <- dplyr::filter(cna_segments, LogR <= logR_homdel) %>% nrow()
  n_cna_gain <- dplyr::filter(cna_segments, LogR >= logR_gain) %>% nrow()
  cna_segments_filtered <- data.frame()
  cna_segments_filtered <- cna_segments %>%
    dplyr::filter(LogR >= logR_gain | LogR <= logR_homdel) %>%
    dplyr::arrange(desc(LogR))
  rlogging::message(paste0("Detected ",nrow(cna_segments_filtered)," segments subject to amplification/deletion (",n_cna_loss," deletions, ",n_cna_gain," gains according to user-defined log(2) ratio thresholds)"))


  onco_ts_sets <- list()
  onco_ts_sets[['oncogene_gain']] <- data.frame()
  onco_ts_sets[['oncogene_gain']] <- dplyr::filter(local_df, ONCOGENE == T & MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct & LogR >= logR_gain)
  onco_ts_sets[['tsgene_loss']] <- data.frame()
  onco_ts_sets[['tsgene_loss']] <- dplyr::filter(local_df, TUMOR_SUPPRESSOR == T & MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct  & LogR <= logR_homdel)
  for(t in c('oncogene_gain','tsgene_loss')){
    if(nrow(onco_ts_sets[[t]]) > 0){
      onco_ts_sets[[t]] <- onco_ts_sets[[t]] %>%
        dplyr::select(-c(segmentID, TUMOR_SUPPRESSOR, ONCOGENE,TARGETED_DRUGS)) %>%
        dplyr::distinct() %>%
        dplyr::left_join(targeted_drugs, by="SYMBOL") %>%
        dplyr::arrange(TARGETED_DRUGS, LogR) %>%
        dplyr::select(CHROMOSOME, SYMBOL, GENE_NAME, KEGG_PATHWAY, TARGETED_DRUGS, dplyr::everything()) %>%
        dplyr::mutate(MEAN_TRANSCRIPT_CNA_OVERLAP = paste0(MEAN_TRANSCRIPT_CNA_OVERLAP,"%")) %>%
        dplyr::mutate(CNA_TYPE = dplyr::if_else(t == 'oncogene_gain','gain','loss'))
      if(t == 'oncogene_gain'){
        rlogging::message(paste0("Detected ",nrow(onco_ts_sets[[t]])," proto-oncogene(s) subject to amplification (log(2) ratio >= ",logR_gain,"): ",paste0(unique(onco_ts_sets[[t]]$SYMBOL),collapse=", ")))
      }else{
        rlogging::message(paste0("Detected ",nrow(onco_ts_sets[[t]])," tumor suppressor gene(s) subject to homozygous deletions (log(2) ratio <= ",logR_homdel,"): ",paste0(unique(onco_ts_sets[[t]]$SYMBOL),collapse=", ")))
      }
    }else{
      if(t == 'tsgene_loss'){
        rlogging::message(paste0("Detected 0 tumor suppressor genes subject to homozygous deletion (log(2) ratio <= ",logR_homdel))
      }else{
        rlogging::message(paste0("Detected 0 proto-oncogenes subject to amplification (log(2) ratio >= ",logR_gain))
      }
    }
  }

  biomarker_hits_any <- pcgrr::get_clinical_associations_cna(onco_ts_sets, pcgr_data, pcgr_config, tumor_type_specificity = 'any_tumortype')
  biomarker_hits_specific <- pcgrr::get_clinical_associations_cna(onco_ts_sets, pcgr_data, pcgr_config, tumor_type_specificity = 'specific_tumortype')

  pcg_report_cna[['eval']] <- T
  pcg_report_cna[['clinical_evidence_item']][['any_tumortype']] <- biomarker_hits_any[['clinical_evidence_item']]
  pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']] <- biomarker_hits_specific[['clinical_evidence_item']]
  pcg_report_cna[['variant_set']][['cna_print']] <- local_df_print_sorted
  pcg_report_cna[['variant_statistic']][['n_cna_gain']] <- n_cna_gain
  pcg_report_cna[['variant_statistic']][['n_cna_loss']] <- n_cna_loss
  pcg_report_cna[['variant_display']][['segment']] <- cna_segments_filtered
  pcg_report_cna[['variant_display']][['oncogene_gain']] <- onco_ts_sets[['oncogene_gain']]
  pcg_report_cna[['variant_display']][['tsgene_loss']] <- onco_ts_sets[['tsgene_loss']]
  pcg_report_cna[['variant_display']][['biomarker']] <- biomarker_hits_specific$cna_biomarkers
  pcg_report_cna <- pcgrr::assign_tier1_tier2_acmg_cna(pcg_report_cna)


  return(pcg_report_cna)
}



# annotate_facets_cna <- function(facets_cna_input_fname, facets_cna_output_fname, pcgr_data, sample_name, transcript_overlap_pct = 50){
#
#   assembly <- 'hg38'
#   ucsc_browser_prefix <- 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position='
#   if(genome_assembly == 'grch37'){
#     assembly <- 'hg19'
#     ucsc_browser_prefix <- 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position='
#   }
#
#   rlogging::message('------')
#   rlogging::message(paste0("Annotating copy number segment file from FACETS -  ",facets_cna_input_fname))
#   cna_df_raw <- read.table(file=facets_cna_input_fname,header = T,stringsAsFactors = F,comment.char="", quote="")
#   for(col in c('Cellular_Fraction','Total_CN','Minor_CN','Start','End','Segment_Mean','Chromosome')){
#     if(!(col %in% colnames(cna_df_raw))){
#       rlogging::stop('Missing required column ',col,' in cna input')
#     }
#   }
#
#   cna_df_raw <- dplyr::rename(cna_df_raw, chromosome = Chromosome, LogR = Segment_Mean, segment_start = Start, segment_end = End, cellular_fraction = Cellular_Fraction, total_cn = Total_CN, minor_cn = Minor_CN) %>% dplyr::distinct()
#   cna_df_raw$chromosome <- stringr::str_replace(cna_df_raw$chromosome,"^chr","")
#
#   ## VALIDATE INPUT CHROMOSOMES
#   cna_df <- pcgrr::get_valid_chromosomes(cna_df_raw, chromosome_column = 'chromosome', bsg = pcgr_data[['assembly']][['seq']])
#   cna_df <- pcgrr::get_valid_chromosome_segments(cna_df, pcgr_data[['assembly']][['grch_name']], bsg = pcgr_data[['assembly']][['seq']])
#   cna_df <- cna_df %>% dplyr::filter(!is.na(LogR))
#   cna_df$LogR <- round(as.numeric(cna_df$LogR),digits=3)
#   cna_df$segmentID <- paste0(cna_df$chromosome,":",cna_df$segment_start,":",cna_df$segment_end)
#   cna_df$sample_id <- sample_name
#
#   ## MAKE GRANGES OBJECT OF INPUT
#   cna_gr <- GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T, seqinfo = pcgr_data[['assembly']][['seqinfo']], seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end', ignore.strand = T, starts.in.df.are.0based = T)
#   cytoband_df <- pcgrr::get_cna_cytoband(cna_gr, pcgr_data[['genomic_ranges']][['cytoband']])
#   cna_df <- dplyr::left_join(cna_df, cytoband_df,by="segmentID")
#
#   cna_segments <- cna_df
#   #cna_segments$segment_link <- paste0("<a href='",paste0(ucsc_browser_prefix,paste0(cna_segments$chromosome,':',cna_segments$segment_start,'-',cna_segments$segment_end)),"' target=\"_blank\">",paste0(cna_segments$chromosome,':',cna_segments$segment_start,'-',cna_segments$segment_end),"</a>")
#   cna_segments$segment_length_Mb <- round((as.numeric((cna_segments$segment_end - cna_segments$segment_start)/1000000)),digits = 4)
#   #cna_segments <- dplyr::rename(cna_segments, SEGMENT_LENGTH_MB = segment_length_Mb, SEGMENT = segment_link)
#   #cna_segments <- dplyr::select(cna_segments, SEGMENT, SEGMENT_LENGTH_MB, cytoband, LogR, event_type) %>% dplyr::distinct()
#
#   cna_gr <- GenomicRanges::makeGRangesFromDataFrame(cna_df, keep.extra.columns = T, seqinfo = pcgr_data[['assembly']][['seqinfo']], seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end', ignore.strand = T, starts.in.df.are.0based = T)
#
#   hits <- GenomicRanges::findOverlaps(cna_gr, pcgr_data[['genomic_ranges']][['gencode_genes']], type="any", select="all")
#   ranges <- pcgr_data[['genomic_ranges']][['gencode_genes']][subjectHits(hits)]
#   mcols(ranges) <- c(mcols(ranges),mcols(cna_gr[queryHits(hits)]))
#
#   local_df <- as.data.frame(mcols(ranges))
#   local_df$segment_start <- start(ranges(cna_gr[queryHits(hits)]))
#   local_df$segment_end <- end(ranges(cna_gr[queryHits(hits)]))
#   local_df$segment_length_Mb <- round((as.numeric((local_df$segment_end - local_df$segment_start)/1000000)),digits = 4)
#
#   local_df$transcript_start <- start(ranges)
#   local_df$transcript_end <- end(ranges)
#   local_df$chrom <- as.character(seqnames(ranges))
#   local_df <- as.data.frame(local_df %>% dplyr::rowwise() %>% dplyr::mutate(transcript_overlap_percent = round(as.numeric((min(transcript_end,segment_end) - max(segment_start,transcript_start)) / (transcript_end - transcript_start)) * 100, digits = 2)))
#   #local_df$segment_link <- paste0("<a href='",paste0(ucsc_browser_prefix,paste0(local_df$chrom,':',local_df$segment_start,'-',local_df$segment_end)),"' target=\"_blank\">",paste0(local_df$chrom,':',local_df$segment_start,'-',local_df$segment_end),"</a>")
#
#   local_df_print <- local_df
#   local_df_print <- dplyr::select(local_df_print,chrom,segment_start,segment_end,segment_length_Mb,sample_id,event_type,cytoband,LogR,cellular_fraction,total_cn,minor_cn,ensembl_gene_id,symbol,ensembl_transcript_id,transcript_start,transcript_end,transcript_overlap_percent,name,biotype,tsgene,p_oncogene,chembl_compound_id,gencode_gene_biotype,gencode_tag,gencode_release)
#
#   chrOrder <- c(as.character(paste0('chr',c(1:22))),"chrX","chrY")
#   local_df_print$chrom <- factor(local_df_print$chrom, levels=chrOrder)
#   local_df_print <- local_df_print[order(local_df_print$chrom),]
#   local_df_print$segment_start <- as.integer(local_df_print$segment_start)
#   local_df_print$segment_end <- as.integer(local_df_print$segment_end)
#
#   local_df_print_sorted <- NULL
#   for(chrom in chrOrder){
#     if(nrow(local_df_print[!is.na(local_df_print$chrom) & local_df_print$chrom == chrom,]) > 0){
#       chrom_regions <- local_df_print[local_df_print$chrom == chrom,]
#       chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(segment_start, segment_end)),]
#       local_df_print_sorted <- rbind(local_df_print_sorted, chrom_regions_sorted)
#     }
#   }
#   write.table(local_df_print_sorted,file=facets_cna_output_fname,col.names = T,row.names = F,quote = F)
#   system(paste0('gzip ',facets_cna_output_fname))
#
#   homdel <- dplyr::filter(local_df_print_sorted, !is.na(minor_cn) & minor_cn == 0 & total_cn == 0)
#   ampl <- dplyr::filter(local_df_print_sorted, ((is.na(minor_cn) & total_cn >= 5) | (!is.na(minor_cn) & total_cn - minor_cn >= 5)))
#
#   homdel_ampl <- data.frame()
#   if(nrow(homdel) > 0 | nrow(ampl) > 0){
#     homdel_ampl <- dplyr::bind_rows(homdel,ampl)
#   }
#   return(homdel_ampl)
#
# }
#
