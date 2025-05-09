#' Function that parses and loads reference data from files in the
#' assembly-specific PCGR bundle directory
#'
#' @param pcgr_db_assembly_dir Assembly-specific root path for data bundle
#' @param genome_assembly grch37/grch38
#' @export
#'
load_reference_data <- function(
    pcgr_db_assembly_dir = NULL,
    genome_assembly = "grch38") {

  pcgr_ref_data <- list()

  log4r_info(paste0(
    "Loading reference datasets - genome assembly: ", genome_assembly))

  pcgr_ref_data[["assembly"]] <- list()
  pcgr_ref_data[["assembly"]][["grch_name"]] <- genome_assembly
  pcgr_ref_data[["assembly"]][["hg_name"]] <- "hg19"
  pcgr_ref_data[["assembly"]][["ref_genome"]] <- "BSgenome.Hsapiens.UCSC.hg19"
  if (genome_assembly == "grch38") {
    pcgr_ref_data[["assembly"]][["grch_name"]] <- genome_assembly
    pcgr_ref_data[["assembly"]][["hg_name"]] <- "hg38"
    pcgr_ref_data[["assembly"]][["ref_genome"]] <- "BSgenome.Hsapiens.UCSC.hg38"
  }

  bsgenome_obj <- pcgrr::get_genome_obj(genome_assembly)
  genome_grch2hg <- c("grch38" = "hg38", "grch37" = "hg19")
  pcgr_ref_data[['assembly']][['seqinfo']] <-
    GenomeInfoDb::Seqinfo(
      seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(bsgenome_obj)),
      seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bsgenome_obj)),
      genome = genome_grch2hg[genome_assembly])
  pcgr_ref_data[['assembly']][['bsg']] <- bsgenome_obj

  chromsizes_fname <- file.path(
    pcgr_db_assembly_dir,
    paste0("chromsize.", genome_assembly, ".tsv"))
  check_file_exists(chromsizes_fname)
  pcgr_ref_data[['assembly']][['chrom_coordinates']] <-
    as.data.frame(
      readr::read_tsv(
        chromsizes_fname, show_col_types = F,
        na = c(".")
      )
    ) |>
    dplyr::mutate(
      genome_end = cumsum(.data$length)) |>
    dplyr::mutate(
      genome_start = .data$genome_end - .data$length
    ) |>
    dplyr::select(
      c("chrom", "genome_start", "genome_end",
      "length"), dplyr::everything())



  pcgr_ref_data[['vcf_infotags']] <- data.frame()
  for(t in c('vep','other')) {
    infotag_fname <- file.path(
      pcgr_db_assembly_dir,
      paste0("vcf_infotags_", t, ".tsv"))
    check_file_exists(infotag_fname)
    pcgr_ref_data[['vcf_infotags']] <-
      pcgr_ref_data[['vcf_infotags']] |>
      dplyr::bind_rows(
        as.data.frame(
          readr::read_tsv(
            infotag_fname, show_col_types = F,
            na = c(".")
          )
        )
    )
  }
  for(cat in c('tcga','clinvar','gwas','gnomad_non_cancer',
               'dbmts','dbnsfp','panel_of_normals')) {
    vcfanno_fname <- file.path(
      pcgr_db_assembly_dir,"variant","vcf",cat,
      paste0(cat,".vcfanno.vcf_info_tags.txt"))

    raw_lines <- readLines(vcfanno_fname)
    for(l in raw_lines) {
      if (startsWith(l,"##INFO")) {
        tag <- stringr::str_replace(
          stringr::str_match(l,"ID=[A-Za-z|_]{1,}")[,1],
          "ID=","")
        number <- NA
        if (stringr::str_detect(l, "Number=1,")) {
          number <- 1
        }
        if (stringr::str_detect(l, "Number=0,")) {
          number <- 0
        }
        type <- "String"
        if (stringr::str_detect(l, "Type=Integer,")) {
          type <- "Integer"
        }
        if (stringr::str_detect(l, "Type=Float,")) {
          type <- "Float"
        }
        if (stringr::str_detect(l, "Type=Flag,")) {
          type <- "Flag"
        }
        description <-
          stringr::str_replace_all(
            stringr::str_match(
              l,
              "Description=.+$")[,1],
            "Description=\\\"|\\\">","")

        category <- "pcgr_cpsr"
        if (cat == "dbmts" | cat == "gnomad_non_cancer") {
          category <- "cpsr"
        }
        if (cat == "panel_of_normals") {
          category <- "pcgr"
        }
        df <- data.frame(
          'tag' = tag,
          'number' = number,
          'type' = type,
          'description' = description,
          'category' = category
        )
        pcgr_ref_data[['vcf_infotags']] <-
          pcgr_ref_data[['vcf_infotags']] |>
          dplyr::bind_rows(df)
      }
    }
  }

  ####---- 1. Gene annotations----####

  pcgr_ref_data[["gene"]] <- list()
  pcgr_ref_data[["gene"]][["panel"]] <- data.frame()
  pcgr_ref_data[["gene"]][["cpg"]] <- data.frame()
  pcgr_ref_data[['gene']][['gene_xref']] <- data.frame()
  pcgr_ref_data[['gene']][['otp_rank']] <- data.frame()

  cpg_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_cpg",
    "gene_cpg.tsv.gz"
  )
  check_file_exists(cpg_tsv_fname)
  pcgr_ref_data[['gene']][['cpg']] <- as.data.frame(
    readr::read_tsv(cpg_tsv_fname, show_col_types = F,
                    comment = "", na = c(".",""))
  ) |>
    dplyr::mutate(
      entrezgene = as.character(.data$entrezgene)
    )
  colnames(pcgr_ref_data[['gene']][['cpg']]) <-
    toupper(colnames(pcgr_ref_data[['gene']][['cpg']]))


  panels_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_virtual_panel",
    "gene_virtual_panel.tsv.gz"
  )
  check_file_exists(panels_tsv_fname)
  pcgr_ref_data[['gene']][['panel']] <- as.data.frame(
    readr::read_tsv(panels_tsv_fname, show_col_types = F,
                    comment = "", na = c("","."))
  ) |>
    dplyr::mutate(
      entrezgene = as.character(.data$entrezgene)
    )
  colnames(pcgr_ref_data[['gene']][['panel']]) <-
    toupper(colnames(pcgr_ref_data[['gene']][['panel']]))


  gene_xref_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_transcript_xref",
    "gene_transcript_xref.tsv.gz"
  )
  check_file_exists(gene_xref_tsv_fname)

  pcgr_ref_data[['gene']][['transcript_biotype']] <- as.data.frame(
    readr::read_tsv(
      gene_xref_tsv_fname, show_col_types = F)) |>
    dplyr::select(
      c("chrom",
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "gencode_transcript_biotype",
        "gene_biotype"
      )
    ) |>
    dplyr::distinct()

  colnames(pcgr_ref_data[['gene']][['transcript_biotype']]) <-
    toupper(colnames(pcgr_ref_data[['gene']][['transcript_biotype']]))

  gene_index_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_transcript_xref",
    "gene_index.tsv.gz"
  )
  check_file_exists(gene_index_tsv_fname)

  pcgr_ref_data[['gene']][['index']] <- as.data.frame(
    readr::read_tsv(
      gene_index_tsv_fname, show_col_types = F))

  otp_rank_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_transcript_xref",
    "otp_rank.tsv.gz"
  )
  check_file_exists(otp_rank_tsv_fname)

  pcgr_ref_data[['gene']][['otp_rank']] <- as.data.frame(
    readr::read_tsv(
      otp_rank_tsv_fname, show_col_types = F,
      na = c('.'))) |>
    dplyr::filter(!is.na(.data$entrezgene)) |>
    dplyr::mutate(entrezgene = as.character(.data$entrezgene)) |>
    dplyr::distinct()

  colnames(pcgr_ref_data[['gene']][['otp_rank']]) <-
    toupper(colnames(pcgr_ref_data[['gene']][['otp_rank']]))

  pcgr_ref_data[['gene']][['gene_xref']] <- as.data.frame(
    readr::read_tsv(
      gene_xref_tsv_fname, na = c("."),
      show_col_types = F, guess_max = 200000)) |>
    dplyr::select(
      c("ensembl_gene_id",
      "gene_biotype",
      "symbol",
      "entrezgene",
      "name",
      "driver",
      "driver_support",
      "tsg",
      "tsg_support",
      "tsg_rank",
      "oncogene",
      "oncogene_support",
      "oncogene_rank",
      "cancergene_evidence")
    ) |>
    dplyr::rename(
      genename = "name"
    ) |>
    dplyr::mutate(
      entrezgene = as.character(.data$entrezgene)
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      oncogene = dplyr::if_else(
        is.na(.data$oncogene),
        FALSE,
        as.logical(.data$oncogene)
      )
    ) |>
    dplyr::mutate(
      tsg = dplyr::if_else(
        is.na(.data$tsg),
        FALSE,
        as.logical(.data$tsg)
      )
    ) |>
    dplyr::mutate(
      driver = dplyr::if_else(
        is.na(.data$driver),
        FALSE,
        as.logical(.data$driver)
      )
    )

  colnames(pcgr_ref_data[['gene']][['gene_xref']]) <-
    toupper(colnames(pcgr_ref_data[['gene']][['gene_xref']]))

  ####--- 2. Variant annotations----####
  pcgr_ref_data[['variant']] <- list()

  #####--A. Sites of pathogenic variants (ClinVar)--#####
  # 1. sites (codons) of pathogenic/likely pathogenic variants
  clinvar_sites_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "variant", "tsv",
      "clinvar", "clinvar_sites.tsv.gz"
    )
  check_file_exists(clinvar_sites_tsv_fname)
  pcgr_ref_data[['variant']][['clinvar_sites']] <- as.data.frame(
    readr::read_tsv(
      clinvar_sites_tsv_fname, show_col_types = F)) |>
    dplyr::mutate(entrezgene = as.character(.data$entrezgene))
  colnames(pcgr_ref_data[['variant']][['clinvar_sites']]) <-
    toupper(colnames(pcgr_ref_data[['variant']][['clinvar_sites']]))

  #####--B. Oncogenic variants (ClinVar)--#####
  # known oncogenic variants
  clinvar_oncogenic_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "variant", "tsv",
      "clinvar", "clinvar_oncogenic.tsv.gz"
    )
  check_file_exists(clinvar_oncogenic_tsv_fname)
  pcgr_ref_data[['variant']][['clinvar_oncogenic']] <- as.data.frame(
    readr::read_tsv(
      clinvar_oncogenic_tsv_fname, show_col_types = F)) |>
    dplyr::mutate(entrezgene = as.character(.data$entrezgene)) |>
    dplyr::mutate(alteration = dplyr::if_else(
      is.na(.data$hgvsp) & !is.na(.data$hgvs_c),
      as.character(.data$hgvs_c),
      .data$hgvsp
    )) |>
    dplyr::select(-c("codon","trait","var_id")) |>
    dplyr::select(c("symbol","alteration","molecular_consequence",
                    "oncogenicity","review_status_oncogenicity"),
                  dplyr::everything())

  colnames(pcgr_ref_data[['variant']][['clinvar_oncogenic']]) <-
    toupper(colnames(pcgr_ref_data[['variant']][['clinvar_oncogenic']]))

  #####--C. Gene-level variant statistics (ClinVar)--#####
  clinvar_gene_varstats_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "variant", "tsv",
      "clinvar", "clinvar_gene_varstats.tsv.gz"
    )
  check_file_exists(clinvar_gene_varstats_tsv_fname)
  pcgr_ref_data[['variant']][['clinvar_gene_stats']] <- as.data.frame(
    readr::read_tsv(
      clinvar_gene_varstats_tsv_fname, show_col_types = F)) |>
    dplyr::mutate(entrezgene = as.character(.data$entrezgene))
  colnames(pcgr_ref_data[['variant']][['clinvar_gene_stats']]) <-
    toupper(colnames(pcgr_ref_data[['variant']][['clinvar_gene_stats']]))



  #####--D. GWAS variants--#####
  gwas_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "variant", "tsv", "gwas", "gwas.tsv.gz"
    )
  check_file_exists(gwas_tsv_fname)
  pcgr_ref_data[['variant']][['gwas']] <- as.data.frame(
    readr::read_tsv(
      gwas_tsv_fname,
      na = c("."),
      show_col_types = F)) |>
    dplyr::filter(!is.na(.data$primary_site))

  colnames(pcgr_ref_data[['variant']][['gwas']]) <-
    toupper(colnames(pcgr_ref_data[['variant']][['gwas']]))


  pcgr_ref_data[['variant']][['varstats']] <- list()
  ## Get variant statistics
  for(vardb in c('clinvar','gwas','tcga',
                 'gnomad_non_cancer','dbmts',
                 'dbnsfp')) {
    varstats_fname <-
      file.path(
        pcgr_db_assembly_dir, "variant", "vcf", vardb,
        paste0(vardb,".vcf_varstats.tsv")
      )

    if (file.exists(varstats_fname)) {
      pcgr_ref_data[['variant']][['varstats']][[vardb]] <-
        as.data.frame(
          readr::read_tsv(
            varstats_fname, show_col_types = F))
    }

  }



  ####---- 3. Phenotype ontologies ----####
  pcgr_ref_data[['phenotype']] <- list()
  phenOncoX_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "phenotype",
      "tsv", "phenotype_onco.tsv.gz"
    )
  check_file_exists(phenOncoX_tsv_fname)
  pcgr_ref_data[['phenotype']][['oncotree']] <- as.data.frame(
    readr::read_tsv(
      phenOncoX_tsv_fname, show_col_types = F,
      na = c(".","NA"))
  )
  pcgr_ref_data[['phenotype']][['cancer_groups']] <-
    unique(pcgr_ref_data[['phenotype']][['oncotree']]$primary_site)
  pcgr_ref_data[['phenotype']][['cancer_groups']] <-
    pcgr_ref_data[['phenotype']][['cancer_groups']][
      !is.na(pcgr_ref_data[['phenotype']][['cancer_groups']])]

  colnames(pcgr_ref_data[['phenotype']][['oncotree']]) <-
    toupper(colnames(pcgr_ref_data[['phenotype']][['oncotree']]))

  umls_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "phenotype", "tsv",
      "phenotype_umls.tsv.gz"
    )
  check_file_exists(umls_tsv_fname)
  pcgr_ref_data[['phenotype']][['umls']] <- as.data.frame(
    readr::read_tsv(
      umls_tsv_fname, show_col_types = F,
      na = ".")
  )
  colnames(pcgr_ref_data[['phenotype']][['umls']]) <-
    toupper(colnames(pcgr_ref_data[['phenotype']][['umls']]))

  ####--- 4. MSI predictor----####
  msi_model_rds <-
    file.path(
      pcgr_db_assembly_dir, "misc", "other",
      "msi_classification",
      "tcga_msi_classifier.rds"
    )
  check_file_exists(msi_model_rds)
  pcgr_ref_data[['msi']] <-
    readRDS(msi_model_rds)


  pcgr_ref_data[['misc']] <- list()
  ####--- 5. Miscellaneous----####
  for(elem in c('tmb',
                'mutational_signature',
                'pathway',
                'oncogenicity',
                'hotspot',
                'protein_domain')) {

    fname_misc <- file.path(
      pcgr_db_assembly_dir, "misc", "tsv", elem,
      paste0(elem,".tsv.gz")
    )
    if(elem == "oncogenicity"){
      fname_misc <- file.path(
        pcgr_db_assembly_dir, "misc", "tsv", elem,
        paste0(elem,".tsv")
      )
    }
    check_file_exists(fname_misc)
    pcgr_ref_data[['misc']][[elem]] <- as.data.frame(
      readr::read_tsv(
        fname_misc, show_col_types = F,
        na = ".", guess_max = 11000)
    )
    colnames(pcgr_ref_data[['misc']][[elem]]) <-
      toupper(colnames(pcgr_ref_data[['misc']][[elem]]))

  }

  tmp = pcgr_ref_data[['misc']][['pathway']]
  pcgr_ref_data[['misc']][['pathway']] <- list()
  pcgr_ref_data[['misc']][['pathway']][['long']] <- tmp
  pcgr_ref_data[['misc']][['pathway']][['wide']] <- as.data.frame(
    tmp |>
    dplyr::group_by(.data$GENE_ID) |>
    dplyr::summarise(LINK = paste(.data$URL_HTML, collapse = ", ")))


  ####--- 6. Drugs ----####
  pcgr_ref_data[['drug']] <- list()
  drug_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "drug",
      "tsv", "drug_targeted.tsv.gz"
    )
  check_file_exists(drug_tsv_fname)
  pcgr_ref_data[['drug']][['targeted']] <- as.data.frame(
    readr::read_tsv(drug_tsv_fname, show_col_types = F, na = ".")
  )
  colnames(pcgr_ref_data[['drug']][['targeted']]) <-
    toupper(colnames(pcgr_ref_data[['drug']][['targeted']]))

  drug_all_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "drug",
      "tsv", "drug_all.tsv.gz"
    )
  check_file_exists(drug_all_tsv_fname)
  pcgr_ref_data[['drug']][['all']] <- as.data.frame(
    readr::read_tsv(drug_all_tsv_fname, show_col_types = F, na = ".")
  )
  colnames(pcgr_ref_data[['drug']][['all']]) <-
    toupper(colnames(pcgr_ref_data[['drug']][['all']]))



  inhibitors_all <-
    dplyr::select(
      pcgr_ref_data[['drug']][['targeted']],
      c("SYMBOL",
        "QUERY_SITE",
        "DRUG_PRIMARY_SITE",
        "ATC_TREATMENT_CATEGORY",
        "ATC_LEVEL3",
        "DRUG_NAME",
        "DRUG_TYPE",
        "DRUG_MAX_PHASE_INDICATION",
        "DRUG_YEAR_FIRST_APPROVAL",
        "DRUG_ACTION_TYPE",
        "DRUG_LINK")) |>
    dplyr::rename(
      DRUG_CLASS = "ATC_LEVEL3") |>
    dplyr::filter(
        .data$ATC_TREATMENT_CATEGORY != "cancer_unclassified" &
        .data$DRUG_TYPE != "Unknown") |>
    dplyr::distinct() |>
    dplyr::mutate(
      DRUG_ACTION_TYPE = stringr::str_to_title(
        .data$DRUG_ACTION_TYPE
      )) |>
    dplyr::filter(
      .data$DRUG_ACTION_TYPE != "Other") |>
    dplyr::arrange(
      .data$SYMBOL,
      dplyr::desc(.data$DRUG_MAX_PHASE_INDICATION),
      dplyr::desc(.data$DRUG_YEAR_FIRST_APPROVAL)) |>
    dplyr::select(-c("DRUG_YEAR_FIRST_APPROVAL",
                     "ATC_TREATMENT_CATEGORY",
                     "DRUG_TYPE")) |>
    dplyr::distinct()

  pcgr_ref_data[['drug']][['inhibitors_on_label']] <-
    inhibitors_all |>
    dplyr::filter(
       .data$DRUG_MAX_PHASE_INDICATION > 2) |>
    dplyr::filter(
      .data$QUERY_SITE == .data$DRUG_PRIMARY_SITE &
      .data$QUERY_SITE != "Any" &
        .data$QUERY_SITE != "Other/Unknown") |>
    dplyr::group_by(
      .data$SYMBOL,
      .data$QUERY_SITE,
      .data$DRUG_CLASS) |>
    dplyr::summarise(
      DRUG_NAME = paste(.data$DRUG_NAME, collapse = "|"),
      DRUG_LINK = paste(.data$DRUG_LINK, collapse = ", "),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      .data$SYMBOL,
      dplyr::desc(.data$n)
    ) |>
    dplyr::mutate(
      DRUG_NAME = paste(
        .data$DRUG_CLASS, .data$DRUG_NAME, sep = ":"),
      DRUG_LINK = paste0(
        "<b>",.data$DRUG_CLASS,": </b>", .data$DRUG_LINK),
    ) |>
    ## ignore proteasome/tubulin/CYP targets -
    dplyr::filter(
      !stringr::str_detect(
        .data$SYMBOL, "^(CYP|PSM|TUB)")
    ) |>
    dplyr::group_by(
      .data$SYMBOL,
      .data$QUERY_SITE
    ) |>
    dplyr::summarise(
      TARGETED_INHIBITORS2 = paste(.data$DRUG_NAME, collapse=", "),
      TARGETED_INHIBITORS = paste(.data$DRUG_LINK, collapse=", "),
      .groups = "drop"
    )

  pcgr_ref_data[['drug']][['inhibitors_any_label']] <-
    as.data.frame(
    inhibitors_all |>
    dplyr::filter(
      .data$DRUG_MAX_PHASE_INDICATION > 2) |>
    dplyr::filter(
        .data$QUERY_SITE != "Any") |>
    dplyr::group_by(
      .data$SYMBOL,
      .data$DRUG_NAME,
      .data$DRUG_LINK,
      .data$DRUG_CLASS) |>
    dplyr::summarise(
      DRUG_PRIMARY_SITE = paste(
        sort(unique(.data$DRUG_PRIMARY_SITE)), collapse=","),
      .groups = "drop"
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      DRUG_LINK = paste0(
        .data$DRUG_LINK, " (", .data$DRUG_PRIMARY_SITE, ")"),
      DRUG_NAME = paste0(
        .data$DRUG_NAME, " (", .data$DRUG_PRIMARY_SITE, ")")
    ) |>
    dplyr::group_by(
      .data$SYMBOL,
      .data$DRUG_CLASS) |>
    dplyr::summarise(
      DRUG_NAME = paste(.data$DRUG_NAME, collapse="|"),
      DRUG_LINK = paste(.data$DRUG_LINK, collapse=", "),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$SYMBOL,
      dplyr::desc(.data$n)
    ) |>
    dplyr::mutate(
      DRUG_NAME = paste(
        .data$DRUG_CLASS,
        .data$DRUG_NAME,
        sep = ":"),
      DRUG_LINK = paste0(
        "<b>",
        .data$DRUG_CLASS,
        ": </b>",
        .data$DRUG_LINK),
    ) |>
    ## ignore proteasome/tubulin/CYP targets -
    dplyr::filter(
      !stringr::str_detect(
        .data$SYMBOL, "^(CYP|PSM|TUB)")
    ) |>
    dplyr::group_by(
      .data$SYMBOL,
    ) |>
    dplyr::summarise(
      TARGETED_INHIBITORS_ALL2 = paste(
        .data$DRUG_NAME, collapse=", "),
      TARGETED_INHIBITORS_ALL = paste(
        .data$DRUG_LINK, collapse=", "),
      .groups = "drop"
    )
  )


  ####--- 7. Biomarkers ----####
  pcgr_ref_data[['biomarker']] <- list()
  for(elem in c('clinical','variant','literature')) {
    pcgr_ref_data[['biomarker']][[elem]] <- data.frame()
    for(db in c('cgi','civic')) {
      fname <-
        file.path(
          pcgr_db_assembly_dir, "biomarker", "tsv",
          paste0(db,".", elem,".tsv.gz")
        )
      check_file_exists(fname)
      bm_data <- as.data.frame(
        readr::read_tsv(
          fname, show_col_types = F,
          na = ".", guess_max = 10000))
      if ("source_id" %in% colnames(bm_data)) {
        bm_data <- bm_data |>
          dplyr::mutate(
            source_id = as.character(.data$source_id))
      }

      if ('entrezgene' %in% colnames(bm_data)) {
        bm_data <- bm_data |>
          dplyr::mutate(
            entrezgene = as.character(.data$entrezgene))
      }
      if ('variant_id' %in% colnames(bm_data)) {
        bm_data <- bm_data |>
          dplyr::mutate(
            variant_id = as.character(.data$variant_id))
      }

      pcgr_ref_data[['biomarker']][[elem]] <- dplyr::bind_rows(
        pcgr_ref_data[['biomarker']][[elem]],
        bm_data
      )
    }
    colnames(pcgr_ref_data[['biomarker']][[elem]]) <-
      toupper(colnames(pcgr_ref_data[['biomarker']][[elem]]))

    # if(elem == 'variant'){
    #   pcgr_ref_data[['biomarker']][[elem]] <-
    #     pcgr_ref_data[['biomarker']][[elem]] |>
    #     dplyr::mutate(VARIANT_NAME_PRIMARY = stringr::str_replace(
    #       .data$VARIANT_NAME_PRIMARY, " FUSION", " Fusion")
    #     ) |>
    #     dplyr::mutate(VARIANT_NAME_PRIMARY = dplyr::if_else(
    #       .data$VARIANT_NAME_PRIMARY == "NA Fusion" &
    #         !is.na(.data$GENE),
    #       paste0(.data$GENE, " Fusion"),
    #       as.character(.data$VARIANT_NAME_PRIMARY)
    #     ))
    # }

  }

  ####-- 8. Metadata ####
  pcgr_ref_data[['metadata']] <- data.frame()
  for(dtype in c('gene','gwas','hotspot','other',
                 'phenotype','biomarker','drug')) {

    fname <- file.path(
      pcgr_db_assembly_dir, ".METADATA", "tsv",
      paste0(dtype,"_metadata.tsv")
    )
    check_file_exists(fname)
    metadata_dtype <- as.data.frame(
      readr::read_tsv(fname, show_col_types = F,
                      guess_max = 100000)) |>
      dplyr::mutate(datatype = dtype) |>
      dplyr::mutate(wflow = dplyr::case_when(
        stringr::str_detect(
          .data$source_abbreviation,
          paste0(
            "^(gepa|cpg_other|maxwell2016|acmg_sf|dbmts|",
            "woods_dnarepair|gerp|tcga_pancan_2018|gwas_catalog|cpic)")) ~ "cpsr",
        stringr::str_detect(
          .data$source_abbreviation,
          paste0("^(cytoband|mitelmandb|tcga|nci|intogen|depmap|treehouse",
                 "|opentargets|dgidb|pubchem)$")) ~ "pcgr",
        TRUE ~ as.character("pcgr_cpsr")
      ))

    pcgr_ref_data[['metadata']] <-
      pcgr_ref_data[['metadata']] |>
      dplyr::bind_rows(metadata_dtype) |>
      dplyr::filter(
        .data$source_abbreviation != "foundation_one" &
          .data$source_abbreviation != "illumina"
      )
  }

  pcgr_ref_data[['expression']] <- list()
  exp_sources <- list.dirs(file.path(
    pcgr_db_assembly_dir, "expression", "tsv"), full.names = F)
  for(source in exp_sources){
    if(source != "" & source != "tcga"){
      metadata_fname <- file.path(
        pcgr_db_assembly_dir, "expression", "tsv",
        source, paste0(source, "_sample_metadata.tsv.gz"))
      if(file.exists(metadata_fname)){
        pcgr_ref_data[['expression']][[source]] <-
          suppressWarnings(readr::read_tsv(
            metadata_fname, show_col_types = F,
            guess_max = 10000, na = ".") |>
              dplyr::mutate(
                EXT_SOURCE = source) |>
              dplyr::mutate(EXT_AGE_AT_DX = as.numeric(
                .data$EXT_AGE_AT_DX))
          )
      }
    }else{
      if(source == "tcga"){
        pcgr_ref_data[['expression']][[source]] <- data.frame()
        metadata_fnames <- list.files(
          file.path(
            pcgr_db_assembly_dir, "expression", "tsv", source),
          pattern = "_sample_metadata.tsv.gz$", full.names = T)
        for(metadata_fname in metadata_fnames){
          if(file.exists(metadata_fname)){
            cohort_metadata <- suppressWarnings(
              readr::read_tsv(
                metadata_fname, show_col_types = F,
                guess_max = 10000, na = ".") |>
                dplyr::mutate(
                  EXT_SOURCE = source) |>
                dplyr::mutate(EXT_AGE_AT_DX = as.numeric(
                  .data$EXT_AGE_AT_DX))
              )

            pcgr_ref_data[['expression']][[source]] <-
              pcgr_ref_data[['expression']][[source]] |>
              dplyr::bind_rows(
                cohort_metadata)
          }
        }
      }
    }
  }

  return(pcgr_ref_data)


}
