#' Function that parses and loads reference data from files in the
#' assembly-specific PCGR bundle directories
#'
#' @param pcgr_db_assembly_dir Assembly-specific root path for data bundle
#' @param genome_assembly grch37/grch38
#' @export
#'
load_reference_data <- function(
    pcgr_db_assembly_dir = NULL,
    genome_assembly = "grch38"){

  pcgr_ref_data <- list()

  pcgr_ref_data[["assembly"]] <- list()
  pcgr_ref_data[["assembly"]][["grch_name"]] <- genome_assembly
  pcgr_ref_data[["assembly"]][["grch_name"]] <- "hg19"
  pcgr_ref_data[["assembly"]][["ref_genome"]] <- "BSgenome.Hsapiens.UCSC.hg19"
  if (genome_assembly == "grch38") {
    pcgr_ref_data[["assembly"]][["grch_name"]] <- genome_assembly
    pcgr_ref_data[["assembly"]][["grch_name"]] <- "hg38"
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


  pcgr_ref_data[['vcf_infotags']] <- data.frame()
  for(t in c('vep','other')){
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
               'dbmts','dbnsfp','panel_of_normals')){
    vcfanno_fname <- file.path(
      pcgr_db_assembly_dir,"variant","vcf",cat,
      paste0(cat,".vcfanno.vcf_info_tags.txt"))

    raw_lines <- readLines(vcfanno_fname)
    for(l in raw_lines){
      if(startsWith(l,"##INFO")){
        tag <- stringr::str_replace(
          stringr::str_match(l,"ID=[A-Za-z|_]{1,}")[,1],
          "ID=","")
        number <- NA
        if(stringr::str_detect(l, "Number=1,")){
          number <- 1
        }
        if(stringr::str_detect(l, "Number=0,")){
          number <- 0
        }
        type <- "String"
        if(stringr::str_detect(l, "Type=Integer,")){
          type <- "Integer"
        }
        if(stringr::str_detect(l, "Type=Float,")){
          type <- "Float"
        }
        if(stringr::str_detect(l, "Type=Flag,")){
          type <- "Flag"
        }
        description <-
          stringr::str_replace_all(
            stringr::str_match(
              l,
              "Description=.+$")[,1],
            "Description=\\\"|\\\">","")

        category <- "pcgr_cpsr"
        if(cat == "dbmts" | cat == "gnomad_non_cancer"){
          category <- "cpsr"
        }
        if(cat == "panel_of_normals"){
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

  ## 1. Gene annotations

  pcgr_ref_data[["gene"]] <- list()
  pcgr_ref_data[["gene"]][["panel"]] <- data.frame()
  pcgr_ref_data[["gene"]][["cpg"]] <- data.frame()
  pcgr_ref_data[['gene']][['gene_xref']] <- data.frame()

  cpg_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_cpg",
    "gene_cpg.tsv.gz"
  )
  check_file_exists(cpg_tsv_fname)
  pcgr_ref_data[['gene']][['cpg']] <- as.data.frame(
    readr::read_tsv(cpg_tsv_fname, show_col_types = F,
                    comment = "", na = c(".",""))
  )


  panels_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_virtual_panel",
    "gene_virtual_panel.tsv.gz"
  )
  check_file_exists(panels_tsv_fname)
  pcgr_ref_data[['gene']][['panel']] <- as.data.frame(
    readr::read_tsv(panels_tsv_fname, show_col_types = F,
                    comment = "", na = c("","."))
  )

  gene_xref_tsv_fname <- file.path(
    pcgr_db_assembly_dir, "gene", "tsv",
    "gene_transcript_xref",
    "gene_transcript_xref.tsv.gz"
  )
  check_file_exists(gene_xref_tsv_fname)
  pcgr_ref_data[['gene']][['gene_xref']] <- as.data.frame(
    readr::read_tsv(gene_xref_tsv_fname, show_col_types = F)) |>
    dplyr::select(
      ensembl_gene_id,
      gene_biotype,
      symbol,
      entrezgene,
      name,
      driver,
      driver_support,
      tsg,
      tsg_support,
      tsg_rank,
      oncogene,
      oncogene_support,
      oncogene_rank,
      cancergene_evidence
    ) |>
    dplyr::rename(
      genename = name
    ) |>
    dplyr::mutate(
      entrezgene = as.character(entrezgene)
    ) |>
    dplyr::filter(gene_biotype == "protein_coding") |>
    dplyr::distinct() |>
    dplyr::mutate(
      oncogene = dplyr::if_else(
        oncogene == ".",
        FALSE,
        as.logical(oncogene)
      )
    ) |>
    dplyr::mutate(
      tsg = dplyr::if_else(
        tsg == ".",
        FALSE,
        as.logical(tsg)
      )
    ) |>
    dplyr::mutate(
      driver = dplyr::if_else(
        driver == ".",
        FALSE,
        as.logical(driver)
      )
    )


  ## 2. Variant annotations

  ## GWAS
  gwas_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "variant", "tsv", "gwas", "gwas.tsv.gz"
    )
  check_file_exists(gwas_tsv_fname)
  pcgr_ref_data[['gwas']] <- as.data.frame(
    readr::read_tsv(
      gwas_tsv_fname, show_col_types = F)
  )


  ## 3. Phenotype ontologies

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


  ## 4. MSI classification
  msi_model_rds <-
    file.path(
      pcgr_db_assembly_dir, "misc", "other",
      "msi_classification",
      "msi_classification.rds"
    )
  check_file_exists(msi_model_rds)
  pcgr_ref_data[['msi']] <-
    readRDS(msi_model_rds)


  ## 5. Miscellaneous
  for(elem in c('tmb',
                'mutational_signature',
                'pathway')){

    fname_misc <- file.path(
      pcgr_db_assembly_dir, "misc", "tsv", elem,
      paste0(elem,".tsv.gz")
    )
    check_file_exists(fname_misc)
    pcgr_ref_data[[elem]] <- as.data.frame(
      readr::read_tsv(
        fname_misc, show_col_types = F,
        na = ".")
    )
  }

  tmp = pcgr_ref_data[['pathway']]
  pcgr_ref_data[['pathway']] <- list()
  pcgr_ref_data[['pathway']][['long']] <- tmp
  pcgr_ref_data[['pathway']][['wide']] <- as.data.frame(
    tmp |>
    dplyr::group_by(gene_id) |>
    dplyr::summarise(link = paste(url_html, collapse = ", ")))


  ## 6. Drugs
  drug_tsv_fname <-
    file.path(
      pcgr_db_assembly_dir, "drug",
      "tsv", "drug.tsv.gz"
    )
  check_file_exists(drug_tsv_fname)
  pcgr_ref_data[['drug']] <- as.data.frame(
    readr::read_tsv(drug_tsv_fname, show_col_types = F, na = ".")
  )


  ## 7. Biomarkers
  pcgr_ref_data[['biomarker']] <- list()
  for(elem in c('clinical','variant','literature')){
    pcgr_ref_data[['biomarker']][[elem]] <- data.frame()
    for(db in c('cgi','civic')){
      fname <-
        file.path(
          pcgr_db_assembly_dir, "biomarker", "tsv",
          paste0(db,".", elem,".tsv.gz")
        )
      check_file_exists(fname)
      data <- as.data.frame(
        readr::read_tsv(fname, show_col_types = F, na = "."))
      if("source_id" %in% colnames(data)){
        data <- data |>
          dplyr::mutate(source_id = as.character(source_id))
      }

      pcgr_ref_data[['biomarker']][[elem]] <- dplyr::bind_rows(
        pcgr_ref_data[['biomarker']][[elem]],
        data
      )
    }
  }

  ## Metadata
  pcgr_ref_data[['metadata']] <- data.frame()
  for(dtype in c('gene','gwas','hotspot','other',
                 'phenotype','biomarker','drug')){

    fname <- file.path(
      pcgr_db_assembly_dir, ".METADATA", "tsv",
      paste0(dtype,"_metadata.tsv")
    )
    check_file_exists(fname)
    metadata_dtype <- as.data.frame(
      readr::read_tsv(fname, show_col_types = F)) |>
      dplyr::mutate(datatype = dtype)

    pcgr_ref_data[['metadata']] <-
      pcgr_ref_data[['metadata']] |>
      dplyr::bind_rows(metadata_dtype)
  }

  return(pcgr_ref_data)


}
