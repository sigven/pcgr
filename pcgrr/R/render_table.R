#' Build column definitions for oncogenicity table with
#' category-aware styling
#'
#' This function generates a list of column definitions for a
#' reactable table that displays oncogenicity annotations.
#' It applies specific styling to the SYMBOL and SAMPLE_ALTERATION
#' columns based on cancer association rank and ClinVar classification,
#' respectively. The ONCOGENICITY column is styled with background colors
#' corresponding to oncogenicity categories. The function also ensures
#' that only relevant columns are shown, while others are hidden, and that
#' styling lookup columns are not displayed in the table.
#'
#' @param data The data frame containing the oncogenicity annotations
#' and styling lookup columns.
#' @param primary_cols Character vector of column names that should be
#' considered primary and kept visible in the table.
#' @param symbol_rank_col Name of the column in `data` that contains the
#' cancer association rank used for styling the SYMBOL column
#' (default: "TISSUE_ASSOC_RANK").
#'
#' @return A list of column definitions that can be passed to
#' the `colDef`
#'
#' @export
#'
build_oncogenicity_col_defs <- function(
    data,
    primary_cols,
    symbol_rank_col) {

  col_defs <- list(
    SYMBOL = reactable::colDef(
      name = "Gene",
      style = render_symbol_assoc_style(
        data, rank_col = symbol_rank_col
      ),
      minWidth = 90
    ),
    ALTERATION = reactable::colDef(
      name = "Alteration",
      style = render_alteration_clinvar_style(data),
      minWidth = 140
    ),
    CONSEQUENCE = reactable::colDef(
      name = "Consequence",
      minWidth = 130
    ),
    ONCOGENICITY = reactable::colDef(
      name = "Oncogenicity",
      cell = render_oncogenicity_cell(data),
      align = "center",
      minWidth = 140
    ),
    VAF_TUMOR = reactable::colDef(
      name = "Allelic Fraction",
      cell = render_bar_cell(fill_color = "#5D5165"),
      minWidth = 120,
      maxWidth = 170
    ),
    DP_TUMOR = reactable::colDef(
      name = "Depth",
      maxWidth = 70
    ),
    SEARCH_INDEX = reactable::colDef(
      show = FALSE,
      searchable = TRUE
    )
  )

  actual_cols <- colnames(data)

  ## Keep only col_defs for columns that exist in data
  col_defs <- col_defs[names(col_defs) %in% actual_cols]

  ## Hide columns used only for styling lookups
  styling_cols <- c(
    "CLINVAR_CLASSIFICATION",
    "TISSUE_ASSOC_RANK",
    "GLOBAL_ASSOC_RANK",
    "CODING_STATUS"
  )
  for (col in styling_cols) {
    if (col %in% actual_cols) {
      col_defs[[col]] <- reactable::colDef(show = FALSE)
    }
  }

  ## Hide any remaining columns not in primary_cols
  for (col in actual_cols) {
    if (!col %in% names(col_defs) && !col %in% primary_cols) {
      col_defs[[col]] <- reactable::colDef(show = FALSE)
    }
  }

  col_defs
}

#' Build column definitions for the germline findings table with
#' category-aware styling
#'
#' Generates a list of column definitions for a reactable table that
#' displays CPSR-classified germline variants. The ALTERATION column is
#' styled by the CPSR/ClinVar-derived CLINICAL_SIGNIFICANCE of the variant,
#' CLINICAL_SIGNIFICANCE itself is rendered as a colored badge, GENOTYPE is
#' flagged when backed by a low-depth control sample, and CPG_MOI (mode of
#' inheritance of the cancer predisposition gene) is shown with a tooltip
#' spelling out the abbreviated code(s). Other columns are hidden from the
#' main row (and thus fall through to the expandable row details).
#'
#' @param data The data frame containing the germline variant annotations
#' and styling lookup columns.
#' @param primary_cols Character vector of column names that should be
#' considered primary and kept visible in the table.
#' @param has_control_depth Whether DP_CONTROL reflects a real
#' (non-sentinel) sequencing depth for the variant set, used to enable the
#' low-depth GENOTYPE styling (default: FALSE).
#'
#' @return A list of column definitions that can be passed to reactable's
#' `columns` argument
#'
#' @export
#'
build_germline_col_defs <- function(
    data,
    primary_cols,
    has_control_depth = FALSE) {

  col_defs <- list(
    SYMBOL = reactable::colDef(
      name = "Gene",
      minWidth = 90
    ),
    ALTERATION = reactable::colDef(
      name = "Alteration",
      style = render_alteration_classification_style(data),
      minWidth = 140
    ),
    CONSEQUENCE = reactable::colDef(
      name = "Consequence",
      minWidth = 130
    ),
    GENOTYPE = reactable::colDef(
      name = "Genotype",
      style = if (has_control_depth) render_genotype_style(data) else NULL,
      align = "center",
      minWidth = 100,
      maxWidth = 120
    ),
    CLINICAL_SIGNIFICANCE = reactable::colDef(
      name = "Classification",
      cell = render_classification_cell(),
      align = "center",
      minWidth = 140
    ),
    ACMG_CODE = reactable::colDef(
      name = "ACMG criteria",
      cell = render_acmg_pills_cell(),
      minWidth = 170
    ),
    CPG_MOI = reactable::colDef(
      name = "Inheritance",
      cell = render_moi_cell(),
      align = "center",
      maxWidth = 110
    )
  )

  actual_cols <- colnames(data)

  ## Keep only col_defs for columns that exist in data
  col_defs <- col_defs[names(col_defs) %in% actual_cols]

  ## Hide any remaining columns not in primary_cols
  for (col in actual_cols) {
    if (!col %in% names(col_defs) && !col %in% primary_cols) {
      col_defs[[col]] <- reactable::colDef(show = FALSE)
    }
  }

  col_defs
}


#' Render one or more source logos from a pipe-separated string
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @return An HTML div element containing the corresponding logos,
#' or an empty string if no valid sources.
#'
#' @export
#'
render_source_logos <- function() {
  function(value) {
    if (is.na(value) || value == "") return("")

    sources <- strsplit(value, "\\|")[[1]]
    icons <- lapply(sources, function(src) {
      cls <- paste0("rt-logo-icon rt-logo-icon--", src)
      src_proper <- switch(
        src,
        civic = "CIViC",
        oncokb = "OncoKB",
        cgi = "Cancer Genome Interpreter",
        src)
      htmltools::span(class = cls, title = src_proper)
    })

    htmltools::div(class = "rt-logos-cell", icons)
  }
}

#' Render actionability tier as a styled badge with
#' tier-specific colors

#'
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param tier_letter Letter prefix for the tier (e.g. "T" for sensitivity, "R" for resistance)
#' @param tier_colors Named list with `t1` and `t2` hex colors for tiers 1 and 2
#' @return A function that can be used as a reactable cell renderer for the ACTIONABILITY_TIER column,
#'
#' @export
#'
render_tier_cell <- function(
    tier_letter = "T",
    tier_colors = NULL) {
  function(value) {
    if (is.na(value)) return("")
    bg <- if (value == 1) {
      tier_colors$t1
    } else if (value == 2) {
      tier_colors$t2
    } else if (value == 3) {
      tier_colors$t3
    } else {
      "transparent"
    }
    htmltools::div(
      class = "rt-badge",
      style = paste0("background-color:", bg, ";"),
      if (!is.na(value)) paste0(tier_letter, value) else ""
    )
  }
}



#' Render a value between 0 and 1 as a horizontal bar
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param fill_color Color of the bar
#' @param bg_color Background track color
#' @param show_value Whether to show the numeric value alongside
#'
#' @return A function that can be used as a reactable cell renderer
#' for numeric columns, which renders the value as a horizontal bar
#' @export
#'
render_bar_cell <- function(
    fill_color = "#5D5165",
    bg_color = "#e8e8e8",
    show_value = TRUE) {
  function(value) {
    if (is.na(value)) return("")
    width_pct <- max(0, min(100, round(value * 100, 1)))

    htmltools::div(
      class = "rt-bar-cell",
      htmltools::div(
        class = "rt-bar-track",
        style = paste0("background:", bg_color, ";"),
        htmltools::div(
          class = "rt-bar-fill",
          style = paste0("width:", width_pct, "%; background:", fill_color, ";")
        )
      ),
      if (show_value) {
        htmltools::span(
          class = "rt-bar-label",
          format(round(value, 2), nsmall = 2)
        )
      }
    )
  }
}


#' Style oncogenicity cell background based on value
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param data The main data frame (passed via closure)
#'
#' @return An HTML span element with background color corresponding
#' to the oncogenicity category, using colors defined in
#' color_palette$oncogenicity.
#'
#' @export
#'
render_oncogenicity_cell <- function(data, color_palette = pcgrr::color_palette) {
  function(value, index) {
    if (is.na(value)) return(value)

    levels <- color_palette$oncogenicity$levels
    colors <- color_palette$oncogenicity$values
    bg <- colors[match(value, levels)]
    if (is.na(bg)) bg <- "transparent"

    badge <- htmltools::span(
      class = "rt-badge rt-badge--sm",
      style = paste0("background-color:", bg, ";"),
      value
    )

    hotspot <- data$MUTATION_HOTSPOT_CANCERTYPE[index]
    show_hotspot <- isTRUE(
      length(hotspot) == 1L &&
        !is.na(hotspot) &&
        nchar(trimws(hotspot)) > 0
    )

    if (show_hotspot) {
      htmltools::div(
        class = "rt-hotspot-wrap",
        badge,
        htmltools::span(
          class = "rt-hotspot-icon rt-hotspot-icon--cancerhotspots",
          title = hotspot
        )
      )
    } else {
      badge
    }
  }
}

#' Style alteration cell background based on ClinVar classification
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param data The main data frame (passed via closure)
#'
#' @return A function that can be used as a reactable cell renderer
#' for the SAMPLE_ALTERATION column, which checks the CLINVAR_CLASSIFICATION
#' column in the same row to determine background color and text styling.
#'
#' @export
#'
render_alteration_clinvar_style <- function(data, color_palette = pcgrr::color_palette) {
  function(value, index) {
    clinvar <- data$CLINVAR_CLASSIFICATION[index]
    levels <- color_palette$pathogenicity_onc$levels
    colors <- color_palette$pathogenicity_onc$values
    bg <- if (!is.na(clinvar)) {
      matched <- colors[match(clinvar, levels)]
      if (!is.na(matched)) matched else "transparent"
    } else {
      "transparent"
    }
    txt_color <- if (bg != "transparent") "white" else "inherit"
    list(
      background = bg,
      color = txt_color,
      fontWeight = if (bg != "transparent") "bold" else "normal",
      padding = "4px 8px"
    )
  }
}

#' Style germline alteration cell background based on the
#' CPSR/ClinVar-derived clinical significance of the variant
#' (Generated by Claude Sonnet 5)
#'
#' @param data The main data frame (passed via closure)
#' @param class_col Name of the column holding the classification
#'   used to determine background color (default: "CLINICAL_SIGNIFICANCE")
#' @param color_palette Named list of color palettes
#'   (default: pcgrr::color_palette), keyed off the "pathogenicity" set
#'   (Pathogenic/Likely Pathogenic/VUS/Likely Benign/Benign)
#'
#' @return A function that can be used as a reactable cell renderer
#' for the ALTERATION column of the germline findings table
#'
#' @export
#'
render_alteration_classification_style <- function(
    data,
    class_col = "CLINICAL_SIGNIFICANCE",
    color_palette = pcgrr::color_palette) {
  levels <- color_palette$pathogenicity$levels
  colors <- color_palette$pathogenicity$values
  function(value, index) {
    cls <- data[[class_col]][index]
    bg <- if (!is.na(cls)) {
      matched <- colors[match(cls, levels)]
      if (!is.na(matched)) matched else "transparent"
    } else {
      "transparent"
    }
    list(
      background = bg,
      color = if (bg != "transparent") "white" else "inherit",
      fontWeight = if (bg != "transparent") "bold" else "normal",
      padding = "4px 8px"
    )
  }
}

#' Render germline clinical significance as a colored badge
#' (Generated by Claude Sonnet 5)
#'
#' @param color_palette Named list of color palettes
#'   (default: pcgrr::color_palette), keyed off the "pathogenicity" set
#'
#' @return A function that can be used as a reactable cell renderer
#' for the CLINICAL_SIGNIFICANCE column of the germline findings table
#'
#' @export
#'
render_classification_cell <- function(color_palette = pcgrr::color_palette) {
  function(value) {
    if (is.na(value)) return(value)
    levels <- color_palette$pathogenicity$levels
    colors <- color_palette$pathogenicity$values
    bg <- colors[match(value, levels)]
    if (is.na(bg)) bg <- "transparent"
    htmltools::span(
      class = "rt-badge rt-badge--sm",
      style = paste0("background-color:", bg, ";"),
      value
    )
  }
}

#' Flag germline GENOTYPE cell when backed by a low-depth control sample
#' (Generated by Claude Sonnet 5)
#'
#' Mirrors the DP_CONTROL-based flagging previously implemented via
#' \code{DT::formatStyle}: variants with a real (non-sentinel), low control
#' sequencing depth are marked with a colored left border and bold text, so
#' that their genotype call is interpreted with appropriate caution. DP_CONTROL
#' is set to -1 (sentinel) when CPSR is run without a matched control sample,
#' in which case no styling is applied.
#'
#' @param data The main data frame (passed via closure)
#' @param low_depth_threshold Sequencing depth (control sample) below
#'   which the GENOTYPE cell is flagged (default: 10)
#'
#' @return A function that can be used as a reactable cell style renderer
#' for the GENOTYPE column of the germline findings table
#'
#' @export
#'
render_genotype_style <- function(data, low_depth_threshold = 10) {
  function(value, index) {
    if (!"DP_CONTROL" %in% colnames(data)) {
      return(list())
    }
    dp <- data$DP_CONTROL[index]
    if (is.na(dp) || dp < 0 || dp >= low_depth_threshold) {
      return(list())
    }
    list(
      `border-left` = "4px solid #E69F00",
      fontWeight = "bold"
    )
  }
}

#' Render CPG_MOI (cancer predisposition gene - mode of inheritance) with
#' a tooltip spelling out the abbreviated code(s)
#' (Generated by Claude Sonnet 5)
#'
#' CPG_MOI is a gene-level annotation, aggregated across the disease
#' phenotype(s) linked to that gene in the CPSR gene panel. For most genes
#' this is a single mode (e.g. "AD" for TP53, "AR" for MUTYH), but for a few
#' genes (e.g. the mismatch repair genes, BRCA1) it combines two modes
#' (e.g. "AD/AR"), reflecting that monoallelic vs. biallelic pathogenic
#' variants in that gene cause distinct conditions with different
#' inheritance patterns (e.g. Lynch syndrome vs. CMMRD for PMS2/MSH6). The
#' cell keeps the short code(s) for a narrow column, with the full name(s)
#' available as a hover tooltip.
#'
#' @return A function that can be used as a reactable cell renderer for
#' the CPG_MOI column of the germline findings table
#'
#' @export
#'
render_moi_cell <- function() {
  moi_labels <- c(
    "AD"  = "Autosomal dominant",
    "AR"  = "Autosomal recessive",
    "XL"  = "X-linked",
    "XLD" = "X-linked dominant",
    "XLR" = "X-linked recessive",
    "MT"  = "Mitochondrial"
  )
  function(value) {
    if (is.na(value) || value == "." || value == "") return("")
    codes <- strsplit(value, "/")[[1]]
    full <- vapply(codes, function(x) {
      if (x %in% names(moi_labels)) moi_labels[[x]] else x
    }, character(1))
    htmltools::span(title = paste(full, collapse = " / "), value)
  }
}

## Condensed descriptions for CPSR/ACMG evidence codes, keyed by the
## display-form code used in the ACMG_CODE column (e.g. "PM2_supporting",
## "PVS1_strong") - a self-contained copy of the code/description pairs
## from cpsr::acmg$evidence_codes (pcgrr cannot import the cpsr package
## here, since cpsr itself depends on pcgrr). Used to power the hover
## tooltip on each ACMG criteria pill; keep in sync if the ACMG scoring
## scheme in cpsr changes.
acmg_evidence_descriptions <- c(
  PM1 = "Variant located in a mutational hotspot or critical functional domain without benign variation.",
  PM1_supporting = "Supporting evidence that a variant lies in a known functional hotspot or critical domain (supporting strength).",
  PM2_supporting = "Supporting evidence that the variant is absent or extremely rare in population databases (supporting strength).",
  BA1 = "Stand-alone evidence that the variant's allele frequency is too high for a pathogenic classification.",
  BP1 = "Supporting evidence that a missense variant occurs in a gene where truncating variants are predominantly known to cause disease.",
  BP4 = "Supporting evidence that multiple computational tools predict a benign effect on the gene or gene product.",
  BP7 = "Supporting evidence that a silent (synonymous) variant has no predicted impact on splicing or gene function.",
  BS1 = "Strong evidence that the variant's allele frequency is greater than expected for a disorder.",
  BS1_supporting = "Supporting evidence that the variant's frequency is slightly higher than expected for a pathogenic variant.",
  PVS1 = "Very strong evidence that a null (loss-of-function) variant occurs in a gene where loss of function is a known disease mechanism.",
  PVS1_strong = "Strong evidence for a predicted loss-of-function variant (reduced strength from PVS1).",
  PVS1_moderate = "Moderate evidence for a predicted loss-of-function variant (further reduced strength from PVS1).",
  PS1 = "Strong evidence that the variant causes the same amino acid change as a previously established pathogenic variant but via a different nucleotide change.",
  PP3 = "Supporting evidence that multiple computational tools predict a deleterious effect on the gene or gene product.",
  PM5 = "Evidence that the variant causes a novel amino acid change at a residue where another pathogenic missense change has been seen.",
  PM4 = "Evidence that the variant results in protein length changes due to in-frame deletions/insertions in a non-repeat region or stop-loss in functional protein domains.",
  PM4_supporting = "As PM4, but for single amino acid changes (supporting strength).",
  PP2 = "Supporting evidence that a missense variant occurs in a gene with low benign missense variation and where missense variants are a common disease mechanism."
)

#' Render ACMG criteria as individual colored pills
#' (Generated by Claude Sonnet 5)
#'
#' Splits a pipe-separated ACMG evidence code string (e.g.
#' "PM2_supporting|PVS1|PP3") into one pill per code. Pathogenic evidence
#' codes (PVS/PS/PM/PP*) are colored to match the "Pathogenic" tone of
#' \code{color_palette$pathogenicity}, benign evidence codes (BA/BS/BP*) are
#' colored to match its "Benign" tone, and any other/unrecognized code
#' falls back to a neutral gray pill. Each pill carries a hover tooltip
#' with a condensed description of the evidence code, sourced from
#' \code{acmg_evidence_descriptions}.
#'
#' @return A function that can be used as a reactable cell renderer for
#' the ACMG_CODE column of the germline findings table
#'
#' @export
#'
render_acmg_pills_cell <- function() {
  function(value) {
    if (is.na(value) || value == "." || value == "") return("")
    codes <- strsplit(value, "\\|")[[1]]
    pills <- lapply(codes, function(code) {
      variant_class <- if (startsWith(code, "P")) {
        "path"
      } else if (startsWith(code, "B")) {
        "benign"
      } else {
        "other"
      }
      desc <- unname(acmg_evidence_descriptions[code])
      htmltools::span(
        class = paste0("rt-acmg-pill rt-acmg-pill--", variant_class),
        title = if (!is.na(desc)) desc else NULL,
        code
      )
    })
    htmltools::div(class = "rt-acmg-cell", pills)
  }
}

#' Style alteration cell background based on cancer association rank
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param data The main data frame (passed via closure)
#' @param rank_col Name of the column containing the cancer association rank
#'   value used to determine cell background color (default: "TISSUE_ASSOC_RANK")
#'
#' @return A function that can be used as a reactable cell renderer
#' for the SAMPLE_ALTERATION column, which checks the TISSUE_ASSOC_RANK
#' column in the same row to determine background color based on predefined
#' breaks and colors in color_palette$cancer_assoc.
#'
#' @export
#'
render_symbol_assoc_style <- function(data, rank_col = "TISSUE_ASSOC_RANK", color_palette = pcgrr::color_palette) {
  breaks <- color_palette$cancer_assoc$breaks
  colors <- color_palette$cancer_assoc$values
  function(value, index) {
    rank_val <- data[[rank_col]][index]
    if (is.na(rank_val)) {
      return(list())
    }
    bin <- findInterval(rank_val, breaks) + 1
    bg <- colors[min(bin, length(colors))]
    list(
      background = bg,
      color = "white",
      padding = "4px 8px"
    )
  }
}



#' Render sample alteration with confidence icon
#' Uses filled circle + bold for high confidence,
#' hollow circle + muted for medium confidence
#'
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param data The main data frame (passed via closure)
#'
#' @return A function that can be used as a reactable
#' cell renderer for the
#' SAMPLE_ALTERATION column, which checks the BM_TOP_MAPPING_CONFIDENCE
#' column in the same row to determine styling.
#'
#' High confidence: filled circle (e.g. #37303A) + bold text
#' Medium confidence: hollow circle (border #999999) + normal weight + muted text color
#'
#' @export
#'
render_alteration_cell <- function(data) {
  function(value, index) {
    confidence <- data$BM_TOP_MAPPING_CONFIDENCE[index]

    is_high <- isTRUE(!is.na(confidence) && confidence == "high")
    htmltools::div(
      class = "rt-alteration-cell",
      htmltools::span(class = if (is_high) "rt-dot rt-dot--high" else "rt-dot rt-dot--medium"),
      htmltools::span(class = if (is_high) "rt-alteration-text rt-alteration-text--high"
                               else        "rt-alteration-text rt-alteration-text--medium",
                      value)
    )
  }
}

#' Render evidence level badge
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @return An HTML div element styled as a badge with background
#' and text color
#'
#' @export
#'
render_evidence_level_cell <- function(color_palette){
  function(value) {
    if (is.na(value) ||
        !value %in% names(color_palette$elevel)) {
      return(value)
    }
    colors <- color_palette$elevel[[value]]
    htmltools::span(
      class = "rt-badge",
      style = paste0("background-color:", colors$bg, "; color:", colors$fg, ";"),
      value
    )
  }
}


#' Render evidence description with truncation and tooltip
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @return An HTML div element that shows the first 200
#' characters of the description,
#' with a tooltip that displays the full text on hover
#' if it exceeds 200 characters.
#' @export
#'
render_evidence_desc_cell <- function(){

  function(value) {
    if (is.na(value)) return("")
    if (nchar(value) > 200) {
      htmltools::div(class = "rt-evidence-desc", title = value,
                     paste0(substr(value, 1, 200), "..."))
    } else {
      htmltools::div(class = "rt-evidence-desc", value)
    }
  }

}

#' Render therapy match cell with background color
#' based on actionability tier
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param data The main data frame (passed via closure)
#' @param tier_colors Named list with `t1` and `t2` hex
#' colors for tiers 1 and 2
#'
#' @export
#'
render_therapy_style <- function(data, tier_colors) {
  function(value, index) {
    tier <- data$ACTIONABILITY_TIER[index]
    bg <- if (!is.na(tier) && tier == 1) {
      tier_colors$t1
    } else if (!is.na(tier) && tier == 2) {
      tier_colors$t2
    } else if (!is.na(tier) && tier == 3) {
      tier_colors$t3
    }
    list(
      background = bg,
      color = "white",
      fontWeight = "bold",
      padding = "4px 8px"
    )
  }
}

#' Render prognostic outcome cell with background color
#' based on actionability tier
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param data The main data frame (passed via closure)
#' @param tier_colors Named list with `t1` and `t2` hex
#' colors for tiers 1 and 2
#'
#' @export
#'
render_prognostic_outcome <- function(data, tier_colors) {
  function(value, index) {
    tier <- data$ACTIONABILITY_TIER[index]
    bg <- if (!is.na(tier) && tier == 1) {
      tier_colors$t1
    } else {
      tier_colors$t2
    }
    list(
      background = bg,
      color = "white",
      fontWeight = "bold",
      padding = "4px 8px"
    )
  }
}

#' Render diagnosis cell with background color
#' based on actionability tier
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' @param data The main data frame (passed via closure)
#' @param tier_colors Named list with `t1` and `t2` hex
#' colors for tiers 1 and 2
#'
#' @export
#'
render_diagnosis <- function(data, tier_colors) {
  function(value, index) {
    tier <- data$ACTIONABILITY_TIER[index]
    bg <- if (!is.na(tier) && tier == 1) {
      tier_colors$t1
    } else {
      tier_colors$t2
    }
    list(
      background = bg,
      color = "white",
      fontWeight = "bold",
      padding = "4px 8px"
    )
  }
}

#' JS function for reactable row details with exclusions and card styling
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' Build JavaScript function for reactable row details
#' Generates a JS function that can be used in reactable's .details argument to
#' (optionally) render additional details for each row when expanded.
#' The function excludes specified primary columns, styling columns, and
#' reactable internals from the details view, and formats the remaining
#' key-value pairs as stacked cards with uppercase labels and values below.
#'
#' @param primary_cols Character vector of primary column names to
#' exclude from details
#' @param styling_cols Character vector of column names used for
#' styling (e.g. "CLINVAR_CLASSIFICATION", "TISSUE_ASSOC_RANK", etc.)
#' to exclude from details
#' @param font_size CSS font-size string applied to the detail card values
#'   (default: "0.94em")
#'
#' @return A JavaScript function as a string that can be passed to reactable's
#' .details argument to render row details with the specified exclusions and
#' formatting.
#'
#' @export
#'
build_rt_row_details <- function(
    primary_cols,
    styling_cols = c(
      "CLINVAR_CLASSIFICATION",
      "TISSUE_ASSOC_RANK",
      "GLOBAL_ASSOC_RANK",
      "CODING_STATUS"),
    font_size = "0.94em") {
  exclude_cols <-
    c(".details", primary_cols, styling_cols)
  exclude_js <-
    paste0("'", exclude_cols, "'", collapse = ", ")

  reactable::JS(paste0(
    "function(rowInfo) {",
    "  var exclude = [", exclude_js, "];",
    "  var items = Object.entries(rowInfo.row)",
    "    .filter(function(e) {",
    "      return exclude.indexOf(e[0]) === -1",
    "        && e[1] !== null",
    "        && e[1] !== ''",
    "        && String(e[1]).length > 0;",
    "    })",
    "    .map(function(e) {",
    "      return '<div style=\"margin-bottom:8px;\">'",
    "           + '<div style=\"font-weight:600;color:#888;",
    "font-size:",font_size,";text-transform:uppercase;",
    "letter-spacing:0.03em;margin-bottom:1px;\">'",
    "           + e[0] + '</div>'",
    "           + '<div style=\"color:#333;\">' + e[1] + '</div>'",
    "           + '</div>';",
    "    })",
    "    .join('');",
    "  if (!items) return null;",
    "  return '<div style=\"padding:10px 16px 12px 40px;background:#f9f9fb;",
    "display:flex;flex-wrap:wrap;gap:4px 24px;",
    "font-size:",font_size,";\">' + items + '</div>';",
    "}"
  ))
}

#' Predefined reactable theme for PCGR/CPSR tables
#' (Generated by Claude Opus 4.6 with some manual tweaks)
#'
#' Defines a consistent look for reactable tables in PCGR/CPSR reports, with
#' custom header styling (dark background, white bold text),
#' cell borders, and search input styling that matches the report's color palette.
#' This theme can be applied to all reactable tables in the report for a
#' cohesive appearance.
#'
#' @export
#'
rt_theme <- function(color_palette = pcgrr::color_palette){
  theme <-
    reactable::reactableTheme(
      style = list(fontFamily = "inherit"),
      headerStyle = list(
        background = color_palette$bg_dark,
        color = "white",
        fontWeight = "bold",
        fontSize = "0.97em",
        padding = "10px 14px",
        borderRight = "1px solid rgba(255,255,255,0.3)",
        display = "flex",
        alignItems = "center"
      ),
      cellStyle = list(
        borderRight = "1px solid #e8e8e8"
      ),
      searchInputStyle = list(
        borderColor = color_palette$bg_dark
      )
  )
  return(theme)
}



#' Build biomarker reactable with category-aware styling
#' Combines tier 1 and tier 2 records in one table.
#' Header uses tier 1 color; THERAPY_MATCH cell background
#' reflects the row's tier (1 or 2).
#' @param rctbl_recs List with $main and $nested data frames.
#'   $main must contain ACTIONABILITY_TIER with values 1 and 2.
#' @param variant_category One of "snv_indel", "cnv", "fusion"
#' @param clnsig One of "therapeutic_sensitivity" or
#' "therapeutic_resistance",
#'
#' @return A reactable object with the biomarker table
#' @export
#'
render_actble_bm_table <- function(
    rctbl_recs = NULL,
    variant_category = "snv_indel",
    clnsig = "therapeutic_sensitivity",
    color_palette = pcgrr::color_palette) {

  tier_colors <-
    list(
      t1 = color_palette$tier_sensitivity$values[1],
      t2 = color_palette$tier_sensitivity$values[2],
      t3 = color_palette$tier_sensitivity$values[3]
    )

  tier_letter <- "T"
  if(clnsig == "therapeutic_resistance"){
    tier_letter <- "R"
    tier_colors <- list(
      t1 = color_palette$tier_resistance$values[1],
      t2 = color_palette$tier_resistance$values[2],
      t3 = color_palette$tier_resistance$values[3]
    )
  }

  ## check that rctbl_recs is
  ## 1. non-null
  ## 2. is a list object that contains two elements
  ## 3. both elements are data frames
  ## 4. main data frame contains required columns
  ##.   (pending upon variant_category)
  if (is.null(rctbl_recs) ||
      !is.list(rctbl_recs) ||
      !all(c("main", "nested") %in% names(rctbl_recs)) ||
      !is.data.frame(rctbl_recs$main) ||
      !is.data.frame(rctbl_recs$nested)) {
    log4r_fatal(
      "rctbl_recs must be a list with 'main' and 'nested' data frames")
  }

  if(variant_category == "snv_indel"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "BM_TOP_MAPPING_CONFIDENCE",
        "VAF_TUMOR",
        "MUTATION_HOTSPOT",
        "ONCOGENICITY")
  } else if(variant_category == "cna"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "BM_TOP_MAPPING_CONFIDENCE",
        "CN_TOTAL")
  } else if(variant_category == "fusion"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "BM_TOP_MAPPING_CONFIDENCE")
  } else {
    log4r_fatal(
      paste0("Invalid variant_category: ", variant_category,
             ". Must be one of 'snv_indel', 'cna', or 'fusion'."))
  }

  if (NROW(rctbl_recs$main) == 0) {
    return(invisible(NULL))
  }

  assertable::assert_colnames(
    rctbl_recs$main,
    required_cols,
    only_colnames = FALSE,
    quiet = TRUE
  )

  assertable::assert_colnames(
    rctbl_recs$nested,
    c("VAR_ID",
      "ENTREZGENE",
      "ACTIONABILITY_TIER",
      "BM_MOLECULAR_PROFILE",
      "BM_REFERENCE",
      "BM_SOURCE_DB",
      "BM_CANCER_TYPE",
      "BM_THERAPEUTIC_CONTEXT",
      "BM_EVIDENCE_LEVEL",
      "BM_EVIDENCE_DESCRIPTION"),
    only_colnames = FALSE,
    quiet = TRUE
  )


  hdr_style <- list(
    background = tier_colors$t1,
    color = "white",
    fontFamily = "inherit",
    fontWeight = "bold",
    fontSize = "1.1em",
    #padding = "24px 18px",
    borderRight = "1px solid rgba(255,255,255,0.3)"
  )

  if(variant_category == "fusion"){
    main_cols <-  list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 100
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
      ),
      THERAPY_MATCH = reactable::colDef(
        name = switch(
          clnsig,
          therapeutic_sensitivity =
            "Therapy Match",
          therapeutic_resistance =
            "Resistance Match",
          "Therapy Match"
        ),
        html = TRUE,
        style = render_therapy_style(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER = reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE = reactable::colDef(show = FALSE)
    )
  }

  if(variant_category == "snv_indel"){
    main_cols = list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 140,       # icon + monospace text like "BCR::ABL1 fusion"
        maxWidth = 200
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
      ),
      ONCOGENICITY = reactable::colDef(
        name = "Oncogenicity",
        html = TRUE,
        align = "center",
        minWidth = 100
      ),
      VAF_TUMOR = reactable::colDef(
        name = "Allelic Fraction",
        cell = render_bar_cell(
          fill_color = "#5D5165"),
        minWidth = 100,
        maxWidth = 140
      ),
      THERAPY_MATCH = reactable::colDef(
        name = switch(
          clnsig,
          therapeutic_sensitivity =
            "Therapy Match",
          therapeutic_resistance =
            "Resistance Match",
          "Therapy Match"
        ),
        html = TRUE,
        minWidth = 200,
        style = render_therapy_style(
          rctbl_recs$main,
          tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      MUTATION_HOTSPOT = reactable::colDef(show = FALSE),
      ONCOGENICITY_CODE = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER = reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE = reactable::colDef(show = FALSE)
    )
  }

  if(variant_category == "cna"){
    main_cols = list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 140
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
        #maxWidth = 105
      ),
      CN_TOTAL = reactable::colDef(
        name = "Total Copy Number",
        align = "center",
        minWidth = 100
      ),
      THERAPY_MATCH = reactable::colDef(
        name = switch(
          clnsig,
          therapeutic_sensitivity =
            "Therapy Match",
          therapeutic_resistance =
            "Resistance Match",
          "Therapy Match"
        ),
        html = TRUE,
        minWidth = 200,
        style = render_therapy_style(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER =
        reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE =
        reactable::colDef(show = FALSE)
    )
  }

  result_table <- reactable::reactable(
    rctbl_recs$main,
    columns = main_cols,
    details = function(index) {
      eg <- rctbl_recs$main$ENTREZGENE[index]
      vid <- rctbl_recs$main$VAR_ID[index]
      tier <- rctbl_recs$main$ACTIONABILITY_TIER[index]
      nested <- rctbl_recs$nested[
        rctbl_recs$nested$ENTREZGENE == eg &
          rctbl_recs$nested$ACTIONABILITY_TIER == tier &
          rctbl_recs$nested$VAR_ID == vid,
      ]
      if (nrow(nested) == 0) return(NULL)
      display_cols <- setdiff(
        names(nested),
        c("VAR_ID",
          "ENTREZGENE",
          if (clnsig != "therapeutic_sensitivity") "ACTIONABILITY_TIER")
      )

      htmltools::div(
        style = "padding: 4px 40px 12px 40px; background: #f9f9f9;",
        reactable::reactable(
          nested[, display_cols],
          columns = c(
            list(
              BM_MOLECULAR_PROFILE = reactable::colDef(
                name = "Molecular Profile",
                html = TRUE,
                minWidth = 140
              ),
              BM_REFERENCE = reactable::colDef(
                name = "Reference",
                html = TRUE,
                minWidth = 180
              ),
              BM_SOURCE_DB = reactable::colDef(
                name = "Source",
                maxWidth = 100
              ),
              BM_CANCER_TYPE = reactable::colDef(
                name = "Cancer Type",
                minWidth = 130
              )
            ),
            if (clnsig == "therapeutic_sensitivity") list(
              ACTIONABILITY_TIER = reactable::colDef(
                name = "Tier",
                align = "center",
                cell = render_tier_cell(
                  tier_letter = tier_letter,
                  tier_colors = tier_colors
                ),
                minWidth = 80
              )
            ),
            list(
              BM_THERAPEUTIC_CONTEXT = reactable::colDef(
                name = "Therapy",
                minWidth = 120
              ),
              BM_EVIDENCE_LEVEL = reactable::colDef(
                name = "Evidence level",
                cell = render_evidence_level_cell(color_palette),
                align = "center",
                minWidth = 100
              ),
              BM_EVIDENCE_DESCRIPTION = reactable::colDef(
                name = "Evidence Description",
                cell = render_evidence_desc_cell(),
                minWidth = 350
              )
            )
          ),
          outlined = TRUE,
          compact = TRUE,
          highlight = TRUE,
          wrap = TRUE,
          defaultPageSize = 5,
          theme = reactable::reactableTheme(
            backgroundColor = "#f9f9f9"
          )
        )
      )
    },
    searchable = TRUE,
    striped = TRUE,
    highlight = TRUE,
    compact = TRUE,
    filterable = TRUE,
    defaultPageSize = 5,
    theme = reactable::reactableTheme(
      style = list(fontFamily = "inherit"),
      headerStyle = hdr_style,
      cellStyle = list(
        borderRight = "1px solid #e0e0e0",
        display = "flex",
        alignItems = "center"),
      searchInputStyle = list(
        borderColor = tier_colors$t1
      )
    )
  )


  return(result_table)
}


#' Build biomarker reactable with category-aware styling
#' Combines tier 1 and tier 2 records in one table.
#' Header uses tier 1 color; THERAPY_MATCH cell background
#' reflects the row's tier (1 or 2).
#' @param rctbl_recs List with $main and $nested data frames.
#'   $main must contain ACTIONABILITY_TIER with values 1 and 2.
#' @param variant_category One of "snv_indel", "cnv", "fusion"
#' @param clnsig One of "therapeutic_sensitivity" or
#' "therapeutic_resistance",
#'
#' @return A reactable object with the biomarker table
#' @export
#'
render_progn_bm_table <- function(
    rctbl_recs = NULL,
    variant_category = "snv_indel",
    clnsig = "prognostic_poor",
    color_palette = pcgrr::color_palette) {

  tier_colors <-
    list(
      t1 = color_palette$prognosis$values[2],
      t2 = color_palette$prognosis$values[2]
    )

  tier_letter <- "PP"
  if(clnsig == "prognostic_better"){
    tier_letter <- "PB"
    tier_colors <- list(
      t1 = color_palette$prognosis$values[1],
      t2 = color_palette$prognosis$values[1]
    )
  }

  ## check that rctbl_recs is
  ## 1. non-null
  ## 2. is a list object that contains two elements
  ## 3. both elements are data frames
  ## 4. main data frame contains required columns
  ##.   (pending upon variant_category)
  if (is.null(rctbl_recs) ||
      !is.list(rctbl_recs) ||
      !all(c("main", "nested") %in% names(rctbl_recs)) ||
      !is.data.frame(rctbl_recs$main) ||
      !is.data.frame(rctbl_recs$nested)) {
    log4r_fatal(
      "rctbl_recs must be a list with 'main' and 'nested' data frames")
  }

  if(variant_category == "snv_indel"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "BM_TOP_MAPPING_CONFIDENCE",
        "VAF_TUMOR",
        "MUTATION_HOTSPOT",
        "PROGNOSTIC_OUTCOME",
        "ONCOGENICITY")
  } else if(variant_category == "cna"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "PROGNOSTIC_OUTCOME",
        "BM_TOP_MAPPING_CONFIDENCE",
        "CN_TOTAL")
  } else if(variant_category == "fusion"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "PROGNOSTIC_OUTCOME",
        "BM_TOP_MAPPING_CONFIDENCE")
  } else {
    log4r_fatal(
      paste0("Invalid variant_category: ", variant_category,
             ". Must be one of 'snv_indel', 'cna', or 'fusion'."))
  }

  assertable::assert_colnames(
    rctbl_recs$main,
    required_cols,
    only_colnames = FALSE,
    quiet = TRUE
  )

  assertable::assert_colnames(
    rctbl_recs$nested,
    c("VAR_ID",
      "ENTREZGENE",
      "ACTIONABILITY_TIER",
      "BM_MOLECULAR_PROFILE",
      "BM_REFERENCE",
      "BM_SOURCE_DB",
      "BM_CANCER_TYPE",
      "BM_EVIDENCE_LEVEL",
      "BM_EVIDENCE_DESCRIPTION"),
    only_colnames = FALSE,
    quiet = TRUE
  )


  hdr_style <- list(
    background = tier_colors$t1,
    color = "white",
    fontFamily = "inherit",
    fontWeight = "bold",
    fontSize = "1.1em",
    #padding = "24px 18px",
    borderRight = "1px solid rgba(255,255,255,0.3)"
  )

  if(variant_category == "fusion"){
    main_cols <-  list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 100
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
      ),
      PROGNOSTIC_OUTCOME = reactable::colDef(
        name = "Prognostic Outcome",
        html = TRUE,
        style = render_prognostic_outcome(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER = reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE = reactable::colDef(show = FALSE)
    )
  }

  if(variant_category == "snv_indel"){
    main_cols = list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 140,   # icon + monospace text like "BCR::ABL1 fusion"
        maxWidth = 200
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
      ),
      ONCOGENICITY = reactable::colDef(
        name = "Oncogenicity",
        html = TRUE,
        align = "center",
        minWidth = 100
      ),
      VAF_TUMOR = reactable::colDef(
        name = "Allelic Fraction",
        cell = render_bar_cell(
          fill_color = "#5D5165"),
        minWidth = 100,
        maxWidth = 140
      ),
      PROGNOSTIC_OUTCOME = reactable::colDef(
        name = "Prognostic Outcome",
        html = TRUE,
        style = render_prognostic_outcome(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      MUTATION_HOTSPOT = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER = reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE = reactable::colDef(show = FALSE)
    )
  }

  if(variant_category == "cna"){
    main_cols = list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 140
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
        #maxWidth = 105
      ),
      CN_TOTAL = reactable::colDef(
        name = "Total Copy Number",
        align = "center",
        minWidth = 100
      ),
      PROGNOSTIC_OUTCOME = reactable::colDef(
        name = "Prognostic Outcome",
        html = TRUE,
        style = render_prognostic_outcome(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER =
        reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE =
        reactable::colDef(show = FALSE)
    )
  }

  result_table <- reactable::reactable(
    rctbl_recs$main,
    columns = main_cols,
    details = function(index) {
      eg <- rctbl_recs$main$ENTREZGENE[index]
      vid <- rctbl_recs$main$VAR_ID[index]
      tier <- rctbl_recs$main$ACTIONABILITY_TIER[index]
      nested <- rctbl_recs$nested[
        rctbl_recs$nested$ENTREZGENE == eg &
          rctbl_recs$nested$ACTIONABILITY_TIER == tier &
          rctbl_recs$nested$VAR_ID == vid,
      ]
      if (nrow(nested) == 0) return(NULL)
      display_cols <- setdiff(
        names(nested),
        c("VAR_ID",
          "ACTIONABILITY_TIER",
          "ENTREZGENE")
      )
      htmltools::div(
        style = "padding: 4px 40px 12px 40px; background: #f9f9f9;",
        reactable::reactable(
          nested[, display_cols],
          columns = list(
            BM_MOLECULAR_PROFILE = reactable::colDef(
              name = "Molecular Profile",
              html = TRUE,
              minWidth = 140
            ),
            BM_REFERENCE = reactable::colDef(
              name = "Reference",
              html = TRUE,
              minWidth = 180
            ),
            BM_SOURCE_DB = reactable::colDef(
              name = "Source",
              maxWidth = 100
            ),
            BM_CANCER_TYPE = reactable::colDef(
              name = "Cancer Type",
              minWidth = 130
            ),
            BM_EVIDENCE_LEVEL = reactable::colDef(
              name = "Evidence level",
              cell = render_evidence_level_cell(color_palette),
              align = "center",
              minWidth = 100
            ),
            BM_EVIDENCE_DESCRIPTION = reactable::colDef(
              name = "Evidence Description",
              cell = render_evidence_desc_cell(),
              minWidth = 350
            ),
            BM_CLINICAL_SIGNIFICANCE = reactable::colDef(
              name = "Outcome",
              align = "center",
              minWidth = 120)
          ),
          outlined = TRUE,
          compact = TRUE,
          highlight = TRUE,
          wrap = TRUE,
          defaultPageSize = 5,
          theme = reactable::reactableTheme(
            backgroundColor = "#f9f9f9"
          )
        )
      )
    },
    searchable = TRUE,
    striped = TRUE,
    highlight = TRUE,
    filterable = TRUE,
    compact = TRUE,
    defaultPageSize = 5,
    theme = reactable::reactableTheme(
      style = list(fontFamily = "inherit"),
      headerStyle = hdr_style,
      cellStyle = list(
        borderRight = "1px solid #e0e0e0",
        display = "flex",
        alignItems = "center"),
      searchInputStyle = list(
        borderColor = tier_colors$t1
      )
    )
  )


  return(result_table)
}


#' Build biomarker reactable with category-aware styling
#' Combines tier 1 and tier 2 records in one table.
#' Header uses tier 1 color; THERAPY_MATCH cell background
#' reflects the row's tier (1 or 2).
#' @param rctbl_recs List with $main and $nested data frames.
#'   $main must contain ACTIONABILITY_TIER with values 1 and 2.
#' @param variant_category One of "snv_indel", "cnv", "fusion"
#' @param clnsig One of "therapeutic_sensitivity" or
#' "therapeutic_resistance",
#'
#' @return A reactable object with the biomarker table
#' @export
#'
render_diagn_bm_table <- function(
    rctbl_recs = NULL,
    variant_category = "snv_indel",
    clnsig = "diagnostic_positive",
    color_palette = pcgrr::color_palette) {

  tier_colors <-
    list(
      t1 = color_palette$diagnosis,
      t2 = color_palette$diagnosis
    )

  tier_letter <- "D"
  # if(clnsig == "prognostic_better"){
  #   tier_letter <- "PB"
  #   tier_colors <- list(
  #     t1 = color_palette$prognosis$values[1],
  #     t2 = color_palette$prognosis$values[1]
  #   )
  # }

  ## check that rctbl_recs is
  ## 1. non-null
  ## 2. is a list object that contains two elements
  ## 3. both elements are data frames
  ## 4. main data frame contains required columns
  ##.   (pending upon variant_category)
  if (is.null(rctbl_recs) ||
      !is.list(rctbl_recs) ||
      !all(c("main", "nested") %in% names(rctbl_recs)) ||
      !is.data.frame(rctbl_recs$main) ||
      !is.data.frame(rctbl_recs$nested)) {
    log4r_fatal(
      "rctbl_recs must be a list with 'main' and 'nested' data frames")
  }

  if(variant_category == "snv_indel"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "BM_TOP_MAPPING_CONFIDENCE",
        "VAF_TUMOR",
        "MUTATION_HOTSPOT",
        "DIAGNOSTIC_EVIDENCE",
        "ONCOGENICITY")
  } else if(variant_category == "cna"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "DIAGNOSTIC_EVIDENCE",
        "BM_TOP_MAPPING_CONFIDENCE",
        "CN_TOTAL")
  } else if(variant_category == "fusion"){
    required_cols <-
      c("VAR_ID",
        "ENTREZGENE",
        "BM_SOURCES",
        "ACTIONABILITY_TIER",
        "DIAGNOSTIC_EVIDENCE",
        "BM_TOP_MAPPING_CONFIDENCE")
  } else {
    log4r_fatal(
      paste0("Invalid variant_category: ", variant_category,
             ". Must be one of 'snv_indel', 'cna', or 'fusion'."))
  }

  assertable::assert_colnames(
    rctbl_recs$main,
    required_cols,
    only_colnames = FALSE,
    quiet = TRUE
  )

  assertable::assert_colnames(
    rctbl_recs$nested,
    c("VAR_ID",
      "ENTREZGENE",
      "ACTIONABILITY_TIER",
      "BM_MOLECULAR_PROFILE",
      "BM_REFERENCE",
      "BM_SOURCE_DB",
      "BM_CANCER_TYPE",
      "BM_EVIDENCE_LEVEL",
      "BM_EVIDENCE_DESCRIPTION"),
    only_colnames = FALSE,
    quiet = TRUE
  )


  hdr_style <- list(
    background = tier_colors$t1,
    color = "white",
    fontFamily = "inherit",
    fontWeight = "bold",
    fontSize = "1.1em",
    #padding = "24px 18px",
    borderRight = "1px solid rgba(255,255,255,0.3)"
  )

  if(variant_category == "fusion"){
    main_cols <-  list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 100
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
      ),
      DIAGNOSTIC_EVIDENCE = reactable::colDef(
        name = "Diagnostic Evidence",
        html = TRUE,
        style = render_diagnosis(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER = reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE = reactable::colDef(show = FALSE)
    )
  }

  if(variant_category == "snv_indel"){
    main_cols = list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 140,   # icon + monospace text like "BCR::ABL1 fusion"
        maxWidth = 200
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
      ),
      ONCOGENICITY = reactable::colDef(
        name = "Oncogenicity",
        html = TRUE,
        minWidth = 100
      ),
      VAF_TUMOR = reactable::colDef(
        name = "Allelic Fraction",
        cell = render_bar_cell(
          fill_color = "#5D5165"),
        minWidth = 100,
        maxWidth = 140
      ),
      DIAGNOSTIC_EVIDENCE = reactable::colDef(
        name = "Diagnostic Evidence",
        html = TRUE,
        style = render_diagnosis(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      MUTATION_HOTSPOT = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER = reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE = reactable::colDef(show = FALSE)
    )
  }

  if(variant_category == "cna"){
    main_cols = list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = render_alteration_cell(
          rctbl_recs$main),
        minWidth = 140
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = render_source_logos(),
        align = "center",
        minWidth = 70
        #maxWidth = 105
      ),
      CN_TOTAL = reactable::colDef(
        name = "Total Copy Number",
        align = "center",
        minWidth = 100
      ),
      DIAGNOSTIC_EVIDENCE = reactable::colDef(
        name = "Diagnostic Evidence",
        html = TRUE,
        style = render_diagnosis(
          rctbl_recs$main, tier_colors)
      ),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      ACTIONABILITY_TIER =
        reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE =
        reactable::colDef(show = FALSE)
    )
  }

  result_table <- reactable::reactable(
    rctbl_recs$main,
    columns = main_cols,
    details = function(index) {
      eg <- rctbl_recs$main$ENTREZGENE[index]
      vid <- rctbl_recs$main$VAR_ID[index]
      tier <- rctbl_recs$main$ACTIONABILITY_TIER[index]
      nested <- rctbl_recs$nested[
        rctbl_recs$nested$ENTREZGENE == eg &
          rctbl_recs$nested$ACTIONABILITY_TIER == tier &
          rctbl_recs$nested$VAR_ID == vid,
      ]
      if (nrow(nested) == 0) return(NULL)
      display_cols <- setdiff(
        names(nested),
        c("VAR_ID",
          "ACTIONABILITY_TIER",
          "ENTREZGENE")
      )
      htmltools::div(
        style = "padding: 4px 40px 12px 40px; background: #f9f9f9;",
        reactable::reactable(
          nested[, display_cols],
          columns = list(
            BM_CLINICAL_SIGNIFICANCE = reactable::colDef(
              name = "Diagnostic Evidence",
              align = "center",
              minWidth = 140),
            BM_MOLECULAR_PROFILE = reactable::colDef(
              name = "Molecular Profile",
              html = TRUE,
              minWidth = 140
            ),
            BM_REFERENCE = reactable::colDef(
              name = "Reference",
              html = TRUE,
              minWidth = 180
            ),
            BM_SOURCE_DB = reactable::colDef(
              name = "Source",
              maxWidth = 100
            ),
            BM_CANCER_TYPE = reactable::colDef(
              name = "Cancer Type",
              minWidth = 130
            ),
            BM_EVIDENCE_LEVEL = reactable::colDef(
              name = "Evidence level",
              cell = render_evidence_level_cell(color_palette),
              align = "center",
              minWidth = 100
            ),
            BM_EVIDENCE_DESCRIPTION = reactable::colDef(
              name = "Evidence Description",
              cell = render_evidence_desc_cell(),
              minWidth = 350
            )
          ),
          outlined = TRUE,
          compact = TRUE,
          highlight = TRUE,
          wrap = TRUE,
          defaultPageSize = 5,
          theme = reactable::reactableTheme(
            backgroundColor = "#f9f9f9"
          )
        )
      )
    },
    searchable = TRUE,
    striped = TRUE,
    highlight = TRUE,
    filterable = TRUE,
    compact = TRUE,
    defaultPageSize = 5,
    theme = reactable::reactableTheme(
      style = list(fontFamily = "inherit"),
      headerStyle = hdr_style,
      cellStyle = list(
        borderRight = "1px solid #e0e0e0",
        display = "flex",
        alignItems = "center"),
      searchInputStyle = list(
        borderColor = tier_colors$t1
      )
    )
  )


  return(result_table)
}


#' Emit the shared "biomarker types and report scope" callout note
#'
#' @param variant_label Variant-type phrase for the Tier III sentence,
#'   e.g. "SNV/InDel variants", "copy number aberrations", "RNA fusions".
#' @export
callout_biomarker_scope <- function(variant_label = "variants") {
  cat(paste0(
    "::: {.callout-note collapse=\"false\"}\n",
    "## Note — biomarker types and report scope\n\n",
    "The ultimate _variant tier classification_ in PCGR is driven by **therapeutic sensitivity** biomarkers ",
    "(i.e. those related to treatment response). Evidence items where variants match ",
    "**resistance**, **prognostic**, or **diagnostic** biomarkers are also shown — ",
    "all restricted to biomarkers matching the tumor type of the query sample.\n\n",
    "The complete set of **Tier III** ", variant_label,
    " (uncertain/unknown clinical significance) is only available in the TSV and Excel workbook outputs.\n\n",
    ":::\n"
  ))
}
