context("Biomarker check")
pcgr_data <- readRDS(file='/Users/sigven/research/docker/pcgr/data/grch37/rds/pcgr_data.rds')

eitems_raw <- pcgr_data[["biomarkers"]]

test_that("Biomarker load returns correct origin/type", {
  expect_error(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "",
                                             origin = "Somatic"))
  expect_error(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                      alteration_type = "MUT",
                                      origin = ""))
  expect_error(pcgrr::load_all_eitems(eitems_raw = NULL,
                                      alteration_type = "MUT",
                                      origin = "Somatic"))
  expect_equal(unique(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "MUT",
                                             origin = "Somatic")$VARIANT_ORIGIN),
               c('Somatic','Somatic Mutation'))
  expect_equal(unique(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "MUT",
                                             origin = "Germline")$VARIANT_ORIGIN),
               "Germline")
  expect_equal(unique(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "CNA",
                                             origin = "Somatic")$VARIANT_ORIGIN),
               c("Somatic", "Somatic Mutation"))
  expect_equal(unique(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "MUT",
                                             origin = "Somatic")$ALTERATION_TYPE),
               "MUT")
  expect_equal(unique(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "CNA",
                                             origin = "Somatic")$ALTERATION_TYPE),
               "CNA")
  expect_equal(unique(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "CNA",
                                             origin = "Germline")$ALTERATION_TYPE),
               "CNA")
  expect_gte(ncol(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                              alteration_type = "MUT",
                              origin = "Somatic")), 19)
  expect_gte(nrow(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                         alteration_type = "MUT_LOF",
                                         origin = "Germline")), 1)
  expect_gte(ncol(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                         alteration_type = "CNA",
                                         origin = "Somatic")), 17)
  expect_equal(unique(pcgrr::load_all_eitems(eitems_raw = eitems_raw,
                                             alteration_type = "MUT",
                                             origin = "Somatic")$BIOMARKER_MAPPING),
               c("exact", "gene", "codon", "exon"))
})
