## Cohort definition for GSE113957 (SRP144355), fixed 2026-09-15 after the
## pre-submission review found that the published 142-donor cohort contains
## nine Hutchinson-Gilford progeria (HGPS) donors and 31 donors under 18, and
## that every donor aged 83 and over comes from one cell repository (Coriell AG).
##
## Source of truth: the SRA sample attributes carried by recount3
##   age;;1|cell id;;AG08498|disease;;Normal|ethnicity;;Asian|Sex;;male|source_name;;...
## Age is read from the attribute, not from the sample title, so paediatric
## samples given in months are handled explicitly.
##
## PRIMARY cohort  : normal donors aged 20 years and over        (adults)
## Covariates      : sequencing depth, cell repository (AG vs other), instrument,
##                   and sex; proliferation score where an adjustment is made
## SENSITIVITY     : all normal donors (children included); normal adults aged
##                   20 to 82 (drops the stratum in which repository and age are
##                   collinear); the published 142-donor cohort
##
## Sourced by the analysis scripts; running it directly writes the donor table.

gse113957_donors <- function(rc = "public_data_tierA/recount3", srp = "SRP144355") {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  md <- data.table::fread(cmd = paste("zcat", file.path(rc, sprintf("sra.sra.%s.MD.gz", srp))))
  att <- function(field) {
    pat <- paste0("(^|\\|)", field, ";;")
    vapply(strsplit(md$sample_attributes, "|", fixed = TRUE), function(kv) {
      hit <- kv[startsWith(tolower(kv), paste0(tolower(field), ";;"))]
      if (!length(hit)) return(NA_character_)
      sub("^[^;]*;;", "", hit[1]) }, character(1))
  }
  age_raw <- att("age")
  ## normal donors carry a whole number of years; the progeria samples are given
  ## as "8yr" or "2yr3mos" (one as "6ys11mos"), so years and months are read out
  yr  <- suppressWarnings(as.numeric(sub("^([0-9]+)(yr|ys|YR).*$", "\\1", age_raw)))
  mos <- suppressWarnings(as.numeric(sub("^.*[^0-9]([0-9]+)mos$", "\\1", age_raw)))
  age <- suppressWarnings(as.numeric(age_raw))
  age[is.na(age)] <- yr[is.na(age)] + ifelse(is.na(mos[is.na(age)]), 0, mos[is.na(age)]) / 12
  cell_id <- att("cell id")
  D <- data.frame(
    rail    = as.character(md$rail_id),
    srr     = md$external_id,
    gsm     = md$sample_name,
    title   = md$sample_title,
    cell_id = cell_id,
    ## Coriell prefixes: AG = NIA Aging Cell Repository, GM = general collection,
    ## HGADFN = Progeria Research Foundation cell bank, PRF = same donor source
    repo    = ifelse(grepl("^AG", cell_id), "AG", "other"),
    disease = att("disease"),
    age     = age,
    sex     = tolower(att("Sex")),
    ethnic  = att("ethnicity"),
    instr   = md$platform_model,
    site    = att("source_name"),
    stringsAsFactors = FALSE)
  D$normal <- D$disease == "Normal"
  D$adult  <- is.finite(D$age) & D$age >= 20
  D$primary <- D$normal & D$adult                    # PRIMARY cohort
  D
}

## keep(D, "primary" | "normal" | "adult2082" | "published")
gse113957_keep <- function(D, cohort = "primary") {
  switch(cohort,
    primary   = D$primary,
    normal    = D$normal,
    adult2082 = D$primary & D$age < 83,
    published = rep(TRUE, nrow(D)),
    stop("unknown cohort: ", cohort))
}

if (sys.nframe() == 0) {
  .libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
  D <- gse113957_donors()
  O <- "public_data_tierA/derived/cohort"; dir.create(O, showWarnings = FALSE, recursive = TRUE)
  write.table(D, file.path(O, "GSE113957_donors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("samples: %d | ages %.0f-%.0f | HGPS %d | under 20 (normal) %d\n",
      nrow(D), min(D$age, na.rm = TRUE), max(D$age, na.rm = TRUE),
      sum(D$disease == "HGPS"), sum(D$normal & !D$adult)))
  for (k in c("published", "normal", "primary", "adult2082")) {
    s <- D[gse113957_keep(D, k), ]
    cat(sprintf("  %-10s n = %3d  ages %2.0f-%2.0f  AG %3d / other %3d  HiSeq %2d\n",
        k, nrow(s), min(s$age), max(s$age), sum(s$repo == "AG"), sum(s$repo != "AG"),
        sum(s$instr != "NextSeq 500")))
  }
  cat("\nrepository by age band, normal donors:\n")
  N <- D[D$normal, ]
  print(table(cut(N$age, c(0, 20, 40, 60, 83, 100), right = FALSE), N$repo))
}
