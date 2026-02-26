# =============================================================================
# 01_dummy_data.R
# Purpose: Load and clean CDISC SDTM example domains from {pharmaversesdtm}.
#
# Pharmaverse packages used:
#   - pharmaversesdtm  Ready-made example SDTM datasets (dm, ae, ex, vs, lb)
#   - dplyr            Data manipulation
#
# About the source data:
#   The dm dataset in pharmaversesdtm is sourced from the CDISC pilot project,
#   a three-arm Alzheimer's disease study. The actual treatment arms are:
#     ARMCD "Xan_Hi"  → ARM "Xanomeline High Dose"
#     ARMCD "Xan_Lo"  → ARM "Xanomeline Low Dose"
#     ARMCD "Plbo"    → ARM "Placebo"
#
#   For this pipeline we treat "Xanomeline High Dose" as the active arm and
#   "Placebo" as the control, keeping the two-arm comparison structure used
#   throughout scripts 02–05.
#
# Exports: load_sdtm_data()
#   Returns a named list: list(dm, ae, ex, vs, lb)
# =============================================================================

library(pharmaversesdtm)
library(dplyr)
library(lubridate)

# -----------------------------------------------------------------------------
#' Load and clean SDTM domains from pharmaversesdtm
#'
#' Reads the built-in CDISC SDTM example datasets, applies minimal cleaning
#' (column selection, arm filtering), and returns them as a named list for use
#' by downstream ADaM derivation functions.
#'
#' The dm dataset originates from the CDISC pilot project (Alzheimer's study).
#' Actual ARMCD values are "Xan_Hi", "Xan_Lo", and "Plbo". This function
#' retains only "Xan_Hi" (Xanomeline High Dose) and "Plbo" (Placebo) to
#' maintain a two-arm comparison throughout the pipeline.
#'
#' No files are written to disk — data is passed directly to the next step
#' via the return value.
#'
#' @return A named list with elements: dm, ae, ex, vs, lb
# -----------------------------------------------------------------------------
load_sdtm_data <- function() {

  message("\n── Step 1: Loading SDTM domains ──────────────────────────────")

  # ── DM: Demographics ──────────────────────────────────────────────────────
  # Actual ARMCD values in the pharmaversesdtm CDISC pilot dm dataset:
  #   "Xan_Hi"  = Xanomeline High Dose  (active arm)
  #   "Xan_Lo"  = Xanomeline Low Dose   (excluded — keeping two-arm structure)
  #   "Plbo"    = Placebo               (control arm)
  #
  # Key variables:
  #   USUBJID           – unique subject identifier
  #   ARM / ARMCD       – treatment arm description / code
  #   AGE, SEX, RACE, ETHNIC
  #   RFSTDTC/RFENDTC   – reference start/end dates (ISO 8601 character)
  #   DTHFL / DTHDTC    – death flag and date
  data("dm", package = "pharmaversesdtm", envir = environment())

  # Diagnostic: show all arms present before filtering so the filter is
  # transparent and any future package changes will be immediately visible.
  message("  DM — all arms in source data:")
  dm %>% count(ARMCD, ARM) %>% as.data.frame() %>%
    apply(1, function(r) message("    ARMCD='", r["ARMCD"], "'  ARM='", r["ARM"],
                                 "'  n=", r["n"]))

  dm_clean <- dm %>%
    filter(ARMCD %in% c("Xan_Hi", "Plbo")) %>%       # High Dose vs Placebo
    select(STUDYID, USUBJID, SUBJID, SITEID,
           ARM, ACTARM, ARMCD, ACTARMCD,
           AGE, AGEU, SEX, RACE, ETHNIC,
           RFSTDTC, RFENDTC, DTHFL, DTHDTC)

  message("  DM after filtering to Xan_Hi + Plbo: ", nrow(dm_clean), " subjects")
  if (nrow(dm_clean) == 0) stop("dm_clean has 0 rows — check ARMCD values above.")

  # ── AE: Adverse Events ────────────────────────────────────────────────────
  # Key variables:
  #   AEDECOD  – MedDRA preferred term
  #   AEBODSYS – MedDRA system organ class (SOC)
  #   AESEV    – severity (MILD / MODERATE / SEVERE)
  #   AESER    – serious AE flag (Y / N)
  #   AEREL    – investigator-assessed relationship to study drug
  #
  # Filter to only subjects retained in dm_clean so domains stay in sync.
  data("ae", package = "pharmaversesdtm", envir = environment())
  ae_clean <- ae %>%
    filter(USUBJID %in% dm_clean$USUBJID) %>%
    select(STUDYID, USUBJID, AESEQ,
           AEDECOD, AEBODSYS, AESEV, AESER,
           AEREL, AEOUT, AESTDTC, AEENDTC)

  message("  AE: ", nrow(ae_clean), " event records")

  # ── EX: Exposure ──────────────────────────────────────────────────────────
  # Key variables: EXTRT, EXDOSE, EXDOSU, EXSTDTC, EXENDTC
  data("ex", package = "pharmaversesdtm", envir = environment())
  ex_clean <- ex %>%
    filter(USUBJID %in% dm_clean$USUBJID) %>%
    select(STUDYID, USUBJID, EXSEQ,
           EXTRT, EXDOSE, EXDOSU, EXDOSFRM,
           EXSTDTC, EXENDTC)

  message("  EX: ", nrow(ex_clean), " exposure records")

  # ── VS: Vital Signs ───────────────────────────────────────────────────────
  # Filtered to: systolic/diastolic BP, pulse, temperature, weight
  data("vs", package = "pharmaversesdtm", envir = environment())
  vs_clean <- vs %>%
    filter(USUBJID %in% dm_clean$USUBJID,
           VSTESTCD %in% c("SYSBP", "DIABP", "PULSE", "TEMP", "WEIGHT")) %>%
    select(STUDYID, USUBJID, VSSEQ, VSTESTCD, VSTEST,
           VSORRES, VSORRESU, VSSTRESN, VSSTRESU,
           VISITNUM, VISIT, VSDTC)

  message("  VS: ", nrow(vs_clean), " vital sign records")

  # ── LB: Laboratory Data ───────────────────────────────────────────────────
  # Filtered to six key analytes: ALT, AST, BILI, CREAT, HGB, WBC
  # LBNRLO / LBNRHI provide the normal range — used in shift table derivations
  data("lb", package = "pharmaversesdtm", envir = environment())
  lb_clean <- lb %>%
    filter(USUBJID %in% dm_clean$USUBJID,
           LBTESTCD %in% c("ALT", "AST", "BILI", "CREAT", "HGB", "WBC")) %>%
    select(STUDYID, USUBJID, LBSEQ, LBTESTCD, LBTEST,
           LBORRES, LBORRESU, LBSTRESN, LBSTRESU,
           LBSTNRLO, LBSTNRHI,      # standard-units normal range (not LBNRLO/LBNRHI)
           LBNRIND,                  # reference range indicator already in source data
           VISITNUM, VISIT, LBDTC)

  message("  LB: ", nrow(lb_clean), " lab records | analytes: ",
          paste(sort(unique(lb_clean$LBTESTCD)), collapse = ", "))
  message("── Step 1 complete ✔\n")

  # Return all domains as a named list — no files written to disk
  list(
    dm = dm_clean,
    ae = ae_clean,
    ex = ex_clean,
    vs = vs_clean,
    lb = lb_clean
  )
}
