# ============================================================
# export_for_dashboard_v2.R
# อัปเดตให้ตรงกับ F_Part3_Analysis.R (REVISED — พ.ค. 2569)
#
# สิ่งที่เปลี่ยนจาก v1:
#   - factor levels ใช้ underscore: "In_Catchment" / "Out_of_Catchment"
#   - continuity_flag: "0_NotContinuous" / "1_Continuous"
#   - ooc_subgroup3 แทน ooc_type
#   - cost model เพิ่ม continuity_flag เป็น covariate
#   - term name: "out_of_catchmentOut_of_Catchment"
#   - covariates ตรงกับ Part 3 ทุก step
#   - ไม่ Winsorize ซ้ำ (ทำใน Part 2 แล้ว — cost_w พร้อมใช้)
#
# วิธีใช้:
#   1. เปิด RStudio → setwd ให้ตรงกับโปรเจกต์
#   2. source("export_for_dashboard_v2.R")
#   3. นำ dashboard_data.json ไปวางใน GitHub repo แทนของเก่า
# ============================================================

setwd("C:/Users/ADMin/Desktop/HT Project")

cat("=========================================\n")
cat("  EXPORT FOR DASHBOARD v2\n")
cat("  (ตรงกับ F_Part3_Analysis.R — REVISED)\n")
cat("=========================================\n\n")

# ── Packages ──────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "sandwich", "lmtest", "lme4",
          "pROC", "ResourceSelection", "jsonlite")
new_pkg <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_pkg) > 0) install.packages(new_pkg)
invisible(lapply(pkgs, library, character.only = TRUE))

OUTPUT_JSON <- "dashboard_data.json"

# ============================================================
# 1. โหลดข้อมูล (จาก checkpoint เดียวกับ Part 3)
# ============================================================
cat("📂 โหลดข้อมูล...\n")
stopifnot(file.exists("checkpoint/DM_FINAL_ready_for_analysis.rds"))
dat_raw <- readRDS("checkpoint/DM_FINAL_ready_for_analysis.rds")
cat("   n รวม:", nrow(dat_raw), "| ตัวแปร:", ncol(dat_raw), "\n")

# ── เตรียมตัวแปร (เหมือน Step A ของ Part 3) ─────────────────
dat <- dat_raw %>%
  filter(
    out_of_catchment %in% c("In_Catchment", "Out_of_Catchment"),
    continuity_flag  %in% c("0_NotContinuous", "1_Continuous")
  ) %>%
  mutate(
    hcode = as.character(hcode),

    out_of_catchment = factor(as.character(out_of_catchment),
                              levels = c("In_Catchment", "Out_of_Catchment")),

    ooc_subgroup3 = factor(as.character(ooc_subgroup3),
                           levels = c("In_Catchment", "OOC_Bangkok", "OOC_Outside")),

    continuity_flag = factor(as.character(continuity_flag),
                             levels = c("0_NotContinuous", "1_Continuous")),

    continuity_flag01 = case_when(
      continuity_flag == "0_NotContinuous" ~ 0L,
      continuity_flag == "1_Continuous"    ~ 1L,
      TRUE ~ NA_integer_
    ),

    bp_controlled = case_when(
      bp_controlled %in% c(0, 1) ~ as.integer(bp_controlled),
      TRUE ~ NA_integer_
    ),

    sex_bin       = factor(as.character(sex_bin),
                           levels = c("0_Female", "1_Male")),
    edu_group     = factor(as.character(edu_group),
                           levels = c("1_Primary_or_less","2_Secondary",
                                      "3_Bachelor_plus","9_Not_Specified")),
    occ_group5    = factor(as.character(occ_group5),
                           levels = c("1_ไม่มีอาชีพ/งานบ้าน","2_รับจ้าง",
                                      "3_ค้าขาย/ธุรกิจ","4_ภาครัฐ/ราชการ",
                                      "5_อื่นๆ/ไม่ระบุ")),
    marital_group = factor(as.character(marital_group),
                           levels = c("1_Married","2_Single",
                                      "3_Formerly_Married","9_Other_or_NS")),
    inscl_group   = factor(as.character(inscl_group),
                           levels = c("1_UCS","2_CSMBS","3_SSS","9_Other_or_NS")),

    log_cost         = log(cost_w + 1),
    staff_per_1000_z = as.numeric(scale(staff_per_1000))
  )

cat("   n หลัง filter:", nrow(dat), "\n")

# ── Dataset สำหรับ BP (กรอง NA) ──────────────────────────────
dat_bp <- dat %>% filter(!is.na(bp_controlled))
cat("   n BP analysis:", nrow(dat_bp), "\n\n")

# ── Covariates (ตรงกับ Part 3) ────────────────────────────────
covariates       <- c("age","sex_bin","inscl_group","edu_group",
                      "occ_group5","marital_group","comorbidity_count")
covariate_center <- "staff_per_1000_z"

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================
icc_logit <- function(var_cluster) {
  var_cluster / (var_cluster + (pi^2 / 3))
}

get_robust_or <- function(model, data) {
  ct  <- coeftest(model, vcov = vcovCL(model, cluster = ~hcode, data = data))
  est <- exp(coef(model))
  ci  <- suppressMessages(exp(confint(model)))
  data.frame(
    variable  = names(est),
    OR        = round(est, 4),
    CI_lower  = round(ci[, 1], 4),
    CI_upper  = round(ci[, 2], 4),
    p_value   = round(ct[, "Pr(>|z|)"], 4),
    sig       = ifelse(ct[, "Pr(>|z|)"] < 0.001, "***",
                ifelse(ct[, "Pr(>|z|)"] < 0.01,  "**",
                ifelse(ct[, "Pr(>|z|)"] < 0.05,  "*", "ns"))),
    stringsAsFactors = FALSE
  )
}

get_robust_pct <- function(model, data) {
  ct  <- coeftest(model, vcov = vcovCL(model, cluster = ~hcode, data = data))
  est <- coef(model)
  ci  <- suppressMessages(confint(model))
  data.frame(
    variable   = names(est),
    beta       = round(est, 4),
    pct_change = round((exp(est) - 1) * 100, 2),
    CI_lower   = round((exp(ci[, 1]) - 1) * 100, 2),
    CI_upper   = round((exp(ci[, 2]) - 1) * 100, 2),
    p_value    = round(ct[, "Pr(>|t|)"], 4),
    sig        = ifelse(ct[, "Pr(>|t|)"] < 0.001, "***",
                 ifelse(ct[, "Pr(>|t|)"] < 0.01,  "**",
                 ifelse(ct[, "Pr(>|t|)"] < 0.05,  "*", "ns"))),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 3. โมเดลที่ 1: CONTINUITY (ตรงกับ STEP D2)
# ============================================================
cat("📊 Model 1: Continuity...\n")

# Crude
m_crude_con <- glm(continuity_flag01 ~ out_of_catchment,
                   data = dat, family = binomial())
res_crude_con <- get_robust_or(m_crude_con, dat)

# Adjusted
f_adj_con <- as.formula(paste(
  "continuity_flag01 ~ out_of_catchment +",
  paste(covariates, collapse = " + "), "+", covariate_center
))
m_adj_con <- glm(f_adj_con, data = dat, family = binomial())
res_adj_con <- get_robust_or(m_adj_con, dat)

# AUC
roc_con <- pROC::roc(dat$continuity_flag01, fitted(m_adj_con), quiet = TRUE)
auc_con <- round(as.numeric(pROC::auc(roc_con)), 3)
cat("   AUC:", auc_con, "\n")

# Hosmer-Lemeshow
hl_con <- tryCatch(
  ResourceSelection::hoslem.test(dat$continuity_flag01, fitted(m_adj_con), g = 10),
  error = function(e) list(p.value = NA)
)

# ICC (null model)
null_con <- tryCatch({
  lme4::glmer(continuity_flag01 ~ 1 + (1 | hcode), data = dat,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
}, error = function(e) NULL)
icc_con <- tryCatch(
  icc_logit(as.numeric(VarCorr(null_con)$hcode)),
  error = function(e) NA
)
cat("   ICC:", round(icc_con, 4), "\n")

# Sensitivity: Multilevel adjusted
ml_con <- tryCatch({
  f_ml <- as.formula(paste(
    "continuity_flag01 ~ out_of_catchment +",
    paste(covariates, collapse = " + "), "+", covariate_center, "+ (1|hcode)"
  ))
  lme4::glmer(f_ml, data = dat, family = binomial,
              control = glmerControl(optimizer = "bobyqa"))
}, error = function(e) { cat("   ⚠️  ML error\n"); NULL })
or_ml_con <- tryCatch(
  round(exp(fixef(ml_con)["out_of_catchmentOut_of_Catchment"]), 4),
  error = function(e) NA
)
cat("   ML OR (OOC):", or_ml_con, "\n")

# ============================================================
# 4. โมเดลที่ 2: BP CONTROL (ตรงกับ STEP D)
# ============================================================
cat("\n📊 Model 2: BP Control...\n")

# Crude
m_crude_bp <- glm(bp_controlled ~ out_of_catchment,
                  data = dat_bp, family = binomial())
res_crude_bp <- get_robust_or(m_crude_bp, dat_bp)

# Adjusted (มี continuity_flag เป็น covariate)
f_adj_bp <- as.formula(paste(
  "bp_controlled ~ out_of_catchment + continuity_flag +",
  paste(covariates, collapse = " + "), "+", covariate_center
))
m_adj_bp <- glm(f_adj_bp, data = dat_bp, family = binomial())
res_adj_bp <- get_robust_or(m_adj_bp, dat_bp)

# AUC
roc_bp <- pROC::roc(dat_bp$bp_controlled, fitted(m_adj_bp), quiet = TRUE)
auc_bp <- round(as.numeric(pROC::auc(roc_bp)), 3)
cat("   AUC:", auc_bp, "\n")

# Hosmer-Lemeshow
hl_bp <- tryCatch(
  ResourceSelection::hoslem.test(dat_bp$bp_controlled, fitted(m_adj_bp), g = 10),
  error = function(e) list(p.value = NA)
)

# ICC
null_bp <- tryCatch({
  lme4::glmer(bp_controlled ~ 1 + (1 | hcode), data = dat_bp,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
}, error = function(e) NULL)
icc_bp <- tryCatch(
  icc_logit(as.numeric(VarCorr(null_bp)$hcode)),
  error = function(e) NA
)
cat("   ICC:", round(icc_bp, 4), "\n")

# ============================================================
# 5. โมเดลที่ 3: LOG COST (ตรงกับ STEP E)
# ============================================================
cat("\n📊 Model 3: Healthcare Cost...\n")

# Crude
m_crude_cost <- lm(log_cost ~ out_of_catchment, data = dat)
res_crude_cost <- get_robust_pct(m_crude_cost, dat)

# Adjusted (มี continuity_flag + n_visits)
f_adj_cost <- as.formula(paste(
  "log_cost ~ out_of_catchment + continuity_flag +",
  paste(covariates, collapse = " + "), "+ n_visits +", covariate_center
))
lm_adj <- lm(f_adj_cost, data = dat)
res_adj_cost <- get_robust_pct(lm_adj, dat)

r2_cost  <- round(summary(lm_adj)$r.squared, 4)
adjr2    <- round(summary(lm_adj)$adj.r.squared, 4)
cat("   R²:", r2_cost, "| Adj R²:", adjr2, "\n")

# Gamma GLM sensitivity
gamma_res <- tryCatch({
  dat_pos <- dat %>% filter(cost_w > 0)
  f_g <- as.formula(paste(
    "cost_w ~ out_of_catchment + continuity_flag +",
    paste(covariates, collapse = " + "), "+ n_visits +", covariate_center
  ))
  m_g   <- glm(f_g, data = dat_pos, family = Gamma(link = "log"))
  ct_g  <- coeftest(m_g, vcov = vcovCL(m_g, cluster = ~hcode, data = dat_pos))
  round((exp(ct_g["out_of_catchmentOut_of_Catchment", "Estimate"]) - 1) * 100, 2)
}, error = function(e) NA)
cat("   Gamma OOC %change:", gamma_res, "\n")

# ============================================================
# 6. สรุปตามกลุ่ม
# ============================================================
cat("\n📊 สรุปตามกลุ่ม...\n")

group_summary <- dat %>%
  group_by(out_of_catchment) %>%
  summarise(
    n              = n(),
    continuity_pct = round(mean(continuity_flag01, na.rm = TRUE) * 100, 2),
    bp_pct         = round(mean(bp_controlled, na.rm = TRUE) * 100, 2),
    cost_median    = round(median(cost_w, na.rm = TRUE), 0),
    cost_mean      = round(mean(cost_w, na.rm = TRUE), 0),
    age_mean       = round(mean(age, na.rm = TRUE), 1),
    .groups = "drop"
  )

group3_summary <- dat %>%
  group_by(ooc_subgroup3) %>%
  summarise(
    n              = n(),
    continuity_pct = round(mean(continuity_flag01, na.rm = TRUE) * 100, 2),
    bp_pct         = round(mean(bp_controlled, na.rm = TRUE) * 100, 2),
    cost_median    = round(median(cost_w, na.rm = TRUE), 0),
    cost_mean      = round(mean(cost_w, na.rm = TRUE), 0),
    age_mean       = round(mean(age, na.rm = TRUE), 1),
    .groups = "drop"
  )

# ============================================================
# 7. สรุปรายศูนย์ 69 ศบส. (ตรงกับ F_Part4_CenterOOC.R)
# ============================================================
cat("\n📊 สรุปรายศูนย์...\n")

# ── โหลด GIS: lat, lon, hc_name จาก Part 1 ──────────────────
gis_data <- tryCatch({
  readRDS("data/hc_resources_clean.rds") %>%
    dplyr::select(hcode, hc_name, lat, lon) %>%
    dplyr::mutate(hcode = as.character(hcode))
}, error = function(e) {
  cat("   ⚠️  ไม่พบ data/hc_resources_clean.rds — ข้ามการ merge GIS\n")
  NULL
})

# ── District mapping (กลุ่มเขต 6 กลุ่ม กทม.) ─────────────────
district_group_map <- tibble::tribble(
  ~hc_name_match,       ~district,          ~district_group,
  # กลุ่มเขตกรุงเทพกลาง
  "ศบส.1",  "พระนคร",          "กรุงเทพกลาง",
  "ศบส.2",  "ป้อมปราบศัตรูพ่าย","กรุงเทพกลาง",
  "ศบส.3",  "สัมพันธวงศ์",      "กรุงเทพกลาง",
  "ศบส.4",  "บางรัก",           "กรุงเทพกลาง",
  "ศบส.5",  "สาทร",             "กรุงเทพกลาง",
  "ศบส.6",  "บางคอแหลม",        "กรุงเทพใต้",
  "ศบส.7",  "ยานนาวา",          "กรุงเทพใต้",
  "ศบส.8",  "บึงกุ่ม",          "กรุงเทพตะวันออก",
  "ศบส.9",  "ลาดกระบัง",        "กรุงเทพตะวันออก",
  "ศบส.10", "มีนบุรี",          "กรุงเทพตะวันออก",
  "ศบส.11", "หนองจอก",          "กรุงเทพตะวันออก",
  "ศบส.12", "คลองสาน",          "กรุงเทพใต้",
  "ศบส.13", "ธนบุรี",           "กรุงธนเหนือ",
  "ศบส.14", "บางกอกน้อย",       "กรุงธนเหนือ",
  "ศบส.15", "บางกอกใหญ่",       "กรุงธนเหนือ",
  "ศบส.16", "ภาษีเจริญ",        "กรุงธนเหนือ",
  "ศบส.17", "หนองแขม",          "กรุงธนใต้",
  "ศบส.18", "ราษฎร์บูรณะ",      "กรุงธนใต้",
  "ศบส.19", "จอมทอง",           "กรุงธนใต้",
  "ศบส.20", "บางแค",            "กรุงธนใต้",
  "ศบส.21", "ทวีวัฒนา",         "กรุงธนใต้",
  "ศบส.22", "ตลิ่งชัน",         "กรุงธนเหนือ",
  "ศบส.23", "บางพลัด",          "กรุงธนเหนือ",
  "ศบส.24", "ดุสิต",            "กรุงเทพกลาง",
  "ศบส.25", "พญาไท",            "กรุงเทพกลาง",
  "ศบส.26", "ราชเทวี",          "กรุงเทพกลาง",
  "ศบส.27", "ดินแดง",           "กรุงเทพกลาง",
  "ศบส.28", "ห้วยขวาง",         "กรุงเทพกลาง",
  "ศบส.29", "วังทองหลาง",       "กรุงเทพตะวันออก",
  "ศบส.30", "สวนหลวง",          "กรุงเทพใต้",
  "ศบส.31", "ประเวศ",           "กรุงเทพใต้",
  "ศบส.32", "พระโขนง",          "กรุงเทพใต้",
  "ศบส.33", "สาทร",             "กรุงเทพกลาง",
  "ศบส.34", "บางนา",            "กรุงเทพใต้",
  "ศบส.35", "สะพานสูง",         "กรุงเทพตะวันออก",
  "ศบส.36", "คันนายาว",         "กรุงเทพตะวันออก",
  "ศบส.37", "ลาดพร้าว",         "กรุงเทพตะวันออก",
  "ศบส.38", "บางเขน",           "กรุงเทพเหนือ",
  "ศบส.39", "สายไหม",           "กรุงเทพเหนือ",
  "ศบส.40", "ดอนเมือง",         "กรุงเทพเหนือ",
  "ศบส.41", "หลักสี่",          "กรุงเทพเหนือ",
  "ศบส.42", "ลาดพร้าว",         "กรุงเทพตะวันออก",
  "ศบส.43", "บึงกุ่ม",          "กรุงเทพตะวันออก",
  "ศบส.44", "บางกะปิ",          "กรุงเทพตะวันออก",
  "ศบส.45", "สาทร",             "กรุงเทพกลาง",
  "ศบส.46", "พระโขนง",          "กรุงเทพใต้",
  "ศบส.47", "คลองเตย",          "กรุงเทพใต้",
  "ศบส.48", "วัฒนา",            "กรุงเทพใต้",
  "ศบส.49", "จตุจักร",          "กรุงเทพเหนือ",
  "ศบส.50", "จตุจักร",          "กรุงเทพเหนือ",
  "ศบส.51", "บางซื่อ",          "กรุงเทพเหนือ",
  "ศบส.52", "บางเขน",           "กรุงเทพเหนือ",
  "ศบส.53", "ลาดกระบัง",        "กรุงเทพตะวันออก",
  "ศบส.54", "มีนบุรี",          "กรุงเทพตะวันออก",
  "ศบส.55", "หนองจอก",          "กรุงเทพตะวันออก",
  "ศบส.56", "คลองสามวา",        "กรุงเทพตะวันออก",
  "ศบส.57", "บางเขน",           "กรุงเทพเหนือ",
  "ศบส.58", "ทุ่งครุ",          "กรุงธนใต้",
  "ศบส.59", "บางบอน",           "กรุงธนใต้",
  "ศบส.60", "บางขุนเทียน",      "กรุงธนใต้",
  "ศบส.61", "หนองแขม",          "กรุงธนใต้",
  "ศบส.62", "บางแค",            "กรุงธนใต้",
  "ศบส.63", "บางแค",            "กรุงธนใต้",
  "ศบส.64", "ลาดกระบัง",        "กรุงเทพตะวันออก",
  "ศบส.65", "ประเวศ",           "กรุงเทพใต้",
  "ศบส.66", "สาทร",             "กรุงเทพกลาง",
  "ศบส.67", "บางรัก",           "กรุงเทพกลาง",
  "ศบส.68", "ดอนเมือง",         "กรุงเทพเหนือ",
  "ศบส.69", "สายไหม",           "กรุงเทพเหนือ"
)

center_summary <- dat %>%
  group_by(hcode) %>%
  summarise(
    n              = n(),
    ooc_pct        = round(mean(out_of_catchment == "Out_of_Catchment") * 100, 2),
    ooc_bkk_pct   = round(mean(ooc_subgroup3 == "OOC_Bangkok", na.rm = TRUE) * 100, 2),
    ooc_out_pct   = round(mean(ooc_subgroup3 == "OOC_Outside", na.rm = TRUE) * 100, 2),
    continuity_pct = round(mean(continuity_flag01, na.rm = TRUE) * 100, 2),
    bp_pct         = round(mean(bp_controlled, na.rm = TRUE) * 100, 2),
    cost_median    = round(median(cost_w, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  mutate(
    ooc_level = case_when(
      ooc_pct >= 50 ~ "สูงมาก (≥50%)",
      ooc_pct >= 30 ~ "ปานกลาง (30–49%)",
      TRUE          ~ "ต่ำ (<30%)"
    )
  ) %>%
  arrange(hcode)

# ── Merge GIS (lat/lon/hc_name) ──────────────────────────────
if (!is.null(gis_data)) {
  center_summary <- center_summary %>%
    dplyr::left_join(gis_data, by = "hcode")
  cat("   ✅ Merge lat/lon/hc_name สำเร็จ\n")
  cat("   มีพิกัด:", sum(!is.na(center_summary$lat)), "จาก", nrow(center_summary), "ศูนย์\n")
} else {
  center_summary <- center_summary %>%
    dplyr::mutate(lat = NA_real_, lon = NA_real_, hc_name = NA_character_)
}

# ── Merge district/district_group จาก mapping table ──────────
# ดึงหมายเลข ศบส. จาก hcode แล้ว match กับ district_group_map
center_summary <- center_summary %>%
  dplyr::mutate(
    hcode_num = as.integer(gsub("[^0-9]", "", hcode)),
    hc_name_match = paste0("ศบส.", hcode_num)
  ) %>%
  dplyr::left_join(
    district_group_map %>% dplyr::select(hc_name_match, district, district_group),
    by = "hc_name_match"
  ) %>%
  dplyr::select(-hc_name_match, -hcode_num)

cat("   จำนวนศูนย์:", nrow(center_summary), "\n")
cat("   มี district:", sum(!is.na(center_summary$district)), "\n")

# ============================================================
# 8. KPI summary
# ============================================================
cat("\n📊 สร้าง KPI...\n")

# ดึง OR/pct จาก adjusted model สำหรับ OOC
or_ooc_con  <- res_adj_con[res_adj_con$variable == "out_of_catchmentOut_of_Catchment", ]
or_ooc_bp   <- res_adj_bp[res_adj_bp$variable   == "out_of_catchmentOut_of_Catchment", ]
pct_ooc_cost <- res_adj_cost[res_adj_cost$variable == "out_of_catchmentOut_of_Catchment", ]
or_cont_bp  <- res_adj_bp[res_adj_bp$variable   == "continuity_flag1_Continuous", ]
pct_cont_cost <- res_adj_cost[res_adj_cost$variable == "continuity_flag1_Continuous", ]

kpi <- list(
  # ขนาดตัวอย่าง
  total_n         = nrow(dat),
  n_centers       = n_distinct(dat$hcode),
  n_bp_analysis   = nrow(dat_bp),
  ooc_n           = sum(dat$out_of_catchment == "Out_of_Catchment"),
  ooc_pct         = round(mean(dat$out_of_catchment == "Out_of_Catchment") * 100, 1),
  ooc_bkk_n       = sum(dat$ooc_subgroup3 == "OOC_Bangkok",  na.rm = TRUE),
  ooc_bkk_pct     = round(mean(dat$ooc_subgroup3 == "OOC_Bangkok", na.rm = TRUE) * 100, 1),
  ooc_out_n       = sum(dat$ooc_subgroup3 == "OOC_Outside", na.rm = TRUE),
  ooc_out_pct     = round(mean(dat$ooc_subgroup3 == "OOC_Outside", na.rm = TRUE) * 100, 1),

  # ผลลัพธ์รวม
  continuity_pct  = round(mean(dat$continuity_flag01, na.rm = TRUE) * 100, 1),
  bp_pct          = round(mean(dat_bp$bp_controlled, na.rm = TRUE) * 100, 1),
  cost_mean_all   = round(mean(dat$cost_w, na.rm = TRUE), 0),
  cost_median_all = round(median(dat$cost_w, na.rm = TRUE), 0),

  # Regression: OOC → Continuity
  or_ooc_cont      = if (nrow(or_ooc_con) > 0) or_ooc_con$OR       else NA,
  or_ooc_cont_low  = if (nrow(or_ooc_con) > 0) or_ooc_con$CI_lower else NA,
  or_ooc_cont_high = if (nrow(or_ooc_con) > 0) or_ooc_con$CI_upper else NA,
  p_ooc_cont       = if (nrow(or_ooc_con) > 0) or_ooc_con$p_value  else NA,
  sig_ooc_cont     = if (nrow(or_ooc_con) > 0) or_ooc_con$sig      else NA,

  # Regression: OOC → BP
  or_ooc_bp        = if (nrow(or_ooc_bp) > 0) or_ooc_bp$OR       else NA,
  or_ooc_bp_low    = if (nrow(or_ooc_bp) > 0) or_ooc_bp$CI_lower else NA,
  or_ooc_bp_high   = if (nrow(or_ooc_bp) > 0) or_ooc_bp$CI_upper else NA,
  p_ooc_bp         = if (nrow(or_ooc_bp) > 0) or_ooc_bp$p_value  else NA,
  sig_ooc_bp       = if (nrow(or_ooc_bp) > 0) or_ooc_bp$sig      else NA,

  # Regression: Continuity → BP
  or_cont_bp       = if (nrow(or_cont_bp) > 0) or_cont_bp$OR       else NA,
  or_cont_bp_low   = if (nrow(or_cont_bp) > 0) or_cont_bp$CI_lower else NA,
  or_cont_bp_high  = if (nrow(or_cont_bp) > 0) or_cont_bp$CI_upper else NA,
  p_cont_bp        = if (nrow(or_cont_bp) > 0) or_cont_bp$p_value  else NA,
  sig_cont_bp      = if (nrow(or_cont_bp) > 0) or_cont_bp$sig      else NA,

  # Regression: OOC → Cost
  pct_ooc_cost     = if (nrow(pct_ooc_cost) > 0) pct_ooc_cost$pct_change else NA,
  pct_ooc_cost_low = if (nrow(pct_ooc_cost) > 0) pct_ooc_cost$CI_lower   else NA,
  pct_ooc_cost_high= if (nrow(pct_ooc_cost) > 0) pct_ooc_cost$CI_upper   else NA,
  p_ooc_cost       = if (nrow(pct_ooc_cost) > 0) pct_ooc_cost$p_value    else NA,
  sig_ooc_cost     = if (nrow(pct_ooc_cost) > 0) pct_ooc_cost$sig        else NA,

  # Regression: Continuity → Cost
  pct_cont_cost     = if (nrow(pct_cont_cost) > 0) pct_cont_cost$pct_change else NA,
  pct_cont_cost_low = if (nrow(pct_cont_cost) > 0) pct_cont_cost$CI_lower   else NA,
  pct_cont_cost_high= if (nrow(pct_cont_cost) > 0) pct_cont_cost$CI_upper   else NA,
  p_cont_cost       = if (nrow(pct_cont_cost) > 0) pct_cont_cost$p_value    else NA,
  sig_cont_cost     = if (nrow(pct_cont_cost) > 0) pct_cont_cost$sig        else NA,

  # Model diagnostics
  auc_con     = auc_con,
  auc_bp      = auc_bp,
  hl_p_con    = round(hl_con$p.value, 4),
  hl_p_bp     = round(hl_bp$p.value, 4),
  icc_con     = round(icc_con, 4),
  icc_bp      = round(icc_bp, 4),
  or_ml_con   = or_ml_con,
  gamma_ooc_pct = gamma_res,
  r2_cost     = r2_cost,
  adjr2_cost  = adjr2,

  # Meta
  fiscal_year  = "2567",
  n_centers_high_ooc = sum(center_summary$ooc_level == "สูงมาก (≥50%)"),
  n_centers_mid_ooc  = sum(center_summary$ooc_level == "ปานกลาง (30–49%)"),
  n_centers_low_ooc  = sum(center_summary$ooc_level == "ต่ำ (<30%)"),
  exported_at  = format(Sys.time(), "%Y-%m-%d %H:%M")
)

# ============================================================
# 9. รวม JSON และ export
# ============================================================
cat("\n💾 Export JSON...\n")

output <- list(
  kpi            = kpi,
  group_summary  = group_summary,
  group3_summary = group3_summary,
  center_summary = center_summary,
  regression     = list(
    continuity_crude = res_crude_con,
    continuity_adj   = res_adj_con,
    bp_crude         = res_crude_bp,
    bp_adj           = res_adj_bp,
    cost_crude       = res_crude_cost,
    cost_adj         = res_adj_cost
  )
)

jsonlite::write_json(output, path = OUTPUT_JSON,
                     auto_unbox = TRUE, pretty = TRUE, na = "null")

cat("✅ Export สำเร็จ:", OUTPUT_JSON, "\n")
cat("   ขนาดไฟล์:", round(file.size(OUTPUT_JSON) / 1024, 1), "KB\n")

cat("\n📋 KPI สรุป:\n")
cat(sprintf("   ผู้รับบริการ        : %d ราย\n",     kpi$total_n))
cat(sprintf("   OOC                : %.1f%%\n",       kpi$ooc_pct))
cat(sprintf("   ความสม่ำเสมอ       : %.1f%%\n",      kpi$continuity_pct))
cat(sprintf("   ควบคุม BP          : %.1f%%\n",       kpi$bp_pct))
cat(sprintf("   OR OOC→Continuity  : %.3f (p=%s)\n",  kpi$or_ooc_cont,  kpi$sig_ooc_cont))
cat(sprintf("   OR OOC→BP          : %.3f (p=%s)\n",  kpi$or_ooc_bp,    kpi$sig_ooc_bp))
cat(sprintf("   %%change OOC→Cost  : %.1f%% (p=%s)\n", kpi$pct_ooc_cost, kpi$sig_ooc_cost))

cat("\n🎯 นำ", OUTPUT_JSON, "ไปวางใน GitHub repo แทนของเก่า\n")
cat("=========================================\n")
cat("  EXPORT COMPLETE\n")
cat("=========================================\n")
