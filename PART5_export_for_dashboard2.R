# ============================================================
# PART5_export_for_dashboard.R
# patch: แก้ hcode + เขต/กลุ่มเขต + รองรับ UI
# ============================================================

setwd("C:/Users/ADMin/Desktop/HT Project")

library(dplyr)
library(readr)
library(stringr)
library(janitor)
library(jsonlite)
library(lme4)
library(sandwich)
library(lmtest)
library(performance)

# ============================================================
# 0. CONFIG
# ============================================================
OUTPUT_JSON <- "dashboard_data.json"

# ============================================================
# helper
# ============================================================
normalize_hcode <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x <- str_remove(x, "^ศบส\\.?\\s*")
  x <- str_extract(x, "\\d+")
  x <- ifelse(is.na(x), NA_character_, str_pad(x, width = 2, side = "left", pad = "0"))
  x
}

pretty_hcode <- function(x) {
  x_num <- suppressWarnings(as.integer(x))
  ifelse(is.na(x_num), NA_character_, paste0("ศบส.", x_num))
}

get_robust_df <- function(model, data, type = "logistic") {
  if (type == "logistic") {
    ct   <- coeftest(model, vcov = vcovCL(model, cluster = ~hcode, data = data))
    est  <- exp(coef(model))
    ci   <- exp(suppressMessages(confint(model)))
    pval <- ct[, "Pr(>|z|)"]
  } else {
    ct   <- coeftest(model, vcov = vcovCL(model, cluster = ~hcode, data = data))
    est  <- coef(model)
    ci   <- suppressMessages(confint(model))
    pval <- ct[, "Pr(>|t|)"]
  }
  
  data.frame(
    variable = names(est),
    estimate = round(est, 4),
    ci_lower = round(ci[, 1], 4),
    ci_upper = round(ci[, 2], 4),
    p_value  = round(pval, 4),
    sig      = ifelse(pval < 0.001, "***",
                      ifelse(pval < 0.01,  "**",
                             ifelse(pval < 0.05,  "*", ""))),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 1. LOAD DATA
# ============================================================
cat("📂 กำลังโหลดข้อมูล...\n")
d <- readRDS("checkpoint/DM_FINAL_ready_for_analysis.rds") %>%
  mutate(
    hcode_raw = as.character(hcode),
    hcode = normalize_hcode(hcode)
  )

cat("   จำนวนแถว:", nrow(d), "| คอลัมน์:", ncol(d), "\n")

# ============================================================
# 2. PREP VARIABLES
# ============================================================
d$out_of_catchment <- factor(
  d$out_of_catchment,
  levels = c("In_Catchment", "Out_of_Catchment"),
  labels = c("In-Catchment", "Out-of-Catchment")
)

d$continuity_flag <- factor(
  d$continuity_flag,
  levels = c("0_NotContinuous", "1_Continuous"),
  labels = c("0", "1")
)

d$bp_controlled <- factor(d$bp_controlled, levels = c(0, 1))

d$ooc_group <- case_when(
  d$out_of_catchment == "In-Catchment" ~ "In-Catchment",
  d$ooc_subgroup == "OOC_Bangkok"      ~ "OOC-Bangkok",
  d$ooc_subgroup == "OOC_Outside"      ~ "OOC-Outside",
  TRUE                                 ~ "OOC-Outside"
)

# Winsorize cost ที่ P99
p99 <- quantile(d$cost_w, 0.99, na.rm = TRUE)
d$cost_w <- pmin(d$cost_w, p99)
d$log_cost <- log(d$cost_w + 1)

# standardize staff
if ("staff_per_1000" %in% names(d)) {
  d$staff_per_1000_z <- as.numeric(scale(d$staff_per_1000))
}

cat("   OOC rate     :", round(mean(d$out_of_catchment == "Out-of-Catchment") * 100, 1), "%\n")
cat("   Continuity   :", round(mean(d$continuity_flag == "1", na.rm = TRUE) * 100, 1), "%\n")
cat("   BP controlled:", round(mean(as.numeric(as.character(d$bp_controlled)), na.rm = TRUE) * 100, 1), "%\n")

# ============================================================
# 3. COVARIATES
# ============================================================
base_cov  <- c("age", "sex_bin", "edu_group", "marital_group",
               "occ_group5", "inscl_group", "staff_per_1000_z",
               "comorbidity_count")
cov_avail <- base_cov[base_cov %in% names(d)]
cov_str   <- paste(cov_avail, collapse = " + ")

cat("   Covariates:", cov_str, "\n")

# ============================================================
# 4. MODELS
# ============================================================
cat("📊 ผลลัพธ์ที่ 1: continuity_flag\n")

m_cont_1 <- glm(continuity_flag ~ out_of_catchment,
                data = d, family = binomial())

f_cont_2 <- as.formula(paste("continuity_flag ~ out_of_catchment +", cov_str))
m_cont_2 <- glm(f_cont_2, data = d, family = binomial())

res_cont_1 <- get_robust_df(m_cont_1, d, "logistic")
res_cont_2 <- get_robust_df(m_cont_2, d, "logistic")

auc_cont <- tryCatch(
  as.numeric(performance::performance_roc(m_cont_2)$AUC),
  error = function(e) NA
)

m_cont_ml <- tryCatch({
  f_ml <- as.formula(paste("continuity_flag ~ out_of_catchment +", cov_str, "+ (1|hcode)"))
  glmer(f_ml, data = d, family = binomial(),
        control = glmerControl(optimizer = "bobyqa"))
}, error = function(e) NULL)

icc_cont <- tryCatch(
  round(as.numeric(performance::icc(m_cont_ml)$ICC_adjusted), 4),
  error = function(e) NA
)

or_ml_cont <- tryCatch(
  round(exp(fixef(m_cont_ml)["out_of_catchmentOut-of-Catchment"]), 4),
  error = function(e) NA
)

cat("   AUC (Model 2):", round(auc_cont, 3), "\n")
cat("   ICC:", icc_cont, "| OR (multilevel):", or_ml_cont, "\n\n")

# ------------------------------------------------------------
cat("📊 ผลลัพธ์ที่ 2: bp_controlled\n")

d_bp <- d %>% filter(!is.na(bp_controlled))
cat("   n:", nrow(d_bp), "\n")

m_bp_1 <- glm(bp_controlled ~ out_of_catchment,
              data = d_bp, family = binomial())

f_bp_2 <- as.formula(paste("bp_controlled ~ out_of_catchment + continuity_flag +", cov_str))
m_bp_2 <- glm(f_bp_2, data = d_bp, family = binomial())

res_bp_1 <- get_robust_df(m_bp_1, d_bp, "logistic")
res_bp_2 <- get_robust_df(m_bp_2, d_bp, "logistic")

auc_bp <- tryCatch(
  as.numeric(performance::performance_roc(m_bp_2)$AUC),
  error = function(e) NA
)

cat("   AUC (Model 2):", round(auc_bp, 3), "\n\n")

# ------------------------------------------------------------
cat("📊 ผลลัพธ์ที่ 3: log_cost\n")

cov_cost <- paste(
  c(cov_avail, "n_visits")[c(cov_avail, "n_visits") %in% names(d)],
  collapse = " + "
)

f_cost <- as.formula(paste("log_cost ~ out_of_catchment +", cov_cost))
m_cost <- lm(f_cost, data = d)

res_cost <- get_robust_df(m_cost, d, "linear")
res_cost$pct_change <- round((exp(res_cost$estimate) - 1) * 100, 2)

r2_cost <- round(summary(m_cost)$r.squared, 4)
adjr2   <- round(summary(m_cost)$adj.r.squared, 4)

cat("   R²:", r2_cost, "| Adj R²:", adjr2, "\n\n")

# ============================================================
# 5. SUMMARY
# ============================================================
cat("📊 สรุปตามกลุ่ม...\n")

group_summary <- d %>%
  group_by(out_of_catchment) %>%
  summarise(
    n              = n(),
    continuity_pct = round(mean(continuity_flag == "1", na.rm = TRUE) * 100, 2),
    bp_pct         = round(mean(as.numeric(as.character(bp_controlled)), na.rm = TRUE) * 100, 2),
    cost_median    = round(median(cost_w, na.rm = TRUE), 0),
    cost_mean      = round(mean(cost_w, na.rm = TRUE), 0),
    age_mean       = round(mean(age, na.rm = TRUE), 1),
    .groups = "drop"
  )

group3_summary <- d %>%
  group_by(ooc_group) %>%
  summarise(
    n              = n(),
    continuity_pct = round(mean(continuity_flag == "1", na.rm = TRUE) * 100, 2),
    bp_pct         = round(mean(as.numeric(as.character(bp_controlled)), na.rm = TRUE) * 100, 2),
    cost_median    = round(median(cost_w, na.rm = TRUE), 0),
    cost_mean      = round(mean(cost_w, na.rm = TRUE), 0),
    age_mean       = round(mean(age, na.rm = TRUE), 1),
    .groups = "drop"
  )

# ============================================================
# 6. HC MAP
# ============================================================
cat("📊 สรุปรายศูนย์...\n")

hc_map <- read_csv("data/HC_locations.csv",
                   locale = locale(encoding = "UTF-8"),
                   show_col_types = FALSE) %>%
  clean_names() %>%
  transmute(
    hcode = normalize_hcode(hcode),
    district = str_squish(str_remove(as.character(hc_district_unit), "^เขต")),
    district_group = str_squish(as.character(hc_district_group))
  ) %>%
  distinct(hcode, .keep_all = TRUE)

center_summary <- d %>%
  group_by(hcode) %>%
  summarise(
    n              = n(),
    ooc_pct        = round(mean(out_of_catchment == "Out-of-Catchment") * 100, 2),
    continuity_pct = round(mean(continuity_flag == "1", na.rm = TRUE) * 100, 2),
    bp_pct         = round(mean(as.numeric(as.character(bp_controlled)), na.rm = TRUE) * 100, 2),
    cost_median    = round(median(cost_w, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  left_join(hc_map, by = "hcode") %>%
  mutate(display_hcode = pretty_hcode(hcode)) %>%
  arrange(hcode)

# ใช้ต่อสำหรับเปรียบเทียบ 3 กลุ่มของศูนย์ที่เลือก
center_group3_summary <- d %>%
  group_by(hcode, ooc_group) %>%
  summarise(
    n              = n(),
    continuity_pct = round(mean(continuity_flag == "1", na.rm = TRUE) * 100, 2),
    bp_pct         = round(mean(as.numeric(as.character(bp_controlled)), na.rm = TRUE) * 100, 2),
    cost_median    = round(median(cost_w, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  left_join(hc_map, by = "hcode") %>%
  mutate(display_hcode = pretty_hcode(hcode)) %>%
  arrange(hcode, ooc_group)

cat("   จำนวนศูนย์:", nrow(center_summary), "\n")
cat("   Missing district:", sum(is.na(center_summary$district)), "\n")
cat("   Missing district_group:", sum(is.na(center_summary$district_group)), "\n")
print(head(center_summary %>% select(hcode, display_hcode, district, district_group), 10))

# ============================================================
# 7. KPI
# ============================================================
ooc_var  <- "out_of_catchmentOut-of-Catchment"
rb_cont  <- get_robust_df(m_cont_2, d, "logistic")
rb_bp    <- get_robust_df(m_bp_2, d_bp, "logistic")
rb_cost  <- get_robust_df(m_cost, d, "linear")

kpi <- list(
  total_n         = nrow(d),
  n_centers       = n_distinct(d$hcode),
  ooc_n           = sum(d$out_of_catchment == "Out-of-Catchment"),
  ooc_pct         = round(mean(d$out_of_catchment == "Out-of-Catchment") * 100, 1),
  continuity_pct  = round(mean(d$continuity_flag == "1", na.rm = TRUE) * 100, 1),
  bp_pct          = round(mean(as.numeric(as.character(d$bp_controlled)), na.rm = TRUE) * 100, 1),
  cost_median_all = round(median(d$cost_w, na.rm = TRUE), 0),
  or_ooc_cont     = round(exp(coef(m_cont_2)[ooc_var]), 3),
  p_ooc_cont      = round(rb_cont$p_value[rb_cont$variable == ooc_var], 3),
  or_ooc_bp       = round(exp(coef(m_bp_2)[ooc_var]), 3),
  p_ooc_bp        = round(rb_bp$p_value[rb_bp$variable == ooc_var], 3),
  beta_ooc_cost   = round(coef(m_cost)[ooc_var], 3),
  p_ooc_cost      = round(rb_cost$p_value[rb_cost$variable == ooc_var], 3),
  icc_cont        = icc_cont,
  or_ml_cont      = or_ml_cont,
  auc_cont        = round(auc_cont, 3),
  auc_bp          = round(auc_bp, 3),
  r2_cost         = r2_cost,
  adjr2_cost      = adjr2,
  fiscal_year     = "2567",
  exported_at     = format(Sys.time(), "%Y-%m-%d %H:%M")
)

# ============================================================
# 8. EXPORT JSON
# ============================================================
cat("💾 กำลัง export JSON...\n")

output <- list(
  kpi                   = kpi,
  group_summary         = group_summary,
  group3_summary        = group3_summary,
  center_summary        = center_summary,
  center_group3_summary = center_group3_summary,
  regression = list(
    continuity_unadj = res_cont_1,
    continuity_adj   = res_cont_2,
    bp_unadj         = res_bp_1,
    bp_adj           = res_bp_2,
    cost_adj         = res_cost
  )
)

write(toJSON(output, auto_unbox = TRUE, pretty = TRUE, na = "null"),
      file = OUTPUT_JSON)

cat("✅ Export สำเร็จ:", OUTPUT_JSON, "\n")
cat("🎯 อัปโหลด dashboard_data.json ใหม่ขึ้น GitHub ได้เลย\n")