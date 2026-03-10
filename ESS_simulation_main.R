
# ============================================================
# Enrollment Strategy Simulation (ESS) Framework
# Monte Carlo simulation for university enrollment planning
# Includes:
#   - 4-year discounted tuition stream
#   - black-and-white publication figures
#   - OLS, quantile, and interaction models
#   - optimization, stress testing, and sensitivity analysis
#   - tuition scenario analysis
# Author: Subir Hait
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
library(dplyr)
library(ggplot2)
library(quantreg)
library(purrr)
library(tidyr)
library(broom)
library(konfound)

set.seed(12345)

# -----------------------------
# 1. settings
# -----------------------------
n_sim <- 10000
total_budget <- 50000
discount_rate <- 0.05
years <- 4

# Baseline allocation ratios
ratios <- c(0.05, 0.10, 0.15)

# Baseline tuition revenue per enrolled student per year
tuition_per_enrollment <- 10000

# Cost assumptions
fixed_cost_mean <- 20000
fixed_cost_sd <- 1000
variable_cost_mean <- 5000
variable_cost_sd <- 800

# Uncertain marketing inputs
cpl_min <- 150
cpl_max <- 500
cr_min <- 0.02
cr_max <- 0.07

# Output folders
dir.create("outputs", showWarnings = FALSE)
dir.create(file.path("outputs", "figures"), showWarnings = FALSE)
dir.create(file.path("outputs", "tables"), showWarnings = FALSE)

# -----------------------------
# 2. Plot theme
# -----------------------------
theme_ess_bw <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# -----------------------------
# 3. Helper functions
# -----------------------------

# Present value of a discounted tuition stream
calc_tuition_stream <- function(tuition_per_enrollment,
                                discount_rate = 0.05,
                                years = 4) {
  tuition_per_enrollment * sum(1 / (1 + discount_rate)^(0:(years - 1)))
}

# Core simulation engine
simulate_ess <- function(ratio,
                         n = 10000,
                         total_budget = 50000,
                         discount_rate = 0.05,
                         tuition_per_enrollment = 10000,
                         years = 4,
                         fixed_cost_mean = 20000,
                         fixed_cost_sd = 1000,
                         variable_cost_mean = 5000,
                         variable_cost_sd = 800,
                         cpl_min = 150,
                         cpl_max = 500,
                         cr_min = 0.02,
                         cr_max = 0.07,
                         scenario_name = "Base",
                         cpl_multiplier = 1,
                         cr_multiplier = 1,
                         budget_multiplier = 1) {

  if (ratio <= 0 || ratio >= 1) stop("ratio must be between 0 and 1.")
  if (n <= 0) stop("n must be positive.")
  if (years < 1) stop("years must be at least 1.")

  effective_budget <- total_budget * budget_multiplier
  marketing_budget <- effective_budget * ratio

  cpl <- runif(n, min = cpl_min, max = cpl_max) * cpl_multiplier
  cr <- runif(n, min = cr_min, max = cr_max) * cr_multiplier
  cr <- pmin(pmax(cr, 0), 1)

  fixed_cost <- pmax(rnorm(n, mean = fixed_cost_mean, sd = fixed_cost_sd), 0)
  variable_cost <- pmax(rnorm(n, mean = variable_cost_mean, sd = variable_cost_sd), 0)

  leads <- marketing_budget / cpl
  enrollments <- leads * cr

  tuition_stream <- calc_tuition_stream(
    tuition_per_enrollment = tuition_per_enrollment,
    discount_rate = discount_rate,
    years = years
  )

  revenue <- enrollments * tuition_stream
  total_cost <- fixed_cost + marketing_budget + variable_cost

  profit <- revenue - total_cost
  roi <- profit / total_cost

  # Because revenue already uses discounted tuition, NPV equals profit here
  npv <- profit

  data.frame(
    scenario = scenario_name,
    ratio = ratio,
    effective_budget = effective_budget,
    marketing_budget = marketing_budget,
    cpl = cpl,
    cr = cr,
    fixed_cost = fixed_cost,
    variable_cost = variable_cost,
    leads = leads,
    enrollments = enrollments,
    tuition_stream = tuition_stream,
    revenue = revenue,
    total_cost = total_cost,
    profit = profit,
    roi = roi,
    npv = npv
  )
}

summarise_simulation <- function(dat) {
  dat %>%
    summarise(
      mean_roi = mean(roi),
      sd_roi = sd(roi),
      mean_npv = mean(npv),
      mean_enrollments = mean(enrollments),
      mean_leads = mean(leads),
      p_positive_roi = mean(roi > 0),
      p_positive_npv = mean(npv > 0),
      var_5 = as.numeric(quantile(roi, 0.05)),
      cvar_5 = mean(roi[roi <= quantile(roi, 0.05)])
    )
}

# -----------------------------
# 4. Baseline simulation
# -----------------------------
sim_base <- map_dfr(ratios, ~simulate_ess(
  ratio = .x,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
  years = years,
  fixed_cost_mean = fixed_cost_mean,
  fixed_cost_sd = fixed_cost_sd,
  variable_cost_mean = variable_cost_mean,
  variable_cost_sd = variable_cost_sd,
  cpl_min = cpl_min,
  cpl_max = cpl_max,
  cr_min = cr_min,
  cr_max = cr_max,
  scenario_name = "Base"
))

sim_base <- sim_base %>%
  mutate(
    ratio_factor = factor(ratio, levels = c(0.10, 0.15, 0.05)),
    ratio_label = factor(
      paste0(round(ratio * 100), "%"),
      levels = c("5%", "10%", "15%")
    )
  )

# -----------------------------
# 5. Table 1: baseline summary
# -----------------------------
table1 <- sim_base %>%
  group_by(ratio, ratio_label) %>%
  summarise(
    mean_roi = mean(roi),
    sd_roi = sd(roi),
    mean_npv = mean(npv),
    mean_enrollments = mean(enrollments),
    mean_leads = mean(leads),
    .groups = "drop"
  ) %>%
  arrange(ratio)

print(table1)

# -----------------------------
# 6. Figure 1: ROI distribution
# -----------------------------
fig1 <- ggplot(
  sim_base,
  aes(x = roi, y = after_stat(density), fill = ratio_label, color = ratio_label)
) +
  geom_histogram(
    bins = 50,
    alpha = 0.15,
    position = "identity",
    linewidth = 0.2
  ) +
  geom_density(linewidth = 1) +
  scale_fill_grey(start = 0.85, end = 0.35) +
  scale_color_grey(start = 0.15, end = 0.55) +
  labs(
    title = "ROI Distribution Across Marketing Allocation Ratios",
    x = "Return on Investment (ROI)",
    y = "Density",
    fill = "Marketing Allocation",
    color = "Marketing Allocation"
  ) +
  theme_ess_bw()

print(fig1)

# -----------------------------
# 7. Figures 2 and 3: scatterplots
# -----------------------------
set.seed(12345)
sim_sample <- dplyr::sample_n(sim_base, min(3000, nrow(sim_base)))

fig2 <- ggplot(sim_sample, aes(x = cr, y = roi)) +
  geom_point(alpha = 0.25, shape = 16, size = 1, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  labs(
    title = "Conversion Rate and ROI",
    x = "Conversion Rate",
    y = "ROI"
  ) +
  theme_ess_bw()

print(fig2)

fig3 <- ggplot(sim_sample, aes(x = cpl, y = roi)) +
  geom_point(alpha = 0.25, shape = 16, size = 1, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  labs(
    title = "Cost per Lead and ROI",
    x = "Cost per Lead (CPL)",
    y = "ROI"
  ) +
  theme_ess_bw()

print(fig3)

# -----------------------------
# 8. Main econometric models
# -----------------------------
ols_model <- lm(roi ~ cr + cpl + ratio_factor, data = sim_base)
qr_model_25 <- rq(roi ~ cr + cpl + ratio_factor, tau = 0.25, data = sim_base)

ols_table <- broom::tidy(ols_model)
qr_table_25 <- broom::tidy(qr_model_25)

print(summary(ols_model))
print(summary(qr_model_25))
print(ols_table)
print(qr_table_25)

# Additional quantile regressions
taus <- c(0.25, 0.50, 0.75)

qr_models <- lapply(taus, function(tau_value) {
  rq(roi ~ cr + cpl + ratio_factor, tau = tau_value, data = sim_base)
})

qr_results <- bind_rows(
  lapply(seq_along(qr_models), function(i) {
    broom::tidy(qr_models[[i]]) %>%
      mutate(tau = taus[i])
  })
)

print(qr_results)

# -----------------------------
# 9. Interaction models
# -----------------------------
ols_model_int <- lm(roi ~ cr * ratio_factor + cpl * ratio_factor, data = sim_base)
qr_model_int <- rq(roi ~ cr * ratio_factor + cpl * ratio_factor, tau = 0.25, data = sim_base)

print(summary(ols_model_int))
print(summary(qr_model_int))

ols_int_table <- broom::tidy(ols_model_int)
qr_int_table <- broom::tidy(qr_model_int)

print(ols_int_table)
print(qr_int_table)

# -----------------------------
# 10. Figure 4: predicted ROI by CR
# -----------------------------
cr_grid <- expand.grid(
  cr = seq(min(sim_base$cr), max(sim_base$cr), length.out = 100),
  cpl = mean(sim_base$cpl),
  ratio_factor = factor(c(0.10, 0.15, 0.05), levels = c(0.10, 0.15, 0.05))
)

cr_grid$pred_roi <- predict(ols_model_int, newdata = cr_grid)
cr_grid$ratio_label <- factor(
  ifelse(cr_grid$ratio_factor == 0.05, "5%",
         ifelse(cr_grid$ratio_factor == 0.10, "10%", "15%")),
  levels = c("5%", "10%", "15%")
)

fig4 <- ggplot(cr_grid, aes(x = cr, y = pred_roi, linetype = ratio_label)) +
  geom_line(color = "black", linewidth = 1) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(
    title = "Predicted ROI by Conversion Rate Across Allocation Scenarios",
    x = "Conversion Rate",
    y = "Predicted ROI",
    linetype = "Allocation Ratio"
  ) +
  theme_ess_bw()

print(fig4)

# -----------------------------
# 11. Figure 5: predicted ROI by CPL
# -----------------------------
cpl_grid <- expand.grid(
  cpl = seq(min(sim_base$cpl), max(sim_base$cpl), length.out = 100),
  cr = mean(sim_base$cr),
  ratio_factor = factor(c(0.10, 0.15, 0.05), levels = c(0.10, 0.15, 0.05))
)

cpl_grid$pred_roi <- predict(ols_model_int, newdata = cpl_grid)
cpl_grid$ratio_label <- factor(
  ifelse(cpl_grid$ratio_factor == 0.05, "5%",
         ifelse(cpl_grid$ratio_factor == 0.10, "10%", "15%")),
  levels = c("5%", "10%", "15%")
)

fig5 <- ggplot(cpl_grid, aes(x = cpl, y = pred_roi, linetype = ratio_label)) +
  geom_line(color = "black", linewidth = 1) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(
    title = "Predicted ROI by Cost per Lead Across Allocation Scenarios",
    x = "Cost per Lead (CPL)",
    y = "Predicted ROI",
    linetype = "Allocation Ratio"
  ) +
  theme_ess_bw()

print(fig5)

# -----------------------------
# 12. Optimization
# -----------------------------
risk_adjusted_objective <- function(ratio) {
  tmp <- simulate_ess(
    ratio = ratio,
    n = 5000,
    total_budget = total_budget,
    discount_rate = discount_rate,
    tuition_per_enrollment = tuition_per_enrollment,
    years = years,
    fixed_cost_mean = fixed_cost_mean,
    fixed_cost_sd = fixed_cost_sd,
    variable_cost_mean = variable_cost_mean,
    variable_cost_sd = variable_cost_sd,
    cpl_min = cpl_min,
    cpl_max = cpl_max,
    cr_min = cr_min,
    cr_max = cr_max,
    scenario_name = "Optimization"
  )

  mu <- mean(tmp$roi)
  sig <- sd(tmp$roi)

  -(mu / sig)
}

opt_res <- optimize(
  f = risk_adjusted_objective,
  interval = c(0.05, 0.20)
)

optimal_ratio <- opt_res$minimum
print(optimal_ratio)

sim_opt <- simulate_ess(
  ratio = optimal_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
  years = years,
  fixed_cost_mean = fixed_cost_mean,
  fixed_cost_sd = fixed_cost_sd,
  variable_cost_mean = variable_cost_mean,
  variable_cost_sd = variable_cost_sd,
  cpl_min = cpl_min,
  cpl_max = cpl_max,
  cr_min = cr_min,
  cr_max = cr_max,
  scenario_name = "Optimal"
)

opt_summary <- summarise_simulation(sim_opt) %>%
  mutate(optimal_ratio = optimal_ratio) %>%
  relocate(optimal_ratio)

print(opt_summary)

# Optimization curve
ratio_grid <- seq(0.05, 0.20, length.out = 60)

opt_curve <- data.frame(
  ratio = ratio_grid,
  objective = sapply(ratio_grid, function(r) {
    tmp <- simulate_ess(
      ratio = r,
      n = 3000,
      total_budget = total_budget,
      discount_rate = discount_rate,
      tuition_per_enrollment = tuition_per_enrollment,
      years = years,
      fixed_cost_mean = fixed_cost_mean,
      fixed_cost_sd = fixed_cost_sd,
      variable_cost_mean = variable_cost_mean,
      variable_cost_sd = variable_cost_sd,
      cpl_min = cpl_min,
      cpl_max = cpl_max,
      cr_min = cr_min,
      cr_max = cr_max,
      scenario_name = "OptimizationCurve"
    )

    mean(tmp$roi) / sd(tmp$roi)
  })
)

fig6 <- ggplot(opt_curve, aes(x = ratio, y = objective)) +
  geom_line(color = "black", linewidth = 1) +
  geom_vline(xintercept = optimal_ratio, linetype = "dashed", color = "black") +
  annotate(
    "text",
    x = optimal_ratio - 0.01,
    y = max(opt_curve$objective) * 0.98,
    label = paste0("Optimal ≈ ", round(optimal_ratio * 100, 1), "%"),
    hjust = 1,
    size = 4
  ) +
  labs(
    title = "Risk-Adjusted Objective Across Marketing Allocation Ratios",
    x = "Marketing Allocation Ratio",
    y = "Risk-Adjusted Objective (Mean ROI / SD of ROI)"
  ) +
  theme_ess_bw()

print(fig6)

# Risk-return frontier
frontier_grid <- seq(0.05, 0.20, length.out = 40)

frontier_data <- data.frame(
  ratio = frontier_grid,
  mean_roi = NA_real_,
  sd_roi = NA_real_
)

for (i in seq_along(frontier_grid)) {
  tmp <- simulate_ess(
    ratio = frontier_grid[i],
    n = 3000,
    total_budget = total_budget,
    discount_rate = discount_rate,
    tuition_per_enrollment = tuition_per_enrollment,
    years = years,
    fixed_cost_mean = fixed_cost_mean,
    fixed_cost_sd = fixed_cost_sd,
    variable_cost_mean = variable_cost_mean,
    variable_cost_sd = variable_cost_sd,
    cpl_min = cpl_min,
    cpl_max = cpl_max,
    cr_min = cr_min,
    cr_max = cr_max,
    scenario_name = "Frontier"
  )

  frontier_data$mean_roi[i] <- mean(tmp$roi)
  frontier_data$sd_roi[i] <- sd(tmp$roi)
}

frontier_opt <- frontier_data[which.min(abs(frontier_data$ratio - optimal_ratio)), ]

fig7 <- ggplot(frontier_data, aes(x = sd_roi, y = mean_roi)) +
  geom_path(color = "black", linewidth = 1) +
  geom_point(color = "black", size = 2) +
  geom_point(
    data = frontier_opt,
    aes(x = sd_roi, y = mean_roi),
    color = "black",
    fill = "white",
    shape = 21,
    size = 3,
    stroke = 1
  ) +
  annotate(
    "text",
    x = frontier_opt$sd_roi,
    y = frontier_opt$mean_roi,
    label = paste0("  Optimal ≈ ", round(optimal_ratio * 100, 1), "%"),
    hjust = 0,
    vjust = -0.5,
    size = 4
  ) +
  labs(
    title = "Risk-Return Frontier for Marketing Allocation Ratios",
    x = "Risk (SD of ROI)",
    y = "Expected Return (Mean ROI)"
  ) +
  theme_ess_bw()

print(fig7)

# -----------------------------
# 13. Policy scenarios
# -----------------------------
policy_ratio <- 0.10

scen_base <- simulate_ess(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
  years = years,
  fixed_cost_mean = fixed_cost_mean,
  fixed_cost_sd = fixed_cost_sd,
  variable_cost_mean = variable_cost_mean,
  variable_cost_sd = variable_cost_sd,
  cpl_min = cpl_min,
  cpl_max = cpl_max,
  cr_min = cr_min,
  cr_max = cr_max,
  scenario_name = "Base"
)

scen_funding_cut <- simulate_ess(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
  years = years,
  fixed_cost_mean = fixed_cost_mean,
  fixed_cost_sd = fixed_cost_sd,
  variable_cost_mean = variable_cost_mean,
  variable_cost_sd = variable_cost_sd,
  cpl_min = cpl_min,
  cpl_max = cpl_max,
  cr_min = cr_min,
  cr_max = cr_max,
  scenario_name = "Funding Cut",
  budget_multiplier = 0.90
)

scen_ad_cost_rise <- simulate_ess(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
  years = years,
  fixed_cost_mean = fixed_cost_mean,
  fixed_cost_sd = fixed_cost_sd,
  variable_cost_mean = variable_cost_mean,
  variable_cost_sd = variable_cost_sd,
  cpl_min = cpl_min,
  cpl_max = cpl_max,
  cr_min = cr_min,
  cr_max = cr_max,
  scenario_name = "Ad Cost Rise",
  cpl_multiplier = 1.30
)

scen_low_demand <- simulate_ess(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
  years = years,
  fixed_cost_mean = fixed_cost_mean,
  fixed_cost_sd = fixed_cost_sd,
  variable_cost_mean = variable_cost_mean,
  variable_cost_sd = variable_cost_sd,
  cpl_min = cpl_min,
  cpl_max = cpl_max,
  cr_min = cr_min,
  cr_max = cr_max,
  scenario_name = "Low Demand",
  cr_multiplier = 0.90
)

policy_data <- bind_rows(
  scen_base,
  scen_funding_cut,
  scen_ad_cost_rise,
  scen_low_demand
)

table3 <- policy_data %>%
  group_by(scenario) %>%
  summarise(
    mean_roi = mean(roi),
    sd_roi = sd(roi),
    p_roi_pos = mean(roi > 0),
    var_5 = as.numeric(quantile(roi, 0.05)),
    cvar_5 = mean(roi[roi <= quantile(roi, 0.05)]),
    .groups = "drop"
  )

print(table3)

table3$scenario <- factor(
  table3$scenario,
  levels = c("Ad Cost Rise", "Low Demand", "Funding Cut", "Base")
)

fig8 <- ggplot(table3, aes(x = scenario, y = mean_roi)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_y_continuous(
    limits = c(min(table3$mean_roi) * 1.1, 0),
    expand = expansion(mult = c(0.05, 0))
  ) +
  coord_flip() +
  labs(
    title = "Policy Scenario Comparison: Mean ROI",
    x = "Scenario",
    y = "Mean ROI"
  ) +
  theme_ess_bw()

print(fig8)

# -----------------------------
# 14. Sensitivity decomposition
# -----------------------------
sens_model <- lm(roi ~ cr + cpl + ratio, data = sim_base)
anova_sens <- anova(sens_model)
anova_sens <- anova_sens %>%
  mutate(
    term = rownames(anova_sens),
    prop_ss = `Sum Sq` / sum(`Sum Sq`, na.rm = TRUE)
  ) %>%
  select(term, everything())
print(anova_sens)

anova_plot_data <- anova_sens %>%
  filter(term != "Residuals") %>%
  mutate(term = factor(term, levels = c("cr", "cpl", "ratio")))

fig9 <- ggplot(anova_plot_data, aes(x = term, y = prop_ss)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  scale_x_discrete(labels = c("cr" = "CR",        # <-- ADD THIS
                              "cpl" = "CPL",
                              "ratio" = "Allocation Ratio")) +
  labs(
    title = "Approximate Variance Decomposition of ROI",
    x = "Predictor",
    y = "Proportion of Sum of Squares"
  ) +
  theme_ess_bw()

print(fig9)

# -----------------------------
# 15. Optional robustness checks
# -----------------------------
print(summary(ols_model))
konfound(ols_model, tested = "cr")
konfound(ols_model, tested = "cpl")

# -----------------------------
# 16. Tuition scenario analysis
# -----------------------------
tuition_grid <- c(10000, 12000, 14000, 16000, 18000, 20000)
policy_ratios <- c(0.05, 0.10, 0.15, 0.20)

positive_roi_table <- expand.grid(
  tuition_per_enrollment = tuition_grid,
  ratio = policy_ratios
)

positive_roi_table$mean_roi <- NA_real_
positive_roi_table$sd_roi <- NA_real_
positive_roi_table$p_roi_pos <- NA_real_

for (i in seq_len(nrow(positive_roi_table))) {
  tmp <- simulate_ess(
    ratio = positive_roi_table$ratio[i],
    n = 10000,
    total_budget = total_budget,
    discount_rate = discount_rate,
    tuition_per_enrollment = positive_roi_table$tuition_per_enrollment[i],
    years = years,
    fixed_cost_mean = fixed_cost_mean,
    fixed_cost_sd = fixed_cost_sd,
    variable_cost_mean = variable_cost_mean,
    variable_cost_sd = variable_cost_sd,
    cpl_min = cpl_min,
    cpl_max = cpl_max,
    cr_min = cr_min,
    cr_max = cr_max,
    scenario_name = "PositiveROI_Tuition"
  )

  positive_roi_table$mean_roi[i] <- mean(tmp$roi)
  positive_roi_table$sd_roi[i] <- sd(tmp$roi)
  positive_roi_table$p_roi_pos[i] <- mean(tmp$roi > 0)
}

positive_roi_table <- positive_roi_table %>%
  mutate(
    ratio_label = factor(
      paste0(round(ratio * 100), "%"),
      levels = c("5%", "10%", "15%", "20%")
    )
  )

print(positive_roi_table)

positive_roi_summary <- positive_roi_table %>%
  mutate(
    mean_roi = round(mean_roi, 3),
    sd_roi = round(sd_roi, 3),
    p_roi_pos = round(p_roi_pos, 3)
  ) %>%
  arrange(tuition_per_enrollment, ratio)

print(positive_roi_summary)

fig10 <- ggplot(
  positive_roi_table,
  aes(x = tuition_per_enrollment, y = mean_roi, linetype = ratio_label)
) +
  geom_line(color = "black", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "twodash")) +
  labs(
    title = "Mean ROI Under Alternative Tuition Revenue Assumptions",
    x = "Tuition Revenue per Enrollment per Year ($)",
    y = "Mean ROI",
    linetype = "Allocation Ratio"
  ) +
  theme_ess_bw()

print(fig10)

