# ============================================================
# Marketing-Finance Simulation (MFS) Framework
# Monte Carlo simulation for higher education marketing analysis
# Author: Subir Hait
# ============================================================

rm(list = ls())

# -----------------------------
# 0. Packages
# -----------------------------
needed <- c("dplyr", "ggplot2", "quantreg", "purrr", "tidyr", "broom")
to_install <- needed[!needed %in% installed.packages()[, "Package"]]
if(length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(ggplot2)
library(quantreg)
library(purrr)
library(tidyr)
library(broom)

set.seed(12345)

# -----------------------------
# 1. Global settings
# -----------------------------
n_sim        <- 10000
total_budget <- 50000
discount_rate <- 0.05

# marketing allocation ratios
ratios <- c(0.05, 0.10, 0.15)

# baseline tuition revenue per enrolled student
tuition_per_enrollment <- 10000

# fixed cost and variable cost assumptions
fixed_cost_mean <- 20000
fixed_cost_sd   <- 1000

variable_cost_mean <- 5000
variable_cost_sd   <- 800

# uncertain marketing inputs
# You can adjust these ranges later if needed
cpl_min <- 150   # cost per lead lower bound
cpl_max <- 500   # cost per lead upper bound

cr_min  <- 0.02  # conversion rate lower bound
cr_max  <- 0.07  # conversion rate upper bound

# -----------------------------
# 2. Core simulation function
# -----------------------------
simulate_mfs <- function(ratio,
                         n = 10000,
                         total_budget = 50000,
                         discount_rate = 0.05,
                         tuition_per_enrollment = 10000,
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

  effective_budget <- total_budget * budget_multiplier
  marketing_budget <- effective_budget * ratio

  cpl <- runif(n, min = cpl_min, max = cpl_max) * cpl_multiplier
  cr  <- runif(n, min = cr_min, max = cr_max) * cr_multiplier

  fixed_cost <- rnorm(n, mean = fixed_cost_mean, sd = fixed_cost_sd)
  variable_cost <- rnorm(n, mean = variable_cost_mean, sd = variable_cost_sd)

  # avoid impossible negative costs
  fixed_cost <- pmax(fixed_cost, 0)
  variable_cost <- pmax(variable_cost, 0)

  leads <- marketing_budget / cpl
  enrollments <- leads * cr

  revenue <- enrollments * tuition_per_enrollment
  total_cost <- fixed_cost + marketing_budget + variable_cost

  profit <- revenue - total_cost
  roi <- profit / total_cost
  npv <- profit / (1 + discount_rate)

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
    revenue = revenue,
    total_cost = total_cost,
    profit = profit,
    roi = roi,
    npv = npv
  )
}

# -----------------------------
# 3. Baseline simulation
# -----------------------------
sim_base <- map_dfr(ratios, ~simulate_mfs(
  ratio = .x,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
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

# -----------------------------
# 4. Summary table by ratio
# -----------------------------
table1 <- sim_base %>%
  group_by(ratio) %>%
  summarise(
    mean_roi = mean(roi),
    sd_roi = sd(roi),
    mean_npv = mean(npv),
    mean_enrollments = mean(enrollments),
    mean_leads = mean(leads),
    .groups = "drop"
  )

print(table1)

# -----------------------------
# 5. Visualizations
# -----------------------------
# ROI density
p1 <- ggplot(sim_base, aes(x = roi, linetype = factor(ratio))) +
  geom_density() +
  labs(
    title = "ROI Distribution by Marketing Allocation Ratio",
    x = "Return on Investment (ROI)",
    y = "Density",
    linetype = "Ratio"
  ) +
  theme_minimal()

print(p1)

# CR vs ROI
p2 <- ggplot(sample_n(sim_base, min(3000, nrow(sim_base))), aes(x = cr, y = roi)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Conversion Rate and ROI",
    x = "Conversion Rate",
    y = "ROI"
  ) +
  theme_minimal()

print(p2)

# CPL vs ROI
p3 <- ggplot(sample_n(sim_base, min(3000, nrow(sim_base))), aes(x = cpl, y = roi)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Cost per Lead and ROI",
    x = "Cost per Lead (CPL)",
    y = "ROI"
  ) +
  theme_minimal()

print(p3)

# -----------------------------
# 6. Econometric models
# -----------------------------
sim_base <- sim_base %>%
  mutate(
    ratio_factor = factor(ratio, levels = c(0.10, 0.15, 0.05))
  )

ols_model <- lm(roi ~ cr + cpl + ratio_factor, data = sim_base)
summary(ols_model)

qr_model <- rq(roi ~ cr + cpl + ratio_factor, tau = 0.25, data = sim_base)
summary(qr_model)

ols_table <- broom::tidy(ols_model)
print(ols_table)

qr_table <- broom::tidy(qr_model)
print(qr_table)

# -----------------------------
# 7. Optimization
# -----------------------------
# Objective: maximize mean(ROI) / sd(ROI)
risk_adjusted_objective <- function(ratio) {
  tmp <- simulate_mfs(
    ratio = ratio,
    n = 5000,
    total_budget = total_budget,
    discount_rate = discount_rate,
    tuition_per_enrollment = tuition_per_enrollment,
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

  # maximize mu/sig -> minimize negative
  -(mu / sig)
}

opt_res <- optimize(
  f = risk_adjusted_objective,
  interval = c(0.05, 0.20)
)

optimal_ratio <- opt_res$minimum
optimal_ratio

sim_opt <- simulate_mfs(
  ratio = optimal_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
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

opt_summary <- sim_opt %>%
  summarise(
    optimal_ratio = optimal_ratio,
    mean_roi = mean(roi),
    sd_roi = sd(roi),
    mean_npv = mean(npv),
    mean_enrollments = mean(enrollments),
    mean_leads = mean(leads),
    p_positive_roi = mean(roi > 0),
    p_positive_npv = mean(npv > 0),
    var_5 = quantile(roi, 0.05),
    cvar_5 = mean(roi[roi <= quantile(roi, 0.05)])
  )

print(opt_summary)

# -----------------------------
# 8. Policy scenarios
# -----------------------------
# Base at chosen ratio (say 10% or optimal ratio)
policy_ratio <- 0.10

scen_base <- simulate_mfs(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
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

scen_funding_cut <- simulate_mfs(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
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

scen_ad_cost_rise <- simulate_mfs(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
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

scen_low_demand <- simulate_mfs(
  ratio = policy_ratio,
  n = n_sim,
  total_budget = total_budget,
  discount_rate = discount_rate,
  tuition_per_enrollment = tuition_per_enrollment,
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
    var_5 = quantile(roi, 0.05),
    cvar_5 = mean(roi[roi <= quantile(roi, 0.05)]),
    .groups = "drop"
  )

print(table3)

# bar chart for mean ROI
p4 <- ggplot(table3, aes(x = scenario, y = mean_roi)) +
  geom_col() +
  labs(
    title = "Policy Scenario Comparison: Mean ROI",
    x = "Scenario",
    y = "Mean ROI"
  ) +
  theme_minimal()

print(p4)

# -----------------------------
# 9. Simple sensitivity decomposition
# -----------------------------
# Approximation: variance explained via linear model ANOVA
sens_model <- lm(roi ~ cr + cpl + ratio, data = sim_base)
anova_sens <- anova(sens_model)

anova_sens <- anova_sens %>%
  mutate(
    term = rownames(anova_sens),
    prop_ss = `Sum Sq` / sum(`Sum Sq`, na.rm = TRUE)
  )

print(anova_sens)

# -----------------------------
# install.packages("konfound")
library(konfound)

summary(ols_model)

konfound(ols_model, tested = "cr")
konfound(ols_model, tested = "cpl")



ols_model_int <- lm(roi ~ cr * ratio_factor + cpl * ratio_factor, data = sim_base)
summary(ols_model_int)

qr_model_int <- rq(roi ~ cr * ratio_factor + cpl * ratio_factor,
                   tau = 0.25, data = sim_base)
summary(qr_model_int)

broom::tidy(ols_model_int)
broom::tidy(qr_model_int)
# Figure 5: Predicted ROI by CR across allocation scenarios
cr_grid <- expand.grid(
  cr = seq(min(sim_base$cr), max(sim_base$cr), length.out = 100),
  cpl = mean(sim_base$cpl),
  ratio_factor = factor(c(0.10, 0.15, 0.05), levels = c(0.10, 0.15, 0.05))
)

cr_grid$pred_roi <- predict(ols_model_int, newdata = cr_grid)

fig5 <- ggplot(cr_grid, aes(x = cr, y = pred_roi, linetype = ratio_factor)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Predicted ROI by Conversion Rate Across Allocation Scenarios",
    x = "Conversion Rate",
    y = "Predicted ROI",
    linetype = "Allocation Ratio"
  ) +
  theme_minimal()

print(fig5)
# Figure 6: Predicted ROI by CPL across allocation scenarios
cpl_grid <- expand.grid(
  cpl = seq(min(sim_base$cpl), max(sim_base$cpl), length.out = 100),
  cr = mean(sim_base$cr),
  ratio_factor = factor(c(0.10, 0.15, 0.05), levels = c(0.10, 0.15, 0.05))
)

cpl_grid$pred_roi <- predict(ols_model_int, newdata = cpl_grid)

fig6 <- ggplot(cpl_grid, aes(x = cpl, y = pred_roi, linetype = ratio_factor)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Predicted ROI by Cost per Lead Across Allocation Scenarios",
    x = "Cost per Lead (CPL)",
    y = "Predicted ROI",
    linetype = "Allocation Ratio"
  ) +
  theme_minimal()

print(fig6)


# -----------------------------
# -----------------------------
# 7b. optimization curve
# -----------------------------
ratio_grid <- seq(0.05, 0.20, length.out = 60)

opt_curve <- data.frame(
  ratio = ratio_grid,
  objective = sapply(ratio_grid, function(r) {
    tmp <- simulate_mfs(
      ratio = r,
      n = 3000,   # increase to 5000 if you want a smoother curve
      total_budget = total_budget,
      discount_rate = discount_rate,
      tuition_per_enrollment = tuition_per_enrollment,
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

fig7 <- ggplot(opt_curve, aes(x = ratio, y = objective)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = optimal_ratio, linetype = "dashed") +
  annotate(
    "text",
    x = optimal_ratio - 0.01,
    y = max(opt_curve$objective) - 0.1,
    label = paste0("Optimal ≈ ", round(optimal_ratio * 100, 1), "%"),
    hjust = 1,
    size = 4
  ) +
  labs(
    title = "Risk-Adjusted Objective Across Marketing Allocation Ratios",
    x = "Marketing Allocation Ratio",
    y = "Risk-Adjusted Objective (Mean ROI / SD[ROI])"
  ) +
  theme_minimal()

print(fig7)




fig7 <- ggplot(opt_curve, aes(x = ratio, y = objective)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = optimal_ratio, linetype = "dashed") +
  annotate(
    "text",
    x = optimal_ratio - 0.015,
    y = max(opt_curve$objective) - 0.1,
    label = paste0("Optimal ≈ ", round(optimal_ratio * 100, 1), "%"),
    hjust = 1,
    size = 4
  ) +
  labs(
    title = "Risk-Adjusted Objective Across Marketing Allocation Ratios",
    x = "Marketing Allocation Ratio",
    y = "Risk-Adjusted Objective (Mean ROI / SD of ROI)"
  ) +
  theme_minimal()

print(fig7)




# -----------------------------
# 7c. Risk-return frontier
# -----------------------------
frontier_grid <- seq(0.05, 0.20, length.out = 40)

frontier_data <- data.frame(
  ratio = frontier_grid,
  mean_roi = NA_real_,
  sd_roi = NA_real_
)

for (i in seq_along(frontier_grid)) {
  tmp <- simulate_mfs(
    ratio = frontier_grid[i],
    n = 3000,   # increase to 5000 for smoother estimates if desired
    total_budget = total_budget,
    discount_rate = discount_rate,
    tuition_per_enrollment = tuition_per_enrollment,
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
  frontier_data$sd_roi[i]   <- sd(tmp$roi)
}

# locate the point closest to the optimized ratio
frontier_opt <- frontier_data[which.min(abs(frontier_data$ratio - optimal_ratio)), ]

fig8 <- ggplot(frontier_data, aes(x = sd_roi, y = mean_roi)) +
  geom_path(linewidth = 1) +
  geom_point(size = 2) +
  geom_point(
    data = frontier_opt,
    aes(x = sd_roi, y = mean_roi),
    size = 3
  ) +
  annotate(
    "text",
    x = frontier_opt$sd_roi - 0.005,
    y = frontier_opt$mean_roi,
    label = paste0("Optimal ≈ ", round(optimal_ratio * 100, 1), "%"),
    hjust = 1,
    size = 4
  ) +
  labs(
    title = "Risk-Return Frontier for Marketing Allocation Ratios",
    x = "Risk (SD of ROI)",
    y = "Expected Return (Mean ROI)"
  ) +
  theme_minimal()

print(fig8)

# -----------------------------
# 10. Positive-ROI scenario analysis
#    Vary tuition revenue per enrollment
# -----------------------------
tuition_grid <- c(10000, 12000, 14000, 16000, 18000, 20000)
policy_ratios <- c(0.05, 0.10, 0.15, 0.20)

positive_roi_table <- expand.grid(
  tuition_per_enrollment = tuition_grid,
  ratio = policy_ratios
)

positive_roi_table$mean_roi <- NA_real_
positive_roi_table$sd_roi   <- NA_real_
positive_roi_table$p_roi_pos <- NA_real_

for (i in seq_len(nrow(positive_roi_table))) {
  tmp <- simulate_mfs(
    ratio = positive_roi_table$ratio[i],
    n = 10000,
    total_budget = total_budget,
    discount_rate = discount_rate,
    tuition_per_enrollment = positive_roi_table$tuition_per_enrollment[i],
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

  positive_roi_table$mean_roi[i]  <- mean(tmp$roi)
  positive_roi_table$sd_roi[i]    <- sd(tmp$roi)
  positive_roi_table$p_roi_pos[i] <- mean(tmp$roi > 0)
}

positive_roi_table
positive_roi_table %>%
  mutate(
    ratio = paste0(round(ratio * 100), "%"),
    mean_roi = round(mean_roi, 3),
    sd_roi = round(sd_roi, 3),
    p_roi_pos = round(p_roi_pos, 3)
  ) %>%
  arrange(tuition_per_enrollment, ratio)
# Figure: Mean ROI by tuition revenue and allocation ratio
fig_positive_roi <- ggplot(
  positive_roi_table,
  aes(x = tuition_per_enrollment, y = mean_roi, linetype = factor(ratio))
) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    title = "Mean ROI Under Alternative Tuition Revenue Assumptions",
    x = "Tuition Revenue per Enrollment ($)",
    y = "Mean ROI",
    linetype = "Allocation Ratio"
  ) +
  theme_minimal()

print(fig_positive_roi)



years <- 4
tuition_stream <- tuition_per_enrollment *
  sum(1 / (1 + discount_rate)^(0:(years-1)))

revenue <- enrollments * tuition_stream
taus <- c(0.25, 0.50, 0.75)

qr_models <- lapply(taus, function(t) {
  rq(roi ~ cr + cpl + ratio_factor, tau = t, data = sim_base)
})

qr_tables <- lapply(qr_models, broom::tidy)

qr_results <- dplyr::bind_rows(qr_tables, .id = "tau")
qr_results
