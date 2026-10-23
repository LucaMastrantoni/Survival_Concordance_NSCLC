# Shiny
library(shiny)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(splines)
library(promises)
library(future)
library(transport)
library(survRM2)
library(gtsummary)
library(survival)
library(survminer)
library(broom.helpers)
library(broom)

retrodesign_HR <- function(HR_assumed, se, alpha = 0.05, alpha_sided = 2, n.sims = 20000) {
  
  D <- log(HR_assumed)
  z <- qnorm(1 - alpha / alpha_sided)
  
  # one-sided test: only the beneficial direction (HR < 1) can be significant
  power_hi <- if (alpha_sided == 1) 0 else pnorm(-z + D / se)
  power_lo <- pnorm(-z - D / se)
  power    <- power_hi + power_lo

  typeS <- power_lo / power

  sims <- rnorm(n.sims, mean = D, sd = se)
  sig  <- if (alpha_sided == 1) sims < -z * se else abs(sims) > z * se
  typeM <- mean(abs(sims[sig]) / abs(D))
  
  data.frame(
    Power = power,
    Type_S_error = typeS,
    Type_M_exaggeration = typeM
  )
}

ui <- navbarPage(
  
  title = "Survival Tools for RCT Concordance",
  
  ## --- TAB 0: Read me ---
  tabPanel(
    "Read me",
    
    fluidPage(
      
      tags$h2("Survival Tools for Randomized Controlled Trials Concordance"),
      
      tags$h3("Overview"),
      p("The app is intended for researchers working with randomized controlled trials, and individual patient data (IPD). It allows you to:"
      ),
      tags$ul(
        tags$li("Evaluate power, Type S error, and Type M exaggeration under alternative design assumptions."),
        tags$li("Visualize KM curves and fit Cox proportional hazards models."),
        tags$li("Compare Bayesian posterior hazard ratio distributions against the original design assumptions, and calculate Wasserstein-1 distance."),
        tags$li("Quantify discrepancies between hazard ratios and restricted mean survival time (RMST)-based treatment effects.")
      ),
      tags$hr(),
      
      tags$h3("Tab description"),
      p("Tabs are intended to be run in order, as information is shared." ),
      
      tags$h4("1. Retrodesign (Single Trial)"),
      p("Implementation of the retrodesign analysis, as proposed by Gelman et al."),
      tags$ul(
        tags$li("Statistical power as a function of the assumed hazard ratio."),
        tags$li("Type S error (probability of sign error)."),
        tags$li("Type M error (exaggeration factor among statistically significant results).")
      ),
      
      tags$h4("2. Survival analysis (IPD)"),
      p("Standard frequentist survival analysis using individual patient data from recontructed KM-curves. Uploaded data must include time, event status, and treatment indicator."
      ),
      tags$ul(
        tags$li("Kaplan-Meier curves with confidence intervals and risk tables."),
        tags$li("Median survival times with 95% confidence intervals."),
        tags$li("Cox proportional hazards regression."),
        tags$li("Assessment of the proportional hazards assumption using Schoenfeld residuals."),
      ),
      
      tags$h4("3. Bayesian Cox (IPD)"),
      p("Comparison of Bayesian posterior distributions of the treatment log-hazard ratio against the original design assumptions."),
      tags$ul(
        tags$li("Upload posterior draws from a Bayesian Cox model (log-HR scale). The app does not directly calculate posterior draws, We provided an R script to do it locally ('bayesian.R'), requiring the same IPD input data which produces the correct 'posterior_b_treat.rds' file."),
        tags$li("Wasserstein distance between posterior and design distributions."),
        tags$li("Density plots for visual comparison.")
      ),
      
      tags$h4("4. RMST (IPD)"),
      p("Restricted mean survival time (RMST) analysis."),
      tags$ul(
        tags$li("RMST estimated directly from Kaplan–Meier curves."),
        tags$li("RMST difference and RMST ratio evaluated over increasing τ."),
        tags$li("Quantification of mismatch between RMST under an exponential proportional hazards assumption and the observed RMST ratio from KM curves."),
        
      ),
      
      tags$hr(),
      
      tags$h4("Input data requirements"),
      tags$ul(
        tags$li("IPD files must contain columns named: 'time', 'status' (0 = censored, 1 = event), and 'treat' (0 = control, 1 = experimental)."),
        tags$li("Bayesian posterior files must be .rds objects containing a column named b_treat. Please use the R script 'bayesian.R'."),
        tags$li("Example files are provided on GitHub, including IPD data (ALEX_KM.xlsx), and posterior draws ('posterior_b_treat.rds'). The default values in the app are referred to the ALEX trial ('Peters S, Camidge DR, Shaw AT, et al. Alectinib versus Crizotinib in Untreated ALK -Positive Non–Small-Cell Lung Cancer. New England Journal of Medicine 2017; 377: 829–38.')."),
        tags$li("This app does not include a method to reconstruct IPD from published KM curves. We refer to https://doi.org/10.1186/1471-2288-12-9' and the r package IPDfromKM (https://CRAN.R-project.org/package=IPDfromKM).")
        ),
      
      tags$hr(),
      tags$h4("Contact"),
      tags$ul(
        tags$li("Correspondence: luca.mastrantoni01@icatt.it"),
        tags$li("GitHub: https://github.com/LucaMastrantoni")),
      
      tags$hr(),
      tags$h4("Affiliation"),
      tags$ul(
        tags$li("Università Cattolica del Sacro Cuore, Rome, Italy."),
        tags$li("Comprehensive Cancer Center, Fondazione Policlinico Universitario Agostino Gemelli, IRCCS Rome, Italy.")
      ),
      div(style = "margin-top: 20px; text-align: left;", tags$img(src = "logo.PNG", width = "600px")),
      
      tags$hr(),
      helpText("Disclaimer. This tool is intended for methodological exploration. Results should be interpreted in the context of trial design assumptions."
      )
    )
  ),
        
  
  ## --- TAB 1: Retrodesign ---
  tabPanel(
    "1. Retrodesign (Single Trial)",
    sidebarLayout(
      sidebarPanel(
        tags$h3("Inputs"),
        helpText("Insert the following inputs for retrodesign analysis."),
        
        numericInput("HR_design", "Design HR", 0.65, min = 0.35, max = 0.80), 
        helpText("HR used for trial design, or HR for which Type S and Type M errors should be calculated."),
        
        numericInput("HR_obs", "Observed HR", 0.47, min = 0.01, max = 3),
        helpText("HR observed in the trial, as reported by trialists or recalculated from IPD."),
        numericInput("HR_obs_low", "Observed HR (Lower CI)", 0.34, min = 0.1, max = 4),
        helpText("Lower bound of confidence interval for observed HR."),
        numericInput("HR_obs_up", "Observed HR (Upper CI)", 0.65, min = 0.1, max = 4),
        helpText("Upper bound of confidence interval for observed HR."),
        
        numericInput("power_design", "Planned Power", 0.8, min = 0.8, max = 0.99),
        helpText("Power used for original sample size calculation."),
        
        numericInput("alpha", "Alpha", 0.05),
        helpText("Alpha used for original sample size calculation."),
        selectInput("sided", "Test", choices = c("Two-sided" = 2, "One-sided" = 1)),
        helpText("Alpha sides used for original sample size calculation."),
        
        #numericInput("nsims", "Monte Carlo simulations", 10000),
        
        actionButton("run_retro", "Run")
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.run_retro > 0",
          
          tags$h3("Retrodesign summary"),
          tableOutput("summary"),
          tags$p("Type S error and Type M error for the HR used for trial design (or other HR selected)."),
          
          tags$hr(),
          tags$h3("Power curve"),
          tags$p("Power for different values of the HR. Black dot is observed HR, with 95% CI shaded in gray. Red dot is the HR used for trial design (or other HR selected)."),
          plotOutput("power_curve", height = "400px"),
          
          tags$hr(),
          tags$h3("Exaggeration ratio curve"),
          tags$p("Exaggeration ratios for different values of the HR. Black dot is observed HR, with 95% CI shaded in gray. Red dot is the HR used for trial design (or other HR selected)."),
          plotOutput("exagg_curve", height = "400px")
        )
      )
    )
  ),
  
  ## --- TAB 2: Survival (KM + Cox) ---
  tabPanel(
    "2. Survival analysis (IPD)",
    sidebarLayout(
      sidebarPanel(
        tags$h3("Inputs"),
        helpText("Upload individual patient data (IPD) for evaluation of Kaplan-Meier reconstruction."),
        helpText("Required column are 'time', 'status' and 'treat'. 'treat' column should be numeric, coded as 0 = 'control arm' and 1 = 'experimental arm'. Follow the ALEX_KM.xlsx example for any doubts."),
        
        fileInput(
          "km_file_surv",
          "Upload IPD (.xlsx)",
          accept = c(".xlsx")
        ),

        actionButton("run_surv", "Run survival analysis")
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.run_surv > 0",
          
          tags$h3("Kaplan-Meier curves"),
          plotOutput("km_curves", height = "1000px"),
          tags$p("Summary of KM curves from reconstructed IPD."),
          
          tags$hr(),
          tags$h3("Cox proportional hazards model"),
          tableOutput("cph_table"),
          tags$p("Summary of Cox regression results."),
          
          tags$hr(),
          tags$h3("Median survival"),
          tableOutput("median_table"),
          tags$p("Summary of median survival times."),
          
          tags$hr(),
          tags$h3("Proportional hazards"),
          plotOutput("cph", height = "500px"),
          tags$p("Results of proportional hazard test."),
          )
        )
      )
    ),
  
  ## --- TAB 3: Bayesian Cox ---
  tabPanel(
    "3. Bayesian Cox (IPD)",
    sidebarLayout(
      sidebarPanel(
        tags$h3("Inputs"),
        helpText("Upload Bayesian analysis results for Wasserstein-1 distance calculation."),
        helpText("Please run Retrodesign first! HR design, alpha and power are taken from the retodesign panel."),
        helpText("An .rds file inlcuding posterior draws is required. Please use the provided 'bayesian.R' script to generate the file. IPD in the correct format are required. An example is provided as 'posterior_b_treat.rds'."),
        
        fileInput(
          "post",
          "Upload posterior draws (.rds)",
          accept = c(".rds")
        ),
        actionButton("run_bayes", "Show Distributions")
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.run_bayes > 0",
          
          tags$h3("Bayesian summary"),
          tableOutput("bayes_summary"),
          tags$p("Results of Bayesian analysis. Hazard Ratio (50%) and credible interval."),
          tags$hr(),
          tableOutput("bayes_summary_p"),
          tags$p("p > 0 = Probability that the observed treatment effect is greater than 0 (the experimental treatment is detrimental)."),
          tags$p("p > HR design = Probability that the observed treatment effect is greater than original HR used for sample size calculations."),
          
          tags$hr(),
          tags$h3("Wasserstein-1"),
          tags$p("Wasserstein-1 distance, indicating the degree of dissimilarity between two probability distributions by quantifying the minimal cost required to transform one distribution into the other."),
          tableOutput("bayes_compare_metrics"),
          
          tags$hr(),
          tags$h3("Expected to observed plot"),
          plotOutput("bayes_vs_design", height = "500px"),
          tags$p("Densities of the HRs for the 'Planned (Design assumption)' and the 'Observed (Bayesian Posterior)' distributions. The solid line at 1 highlights the null effect, and the dashed line the design HR. X-axis is on log(HR) scale, and labels on HR scale."),
        )
      )
    )
  ),
  
  ## --- TAB 4: RMST ---
  tabPanel(
    "4. Restricted Mean Survival Time (IPD)",
    sidebarLayout(
      sidebarPanel(
        tags$h3("Inputs"),
        helpText("Upload individual patient data (IPD) for RMST analysis."),
        
        numericInput("HR_single", "Observed HR", 0.47, min = 0.1, max = 4),
        helpText("HR observed in the trial, as reported by trialists or recalculated from IPD."),
        
        numericInput("median_C", "Median survival in control arm (months)", 10.9),
        helpText("Median survival in the control arm, as reported by trialists or recalculated from IPD."),
        
        numericInput("time_by", "τ interval", 3),
        helpText("Minimum interval time between τ values. Maximum τ is estimated automatically from follow-up time."),
        
        fileInput(
          "km_file_rmst",
          "Upload IPD (.xlsx)",
          accept = c(".xlsx")
        ),
        helpText("Required column are 'time', 'status' and 'treat'. 'treat' column should be numeric, coded as 0 = control arm and 1 = experimental arm. Follow the ALEX_KM.xlsx example for any doubts."),
        
        actionButton("run_rmst", "Compute RMST")
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.run_rmst > 0",
        
          tags$h3("RMST plot"),
          plotOutput("rmst_plot", height = "400px"),
          tableOutput("rmst_single_summary"),
          textOutput("rmst_tau_text"),
          
          tags$hr(),
          tags$h3("RMST across tau"),
          plotOutput("rmst_tau_plot", height = "400px"),
          helpText("RMST ratio and RMST difference across incresing values of τ."),
          
          tags$hr(),
          tags$h3("RMST ratio mismatch"),
          tableOutput("rmst_mismatch"),
          helpText("Mismatch between the RMST ratio expected under a perfect exponential survival distribution and the RMST observed."),
          helpText("The mismatch indicates the extent to which the observed survival departs from the exponential and proportional hazards assumptions, potentially driven by non-constant hazards and unequal or truncated follow-up times."),
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  ## ======================
  ## Retrodesign panel
  ## ======================
  
  res <- eventReactive(input$run_retro, {
    
    effects <- seq(0.1, 0.95, by = 0.01)
    # SE from the reported 95% CI, independent of the design alpha
    z_ci <- qnorm(0.975)
    se <-(log(input$HR_obs_up) - log(input$HR_obs_low)) / (2 * z_ci)
    
    df_plot <- tibble(HR_assumed = effects) %>%
      rowwise() %>%
      mutate(
        metrics = list(
          retrodesign_HR(
            HR_assumed  = HR_assumed,
            se          = se,
            alpha       = input$alpha,
            alpha_sided = as.numeric(input$sided),
            #n.sims      = input$nsims
          )
        )
      ) %>%
      unnest_wider(metrics) %>%
      ungroup()
    
    # computed directly, so HRs outside the plotting grid are not dropped
    point_metrics <- function(hr) {
      tibble(HR_assumed = hr,
             retrodesign_HR(hr, se = se, alpha = input$alpha, alpha_sided = as.numeric(input$sided)))
    }
    df_assumed_pt <- point_metrics(input$HR_design)
    df_obs_pt     <- point_metrics(input$HR_obs)
    ci_band <- data.frame(
      xmin = input$HR_obs_low,
      xmax = input$HR_obs_up,
      ymin = -Inf,
      ymax = Inf)
    
    ###
    df_assumed_pt_plot <- df_assumed_pt %>%
      rename(`Design HR` = HR_assumed, `1 - Type S error (Sign error)` = Type_S_error, `Type M error (Exaggeration ratio)` = Type_M_exaggeration)
    
    list(
      df_plot = df_plot,
      df_assumed_pt = df_assumed_pt,
      df_obs_pt = df_obs_pt,
      ci_band = ci_band,
      df_assumed_pt_plot = df_assumed_pt_plot
    )
  })
  
  output$summary <- renderTable({
    req(res())
    round(res()$df_assumed_pt_plot, 3)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$power_curve <- renderPlot({
    req(res())
    
    ggplot(res()$df_plot, aes(HR_assumed, Power)) +
      geom_rect(
        data = res()$ci_band,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "grey60", alpha = 0.3
      ) +
      geom_line(linewidth = 2) +
      geom_vline(xintercept = input$HR_obs, linetype = "dashed", linewidth = 1) +
      geom_point(
        data = res()$df_assumed_pt,
        aes(color = "Design HR"),
        shape = 21, fill = "white", size = 6, stroke = 2
      ) +
      geom_point(
        data = res()$df_obs_pt,
        aes(color = "Observed HR"),
        shape = 21, fill = "white", size = 6, stroke = 2
      ) +
      scale_color_manual(
        name = "Legend",
        values = c(
          "Design HR"   = "#B24745FF",
          "Observed HR" = "black"
        )
      ) +
      
      labs(x = "Hazard Ratio", y = "Power") +
      theme_bw(base_size = 14)
  })
  
  output$exagg_curve <- renderPlot({
    req(res())
    
    ggplot(res()$df_plot, aes(HR_assumed, Type_M_exaggeration)) +
      geom_rect(
        data = res()$ci_band,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "grey60", alpha = 0.3
      ) +
      geom_line(linewidth = 2) +
      geom_vline(
        xintercept = input$HR_obs,
        linetype = "dashed",
        linewidth = 1
      ) +
      geom_point(
        data = res()$df_assumed_pt,
        aes(color = "Design HR"),
        shape = 21, fill = "white", size = 6, stroke = 2
      ) +
      geom_point(
        data = res()$df_obs_pt,
        aes(color = "Observed HR"),
        shape = 21, fill = "white", size = 6, stroke = 2
      ) +
      scale_color_manual(
        name = "Legend",
        values = c(
          "Design HR"   = "#B24745FF",
          "Observed HR" = "black"
        )
      ) +
      labs(x = "Hazard Ratio", y = "Exaggeration (Type M)") +
      theme_bw(base_size = 14)
  })
  
  ## ======================
  ## Cox regression
  ## ======================
  
  observeEvent(input$run_surv, {
    req(input$km_file_surv)
    
    m <- readxl::read_excel(input$km_file_surv$datapath)
    
    ## ---- validation ----
    required <- c("time", "status", "treat")
    if (!all(required %in% names(m))) {
      showNotification(
        "KM file must contain columns: time, status, treat (0=control, 1=experimental)",
        type = "error",
        duration = NULL
      )
      return(NULL)
    }

    ## ---- CoX ReGrEsSiOn -----
    surv_diff <- survdiff(Surv(time, status) ~ treat, data = m)
    
    ### Median 
    ### Median survival with 95% CI
    ### Median survival with 95% CI
    km <- survfit(Surv(as.numeric(time), status) ~ treat, data = m)
    km_med <- summary(km)$table
    
    median_tbl <- tibble(
      Arm    = rownames(km_med),
      Median = km_med[, "median"],
      LCL    = km_med[, "0.95LCL"],
      UCL    = km_med[, "0.95UCL"]
    ) %>%
      mutate(
        `Median survival (95% CI)` =
          sprintf("%.2f (%.2f–%.2f)", Median, LCL, UCL)
      ) %>%
      select(Arm, `Median survival (95% CI)`) %>%
      mutate(Arm = recode(Arm, "treat=0" = "Ctr", "treat=1" = "Exp"))
    
    output$median_table <- renderTable({median_tbl}, striped = TRUE, bordered = TRUE, spacing = "s")
    
    ### Cox regression
    fit <- coxph(Surv(as.numeric(time), status) ~ treat, data = m)
    s   <- summary(fit)
    
    cox_tbl <- tibble(
      Variable = "Treatment",
      HR       = exp(s$coefficients[, "coef"]),
      LCL      = s$conf.int[, "lower .95"],
      UCL      = s$conf.int[, "upper .95"],
      p_value  = s$coefficients[, "Pr(>|z|)"]
    ) %>%
      mutate(
        `HR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", HR, LCL, UCL),
        `P-value` = ifelse(
          p_value < 0.001,
          "<0.001",
          formatC(p_value, digits = 2, format = "g")
        )
      ) %>%
      select(Variable, `HR (95% CI)`, `P-value`)
    
    output$cph_table <- renderTable({cox_tbl}, striped = TRUE, bordered = TRUE, spacing = "s")

    cz <- cox.zph(fit)
    global_stats <- cz$table["GLOBAL", ]
    pval  <- signif(global_stats["p"], 2)
    
    
    output$km_curves <- renderPlot({
      req(km)
      req(m)
      
      ggsurvplot(km, data=m, 
                 conf.int = T, surv.median.line = "hv",
                 pval = TRUE,
                 fun = "pct",
                 size = 1,
                 linetype = "strata",
                 palette = c("firebrick", "steelblue"),
                 legend.title = "Treatment",
                 legend.labs = c("Control", "Experimental"),
                 legend = "bottom",
                 break.time.by=3,
                 censor.shape = "|",
                 risk.table = c("nrisk_cumcensor"), risk.table.height = 0.15, risk.table.y.text = F, tables.y.text = F,
                 cumevents = T, cumevents.height = 0.15, cumevents.table.y.text = F,
                 ncensor.plot = T, ncensor.plot.height = 0.15)
    })
    
    output$cph <- renderPlot({
      ggcoxzph(cz)[[1]] +
        ggtitle(paste0("Schoenfeld test p=", pval))
    })
  })
    
  ## ======================
  ## Bayesian Cox panel (POSTERIOR FROM RDS)
  ## ======================
  
  bayes_res        <- reactiveVal(NULL)
  bayes_draws_long <- reactiveVal(NULL)
  bayes_metrics_rv <- reactiveVal(NULL)
  
  observeEvent(input$run_bayes, {
    
    req(input$post)   # <-- obbligatorio
    
    ## ---- load posterior draws from uploaded RDS ----
    post <- readRDS(input$post$datapath)
    
    if (!("b_treat" %in% names(post))) {
      showNotification(
        "Uploaded RDS must contain column 'b_treat' (log-HR draws).",
        type = "error",
        duration = NULL
      )
      return(NULL)
    }
    
    ## ---- freeze inputs ----
    HR_design    <- isolate(input$HR_design)
    alpha        <- isolate(input$alpha)
    sided        <- isolate(as.numeric(input$sided))
    power_design <- isolate(input$power_design)
    
    ## ---- posterior log-HR ----
    x <- post[["b_treat"]]
    
    ## ---- design distribution (log-HR scale) ----
    z_alpha <- qnorm(1 - alpha / sided)
    z_power <- qnorm(power_design)
    
    se_log_design <- abs(log(HR_design)) / (z_alpha + z_power)
    
    M <- 48000
    set.seed(1159)
    
    y <- rnorm(M, mean = log(HR_design), sd = se_log_design)
    
    ## ---- metrics ----
    P_post_gt_design <- mean(x > y)
    
    w1 <- tryCatch(
      transport::wasserstein1d(x, y),
      error = function(e) NA_real_
    )
    
    bayes_metrics_rv(
      tibble(
        #P_posterior_gt_design = P_post_gt_design,
        `Wasserstein-1 (logHR)`  = w1,
        `Wasserstein-1`          = exp(-w1) ,
        #`logHR design`     = log(HR_design),
        #`SE logHR design`  = se_log_design,
        `HR design`        = HR_design,
        `Power`          = power_design
      )
    )
    
    ## ---- store draws (LOG scale internally) ----
    bayes_draws_long(
      bind_rows(
        tibble(b_treat = x, source = "Observed (Bayesian Posterior)"),
        tibble(b_treat = y, source = "Planned (Design assumption)")
      )
    )
    
    ## ---- summary ----
    bayes_res(
      tibble(b_treat = x)
    )
  })
  
  # ---- outputs ----
  
  output$bayes_summary <- renderTable({
    req(bayes_res())
    
    bayes_res() %>%
      summarise(
        #mean    = mean(b_treat),
        `2.5%` = exp(quantile(b_treat, 0.025)),
        `25%`  = exp(quantile(b_treat, 0.25)),
        `50%`  = exp(quantile(b_treat, 0.50)),
        `75%`  = exp(quantile(b_treat, 0.75)),
        `97.5%` = exp(quantile(b_treat, 0.975)),
        .groups = "drop"
      ) %>%
      round(3)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$bayes_summary_p <- renderTable({
    req(bayes_res())
    log_HR_design <- log(isolate(input$HR_design))
    HR_design <- isolate(input$HR_design)
    
    bayes_res() %>%
      summarise(
        p1 = mean(.data$b_treat > 0),
        p2 = mean(.data$b_treat > log_HR_design),
        `HR design` = round(HR_design, 3),
        .groups = "drop"
      ) %>%
      transmute(
        `p > 0` = ifelse(p1 < 0.001, "<0.001", formatC(p1, digits = 2, format = "g")),
        `p > HR design` = ifelse(p2 < 0.001, "<0.001", formatC(p2, digits = 2, format = "g")),
        `HR design`
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$bayes_compare_metrics <- renderTable({
    req(bayes_metrics_rv())
    round(bayes_metrics_rv(), 3)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  # ---- density ----
  output$bayes_vs_design <- renderPlot({
    req(bayes_draws_long())
    log_HR_design <- log(isolate(input$HR_design))
    
    dfp <- bayes_draws_long() %>%
      mutate(HR = b_treat)  # show everything on HR scale
    
    ggplot(dfp, aes(x = HR, fill = source)) +
      geom_density(alpha = 0.6) +
      scale_fill_manual(
        values = c(
          "Observed (Bayesian Posterior)" = "#FF968A",
          "Planned (Design assumption)" = "#ABDEE6"
        )
      ) +
      geom_vline(xintercept = 0, linewidth = 0.8) +
      geom_vline(xintercept = log_HR_design, linetype = "dashed", linewidth = 0.8) +
      scale_x_continuous(
        breaks = function(x) pretty(x, n = 8),
        minor_breaks = function(x) pretty(x, n = 8),
        labels = function(x) round(exp(x), 2)
      ) +
      labs(
        x = "Hazard Ratio",
        y = "Density",
        fill = "Distribution"
      ) +
      theme_bw(base_size = 14) +
      theme(legend.position = "right")
  })
  
  ## ======================
  ## RMST (single trial, IPD)
  ## ======================
  
  observeEvent(input$run_rmst, {
    req(input$km_file_rmst)
    
    m <- readxl::read_excel(input$km_file_rmst$datapath)
    
    ## ---- validation ----
    required <- c("time", "status", "treat")
    if (!all(required %in% names(m))) {
      showNotification(
        "KM file must contain columns: time, status, treat (0=control, 1=experimental)",
        type = "error",
        duration = NULL
      )
      return(NULL)
    }

    ## ---- inputs ----
    HR  <- input$HR_single
    tau_rv <- reactiveVal(NULL)
 
    time   <- m$time
    status <- m$status
    arm    <- m$treat
    
    ## ---- RMST from KM (OBSERVED) ----
    rmst_obj <- rmst2(time = time, status = status, arm = arm)
    tau_rv(rmst_obj$tau)
    tau <- rmst_obj$tau
    
    output$rmst_tau_text <- renderText({
      req(tau_rv())
      paste0("Summary of RMST at τ = ", round(tau_rv(), 1), " months.")
    })
    
    rmst_diff <- rmst_obj$unadjusted.result["RMST (arm=1)-(arm=0)", ]
    rmst_ratio <- rmst_obj$unadjusted.result["RMST (arm=1)/(arm=0)", ]
    
    #output$rmst_plot <- renderPlot({plot(rmst_obj, xlab="Time", ylab="Probability", col = "black", col.RMST = "steelblue", col.RMTL = "firebrick") + theme_classic()})
    
    fit <- survfit(Surv(time, status) ~ arm, data = m)
    km_df <- broom::tidy(fit)
    km_rmst <- km_df %>%
      filter(time <= tau)

    output$rmst_plot <- renderPlot({
      
      ggplot() +
        geom_step(data = km_df, aes(time, estimate, color = strata), linewidth = 1) +
        geom_vline(
          xintercept = tau,
          linetype = "dashed",
          color = "grey40"
        ) +
        scale_x_continuous(
          breaks = seq(0, max(km_df$time, na.rm = TRUE), by = input$time_by)
        ) +
        scale_color_manual(
          values = c("arm=0" = "black", "arm=1" = "steelblue"),
          labels = c("Control", "Experimental")
        ) +
        labs(
          x = "Time (months)",
          y = "Survival probability",
          color = "Arm",
        ) +
        theme_classic(base_size = 14)
    })

    time.by <- input$time_by
    if (is.na(time.by) || time.by <= 0 || time.by > tau) {
      showNotification(
        paste0("τ interval must be between 0 and ", round(tau, 1), " months."),
        type = "error",
        duration = NULL
      )
      return(NULL)
    }
    # rmst2 rejects tau beyond the shorter arm's follow-up, so cap the grid at its default tau
    taus  <- seq(time.by, tau, by = time.by)
    
    metrics <- c(delta = 1, ratio = 2)
    
    rmst_grid <- purrr::imap_dfr(metrics, ~{
      row_idx <- .x
      type    <- .y
      
      purrr::map_dfr(taus, function(tau){
        
        unadj <- rmst2(
          time   = m$time,
          status = m$status,
          arm    = m$treat,
          tau    = tau
        )$unadjusted.result
        
        vals <- unadj[row_idx, ]
        
        tibble(
          type    = type,
          tau     = tau,
          est     = vals[1],
          lower95 = vals[2],
          upper95 = vals[3],
          p.value = vals[4]
        )
      })
    })
    
    output$rmst_tau_plot <- renderPlot({
      
      ggplot(rmst_grid, aes(x = tau, y = est)) +
        geom_hline(
          data    = dplyr::filter(rmst_grid, type == "delta"),
          aes(yintercept = 0),
          linetype = "dashed",
          colour   = "grey60"
        ) +
        geom_hline(
          data    = dplyr::filter(rmst_grid, type == "ratio"),
          aes(yintercept = 1),
          linetype = "dashed",
          colour   = "grey60"
        ) +
        ## estimates + CI
        geom_pointrange(
          aes(ymin = lower95, ymax = upper95),
          linewidth = 0.6
        ) +
        facet_wrap(~ type, scales = "free_y",
                   labeller = labeller(
                     type = c(
                       ratio = "RMST ratio (E / C)",
                       delta = "RMST difference (E − C)"
                     )
                   )) +
        scale_x_continuous(breaks = taus) +
        labs(
          x = "τ (months)",
          y = NULL,
          title = ""
        ) +
        theme_classic(base_size = 14) +
        theme(
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold")
        )
    })
    
    # rmst2 reports RMST difference (E - C) and ratio internally,
    # but we extract group-specific RMSTs for clarity
    rmstC_obs <- rmst_obj$RMST.arm0$rmst["Est."]
    rmstE_obs <- rmst_obj$RMST.arm1$rmst["Est."]
    
    RMST_ratio_obs <- rmstE_obs / rmstC_obs
    
    ## ---- Exponential PH expectation ----
    lambdaC <- log(2) / input$median_C
    lambdaE <- HR * lambdaC
    
    rmstC_exp <- (1 - exp(-lambdaC * tau)) / lambdaC
    rmstE_exp <- (1 - exp(-lambdaE * tau)) / lambdaE
    
    RMST_ratio_exp <- rmstE_exp / rmstC_exp
    
    ## ---- mismatch definition (coerente con HR < 1) ----
    HR_inv <- 1 / HR
    
    observed_mismatch <- RMST_ratio_obs / HR_inv
    expected_mismatch <- RMST_ratio_exp / HR_inv
    
    delta_mismatch <- observed_mismatch - expected_mismatch
    
    ## ---- output table ----
    output$rmst_single_summary <- renderTable({
      
      tibble(
        Quantity = c(
          "RMST Control (C)",
          "RMST Experimental (E)",
          "RMST difference (E − C)",
          "RMST ratio (E / C)"
        ),
        Estimate = c(
          rmst_obj$RMST.arm0$rmst["Est."],
          rmst_obj$RMST.arm1$rmst["Est."],
          rmst_diff["Est."],
          rmst_ratio["Est."]
        ),
        `Lower 95% CI` = c(
          rmst_obj$RMST.arm0$rmst["lower .95"],
          rmst_obj$RMST.arm1$rmst["lower .95"],
          rmst_diff["lower .95"],
          rmst_ratio["lower .95"]
        ),
        `Upper 95% CI` = c(
          rmst_obj$RMST.arm0$rmst["upper .95"],
          rmst_obj$RMST.arm1$rmst["upper .95"],
          rmst_diff["upper .95"],
          rmst_ratio["upper .95"]
        )
      ) %>%
        mutate(across(where(is.numeric), round, 3))
      
    }, striped = TRUE, bordered = TRUE, spacing = "s")
    output$rmst_mismatch <- renderTable({
      tibble(
        `RMST ratio observed` = RMST_ratio_obs,
        `RMST ratio expected` = RMST_ratio_exp,
        #`Observed mismatch` = observed_mismatch,
        #`Expected mismatch` = expected_mismatch,
        `Δ mismatch` = delta_mismatch
      ) %>% round(3)
    }, striped = TRUE, bordered = TRUE, spacing = "s")
  })
}  

shinyApp(ui, server)

