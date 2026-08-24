## =========================================================
## SurvSimMicro -- Microbiome Simulation Explorer
##
## Two connected workflows:
##   1. Simulate microbiome counts with SimulateMSeqU().
##   2. Generate survival outcomes from those counts (or an upload)
##      with simulate_tte().
##
## Expected layout:
##   SurvSimMicro/
##   |- app.R
##   |- R/SimulateMSeq_V2.R
##   |- R/utils.R, transform.R, predictor.R, tte.R,
##   |    censoring.R, simulate.R
##   |- params/
##      |- oral_Saliva_v2_para.RData
##      |- vag_Mid_Vagina_v1_para.RData
##      |- ...
##
## Parameter files are named {source}_{site}_v{n}_para.RData
## and hold one EstPara() result. Both .RData and .rds work.
## =========================================================

library(shiny)
library(bslib)
library(MASS)   # rnegbin, used by SimulateMSeqU
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

source("R/SimulateMSeq_V2.R")
source("R/utils.R")
source("R/transform.R")
source("R/predictor.R")
source("R/tte.R")
source("R/censoring.R")
source("R/simulate.R")

PARAM_DIR <- "params"

SOURCE_LABELS <- c(gut = "Gut", oral = "Oral",
                   vag = "Vaginal", vaginal = "Vaginal")

## ---- loading -------------------------------------------

# An EstPara() result is a list with `mu` and `ref.otu.tab`.
is_para <- function(x) {
  is.list(x) && all(c("mu", "ref.otu.tab") %in% names(x))
}

# .RData stores named objects, so pull out whichever one is a
# para object rather than assuming what it was called on save.
load_para <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    return(if (is_para(obj)) obj else NULL)
  }
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  for (nm in ls(e)) {
    obj <- get(nm, envir = e)
    if (is_para(obj)) return(obj)
    # Some files wrap the result one level deep.
    if (is.list(obj)) for (inner in obj) if (is_para(inner)) return(inner)
  }
  NULL
}

# Build the source / site / time point index from file names.
scan_params <- function() {
  if (!dir.exists(PARAM_DIR)) return(NULL)
  
  files <- list.files(PARAM_DIR, pattern = "\\.(RData|rds)$",
                      ignore.case = TRUE)
  if (!length(files)) return(NULL)
  
  stem <- sub("\\.(RData|rds)$", "", files, ignore.case = TRUE)
  stem <- sub("_para$", "", stem, ignore.case = TRUE)
  
  # {source}_{site, may contain _}_v{n}
  m <- regmatches(stem, regexec("^([^_]+)_(.+)_[vV](\\d+)$", stem))
  ok <- lengths(m) == 4L
  if (!any(ok)) return(NULL)
  
  data.frame(
    file   = file.path(PARAM_DIR, files[ok]),
    source = vapply(m[ok], `[`, "", 2L),
    site   = vapply(m[ok], `[`, "", 3L),
    visit  = as.integer(vapply(m[ok], `[`, "", 4L)),
    stringsAsFactors = FALSE
  )
}

source_label <- function(key) {
  if (key %in% names(SOURCE_LABELS)) SOURCE_LABELS[[key]] else key
}

site_label <- function(key) gsub("_", " ", key)

shannon <- function(counts) {            # counts: samples x taxa
  p <- counts / rowSums(counts)
  p[p == 0] <- NA
  -rowSums(p * log(p), na.rm = TRUE)
}

## ---- phylum aggregation --------------------------------

# Taxa are HMP OTU IDs, so the phylum has to come from a lookup exported
# out of the phyloseq tax_table(). See save_phylum_lookup.R.
PHYLUM_LOOKUP <- local({
  f <- file.path(PARAM_DIR, "otu2phylum.rds")
  if (file.exists(f)) readRDS(f) else NULL
})

to_relab <- function(mat) sweep(mat, 2, colSums(mat), "/")

# mat: taxa x samples of relative abundance -> long df, one row per
# phylum x sample, tagged with which side it came from.
phylum_long <- function(mat, src) {
  ph <- PHYLUM_LOOKUP[rownames(mat)]
  ph[is.na(ph)] <- "Unassigned"
  ph_mat <- rowsum(mat, group = ph)
  # check.names would turn HMP's numeric sample IDs into X700013549.
  as.data.frame(ph_mat, check.names = FALSE) %>%
    mutate(PHYLUM = rownames(ph_mat)) %>%
    pivot_longer(-PHYLUM, names_to = "Sample", values_to = "Abundance") %>%
    mutate(source = src)
}

## ---- inline help ---------------------------------------

# A label with a hoverable "?" after it.
labelled <- function(text, help) {
  tagList(text, tooltip(span("?", class = "help-dot"), help,
                        placement = "right"))
}

HELP_SOURCE <- paste(
  "The body region the reference samples came from. HMP collected",
  "specimens from several regions of the same subjects; each region",
  "has a distinct community, so the Dirichlet parameters are estimated",
  "separately and are not interchangeable."
)

HELP_SITE <- paste(
  "The specific location sampled within the region. HMP swabbed",
  "multiple sites per region (oral alone covers saliva, throat,",
  "tongue dorsum, and the plaques), and their compositions differ",
  "enough that each site gets its own parameter estimate."
)

HELP_VISIT <- paste(
  "The HMP visit number the reference samples were drawn from.",
  "Visits are repeat collections on the same subjects, so Visit 1 and",
  "Visit 2 hold different specimens from the same people. Parameters",
  "are estimated per visit; pick the one your study design mirrors."
)

HELP_CONF_DIFF <- paste(
  "Given as a fraction of ALL taxa, not of the differential ones.",
  "The taxa are drawn from the differential set, so this caps out at",
  "the 'Proportion differential' value above: any higher and every",
  "differential taxon is confounded and nothing further changes."
)

HELP_CONF_NONDIFF <- paste(
  "Given as a fraction of ALL taxa, not of the non-differential ones.",
  "The taxa are drawn from the non-differential set, so this plus",
  "'Proportion differential' must stay at or below 1."
)

HELP_CONF_EFF <- paste(
  "Size of the confounder's effect on the taxa flagged as confounded.",
  "The effects are drawn as rnorm(mean, sd), so leaving both this and",
  "Effect SD at 0 makes every effect exactly zero: the taxa are still",
  "flagged, but nothing about the counts changes. Setting Type to",
  "'none' forces both to 0 for the same reason."
)

INDEX <- scan_params()

## ---- ui ------------------------------------------------

ui <- page_navbar(
  title = "SurvSimMicro",
  id = "workflow",
  selected = "1 Microbiome counts",
  theme = bs_theme(version = 5, primary = "#176b87"),
  header = tags$style(HTML("
    .help-dot {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      width: 1.05em;
      height: 1.05em;
      margin-left: .35em;
      border: 1px solid currentColor;
      border-radius: 50%;
      font-size: .75em;
      line-height: 1;
      opacity: .55;
      cursor: help;
      vertical-align: text-top;
    }
    .help-dot:hover, .help-dot:focus { opacity: 1; }
    .workflow-intro {
      padding: .75rem 1rem;
      margin-bottom: 1rem;
      border-left: 4px solid var(--bs-primary);
      background: rgba(var(--bs-primary-rgb), .08);
    }
    .download-row { display: flex; flex-wrap: wrap; gap: .5rem; }
  ")),

  nav_panel(
    "1 Microbiome counts",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,

        h6("Reference data"),
        selectInput("source", labelled("Source", HELP_SOURCE), choices = NULL),
        selectInput("site", labelled("Body site", HELP_SITE), choices = NULL),
        selectInput("visit", labelled("Time point", HELP_VISIT), choices = NULL),
        helpText(textOutput("para_dims", inline = TRUE)),

        hr(),
        h6("Design"),
        numericInput("nSam", "Samples", 100, min = 4, step = 10),
        numericInput("nOTU", "Taxa", 200, min = 2, step = 50),
        numericInput("seed", "Random seed", 1),

        hr(),
        h6("Differential taxa"),
        sliderInput("diff.otu.pct", "Proportion differential",
                    0, 0.5, 0.1, step = 0.01),
        selectInput("diff.otu.direct", "Direction",
                    c("balanced", "unbalanced")),
        selectInput("diff.otu.mode", "Mode",
                    c("abundant", "rare", "mix", "user_specified")),
        conditionalPanel(
          "input['diff.otu.mode'] == 'user_specified'",
          textAreaInput("user_specified_otu", "Taxa names (one per line)",
                        rows = 3)
        ),

        hr(),
        h6("Covariate of interest"),
        selectInput("covariate.type", "Type", c("binary", "continuous")),
        conditionalPanel(
          "input['covariate.type'] == 'binary'",
          numericInput("grp.ratio", "Group ratio", 1, min = 0.1, step = 0.1)
        ),
        numericInput("covariate.eff.mean", "Effect mean", 1, step = 0.1),
        numericInput("covariate.eff.sd", "Effect SD", 0, min = 0, step = 0.1),

        hr(),
        h6("Confounding"),
        selectInput("confounder.type", "Type",
                    c("none", "binary", "continuous", "both")),
        conditionalPanel(
          "input['confounder.type'] != 'none'",
          sliderInput("conf.cov.cor", "Confounder-covariate correlation",
                      0, 0.95, 0.6, step = 0.05),
          sliderInput("conf.diff.otu.pct",
                      labelled("Confounded differential taxa", HELP_CONF_DIFF),
                      0, 0.5, 0, step = 0.01),
          sliderInput("conf.nondiff.otu.pct",
                      labelled("Confounded non-differential taxa", HELP_CONF_NONDIFF),
                      0, 1, 0.1, step = 0.05),
          numericInput("confounder.eff.mean",
                       labelled("Effect mean", HELP_CONF_EFF), 1, step = 0.1),
          numericInput("confounder.eff.sd", "Effect SD", 0, min = 0, step = 0.1)
        ),
        conditionalPanel(
          "input['confounder.type'] == 'none'",
          helpText("No confounder: no taxa are confounded.")
        ),

        hr(),
        h6("Sequencing"),
        numericInput("depth.mu", "Mean depth", 10000, min = 100, step = 1000),
        numericInput("depth.theta", "Depth dispersion (theta)", 5, min = 0.1),
        numericInput("depth.conf.factor", "Depth confounding", 0, step = 0.1),
        numericInput("error.sd", "Error SD", 0, min = 0, step = 0.1),

        hr(),
        actionButton("run", "Run count simulation", class = "btn-primary w-100")
      ),

      div(
        class = "workflow-intro",
        strong("Step 1: Generate microbiome counts."),
        " Inspect or download them here, then continue to Survival outcomes."
      ),
      navset_card_tab(
        nav_panel(
          "Summary",
          verbatimTextOutput("summary"),
          plotOutput("phylum_plot", height = "420px"),
          div(
            class = "download-row",
            downloadButton("dl_rds", "Download .rds"),
            downloadButton("dl_csv", "Download counts (.csv)")
          )
        ),
        nav_panel(
          "Counts",
          helpText("First 15 taxa x 10 samples."),
          tableOutput("preview")
        ),
        nav_panel(
          "Diagnostics",
          plotOutput("depth_plot", height = "260px"),
          plotOutput("shannon_plot", height = "260px"),
          plotOutput("pca_plot", height = "320px")
        )
      )
    )
  ),

  nav_panel(
    "2 Survival outcomes",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,

        h6("Input counts"),
        radioButtons(
          "surv_count_source", NULL,
          choices = c(
            "Use count simulation from Step 1" = "generated",
            "Upload a count matrix" = "upload"
          ),
          selected = "generated"
        ),
        conditionalPanel(
          "input.surv_count_source == 'upload'",
          fileInput("surv_count_file", "CSV or TSV file",
                    accept = c(".csv", ".tsv", ".txt")),
          radioButtons(
            "surv_upload_rows", "Rows represent",
            choices = c("Samples" = "samples", "Taxa" = "taxa"),
            selected = "samples", inline = TRUE
          ),
          helpText("The first column must contain row names.")
        ),
        helpText(textOutput("surv_count_dims", inline = TRUE)),

        hr(),
        h6("Survival design"),
        radioButtons(
          "surv_design", NULL,
          choices = c("One-stage" = "one_stage",
                      "Two-stage mixture cure" = "two_stage"),
          selected = "one_stage"
        ),
        selectInput(
          "surv_model", "Model",
          choices = c("Proportional hazards" = "ph",
                      "Accelerated failure time" = "aft")
        ),
        selectInput("surv_dist", "Distribution",
                    choices = c("Weibull" = "weibull",
                                "Exponential" = "exponential")),
        conditionalPanel(
          "input.surv_dist == 'weibull'",
          numericInput("surv_shape", "Weibull shape", 1.5, min = 0.05,
                       step = 0.1)
        ),
        numericInput("surv_median", "Baseline median survival", 5,
                     min = 0.01, step = 0.5),
        sliderInput("surv_target_event_prop",
                    "Target event / susceptible proportion",
                    min = 0.05, max = 0.95, value = 0.6, step = 0.05),
        conditionalPanel(
          "input.surv_design == 'one_stage'",
          helpText("With administrative censoring, this calibrates the observed event proportion.")
        ),
        conditionalPanel(
          "input.surv_design == 'two_stage'",
          helpText("For two-stage simulation, this sets the expected susceptible (noncured) proportion before censoring.")
        ),

        hr(),
        h6("Microbiome effects"),
        selectizeInput("surv_timing_taxa", "Timing-associated taxa",
                       choices = NULL, multiple = TRUE),
        numericInput("surv_timing_beta", "Timing effect per selected taxon",
                     0.2, step = 0.05),
        conditionalPanel(
          "input.surv_design == 'two_stage'",
          selectizeInput("surv_occurrence_taxa", "Occurrence-associated taxa",
                         choices = NULL, multiple = TRUE),
          numericInput("surv_occurrence_beta",
                       "Occurrence effect per selected taxon", 0.2, step = 0.05)
        ),

        hr(),
        h6("Censoring"),
        selectInput(
          "surv_censor_type", "Type",
          choices = c(
            "Administrative" = "administrative",
            "Random exponential" = "random",
            "Administrative + random" = "administrative_plus_random",
            "None" = "none"
          )
        ),
        conditionalPanel(
          "(input.surv_design == 'two_stage' && input.surv_censor_type == 'administrative') || input.surv_censor_type == 'administrative_plus_random'",
          numericInput("surv_admin_time", "Administrative cutoff", 10,
                       min = 0.01, step = 0.5)
        ),
        conditionalPanel(
          "input.surv_censor_type == 'random' || input.surv_censor_type == 'administrative_plus_random'",
          numericInput("surv_random_rate", "Random censoring rate", 0.1,
                       min = 0.001, step = 0.01)
        ),
        numericInput("surv_seed", "Random seed", 101, step = 1),

        hr(),
        actionButton("run_survival", "Generate survival data",
                     class = "btn-primary w-100")
      ),

      if (!requireNamespace("simsurv", quietly = TRUE))
        div(
          class = "alert alert-warning",
          strong("Survival dependency missing. "),
          "Install the R package ", code("simsurv"),
          " and restart the app. The count simulator remains available."
        ),
      div(
        class = "workflow-intro",
        strong("Step 2: Generate survival outcomes."),
        " Use the Step 1 counts or upload an independent sample-by-taxa matrix."
      ),
      navset_card_tab(
        nav_panel(
          "Overview",
          verbatimTextOutput("surv_summary"),
          plotOutput("surv_km_plot", height = "380px")
        ),
        nav_panel(
          "Data",
          helpText("First 10 subjects. Risk score is the timing linear predictor."),
          tableOutput("surv_preview"),
          div(
            class = "download-row",
            downloadButton("dl_surv_csv", "Download survival outcomes"),
            downloadButton("dl_combined_csv", "Download counts + outcomes"),
            downloadButton("dl_surv_rds", "Download full .rds")
          )
        ),
        nav_panel(
          "Diagnostics",
          plotOutput("surv_time_plot", height = "280px"),
          plotOutput("surv_risk_plot", height = "320px")
        )
      )
    )
  )
)

## ---- server --------------------------------------------

server <- function(input, output, session) {
  
  # No params folder, or nothing parseable in it: say so once, clearly.
  if (is.null(INDEX)) {
    showModal(modalDialog(
      title = "No parameter files found",
      tagList(
        p(paste0("Put your .RData or .rds files in ", PARAM_DIR, "/.")),
        p("Names must follow {source}_{site}_v{n}_para, for example ",
          code("oral_Saliva_v2_para.RData"), ".")
      ),
      easyClose = TRUE
    ))
  }
  
  ## Cascade: source -> site -> visit
  observe({
    req(INDEX)
    keys <- unique(INDEX$source)
    keys <- keys[order(match(keys, names(SOURCE_LABELS)), keys)]
    updateSelectInput(session, "source",
                      choices = stats::setNames(keys, vapply(keys, source_label, "")))
  })
  
  sites <- reactive({
    req(INDEX, input$source)
    sort(unique(INDEX$site[INDEX$source == input$source]))
  })
  
  observeEvent(sites(), {
    s <- sites()
    req(length(s) > 0)
    sel <- if (isTruthy(input$site) && input$site %in% s) input$site else s[1]
    updateSelectInput(session, "site", selected = sel,
                      choices = stats::setNames(s, vapply(s, site_label, "")))
  })
  
  visits <- reactive({
    req(INDEX, input$source, input$site)
    # input$site can lag a source change by one flush; ignore the stale pair.
    req(input$site %in% sites())
    v <- INDEX$visit[INDEX$source == input$source & INDEX$site == input$site]
    sort(unique(v))
  })
  
  # Time points are per site: keep the visit when it exists there too.
  observeEvent(visits(), {
    v <- visits()
    req(length(v) > 0)
    sel <- if (isTruthy(input$visit) && input$visit %in% v) input$visit else v[1]
    updateSelectInput(session, "visit", selected = sel,
                      choices = stats::setNames(v, paste("Visit", v)))
  })
  
  para <- reactive({
    req(INDEX, input$source, input$site, input$visit)
    row <- INDEX$source == input$source &
      INDEX$site   == input$site &
      INDEX$visit  == as.integer(input$visit)
    # Mid-cascade the inputs can name a combination that has no file.
    req(any(row))
    
    path <- INDEX$file[which(row)[1]]
    p    <- load_para(path)
    validate(need(!is.null(p),
                  paste0(basename(path),
                         " holds no EstPara() result (needs mu and ref.otu.tab).")))
    p
  })
  
  output$para_dims <- renderText({
    p <- try(para(), silent = TRUE)
    if (inherits(p, "try-error")) return("")
    sprintf("Reference: %d taxa x %d samples",
            nrow(p$ref.otu.tab), ncol(p$ref.otu.tab))
  })
  
  # Ceilings differ a lot by body site, so clamp rather than let
  # SimulateMSeqU() slice past the end and return NA rows.
  observeEvent(para(), {
    updateNumericInput(session, "nOTU", max = nrow(para()$ref.otu.tab))
  })
  
  sim <- eventReactive(input$run, {
    p <- para()
    
    validate(
      need(input$nSam >= 4, "Use at least 4 samples."),
      need(input$nOTU >= 2, "Use at least 2 taxa."),
      need(input$nOTU <= nrow(p$ref.otu.tab),
           sprintf("%s at visit %s has only %d taxa.",
                   site_label(input$site), input$visit, nrow(p$ref.otu.tab)))
    )
    
    user_otu <- NULL
    if (input$diff.otu.mode == "user_specified") {
      user_otu <- trimws(strsplit(input$user_specified_otu, "\n")[[1]])
      user_otu <- user_otu[nzchar(user_otu)]
      validate(need(length(user_otu) > 0, "List at least one taxon name."))
    }
    
    # "none" zeroes the effect but not the selection, so the function would
    # still flag taxa as confounded by a confounder that does nothing.
    # Force the pools empty here so "none" means none.
    conf_off  <- input$confounder.type == "none"
    pct_cdiff <- if (conf_off) 0 else input$conf.diff.otu.pct
    pct_cnd   <- if (conf_off) 0 else input$conf.nondiff.otu.pct
    
    # SimulateMSeqU() scales both confounding counts by nOTU rather than by
    # the pool each is drawn from, so check them against the real pools.
    n_diff <- if (input$diff.otu.mode == "user_specified") {
      length(user_otu)
    } else {
      round(input$diff.otu.pct * input$nOTU)
    }
    n_conf_diff    <- round(input$nOTU * pct_cdiff)
    n_conf_nondiff <- round(input$nOTU * pct_cnd)
    
    validate(
      # Line 188 samples from the non-differential taxa without a size check.
      need(n_conf_nondiff <= input$nOTU - n_diff,
           sprintf(paste("Confounded non-differential taxa asks for %d of only",
                         "%d non-differential taxa. Lower it below %.2f, or",
                         "lower the proportion differential."),
                   n_conf_nondiff, input$nOTU - n_diff,
                   (input$nOTU - n_diff) / input$nOTU))
    )
    
    if (n_conf_diff > n_diff) {
      showNotification(
        sprintf(paste("Confounded differential taxa (%.2f of all taxa = %d)",
                      "exceeds the %d differential taxa available, so all of",
                      "them are confounded. Values above %.2f change nothing."),
                pct_cdiff, n_conf_diff, n_diff, input$diff.otu.pct),
        type = "warning", duration = 10)
    }
    
    set.seed(input$seed)
    
    # SimulateMSeqU() has no defaults for these two and always uses them.
    epsilon   <- rnorm(input$nSam)
    cont.conf <- rnorm(input$nSam)
    
    withProgress(message = "Simulating", value = 0.5, {
      SimulateMSeqU(
        para                 = p,
        nSam                 = input$nSam,
        nOTU                 = input$nOTU,
        diff.otu.pct         = input$diff.otu.pct,
        diff.otu.direct      = input$diff.otu.direct,
        diff.otu.mode        = input$diff.otu.mode,
        user_specified_otu   = user_otu,
        covariate.type       = input$covariate.type,
        grp.ratio            = input$grp.ratio,
        covariate.eff.mean   = input$covariate.eff.mean,
        covariate.eff.sd     = input$covariate.eff.sd,
        confounder.type      = input$confounder.type,
        conf.cov.cor         = input$conf.cov.cor,
        conf.diff.otu.pct    = pct_cdiff,
        conf.nondiff.otu.pct = pct_cnd,
        confounder.eff.mean  = input$confounder.eff.mean,
        confounder.eff.sd    = input$confounder.eff.sd,
        error.sd             = input$error.sd,
        depth.mu             = input$depth.mu,
        depth.theta          = input$depth.theta,
        depth.conf.factor    = input$depth.conf.factor,
        cont.conf            = cont.conf,
        epsilon              = epsilon
      )
    })
  })
  
  output$summary <- renderPrint({
    s  <- sim()
    ct <- s$otu.tab.sim
    # conf.otu.ind is the union of both pools; report the split, since one
    # total can't be read back against the two sliders that produced it.
    n_cd  <- sum(s$conf.otu.ind & s$diff.otu.ind)
    n_cnd <- sum(s$conf.otu.ind & !s$diff.otu.ind)
    cat("Source        :", source_label(input$source), "\n")
    cat("Body site     :", site_label(input$site), "\n")
    cat("Time point    : Visit", input$visit, "\n")
    cat("Seed          :", input$seed, "\n\n")
    cat("Counts        :", nrow(ct), "taxa x", ncol(ct), "samples\n")
    cat("Differential  :", sum(s$diff.otu.ind), "taxa\n")
    if (input$confounder.type == "none") {
      cat("Confounded    : 0 taxa (confounding off)\n")
    } else {
      cat("Confounded    :", sum(s$conf.otu.ind), "taxa",
          sprintf("(%d differential + %d non-differential)", n_cd, n_cnd), "\n")
      # eta.conf ~ rnorm(mean, sd); both zero => zero matrix => flags only.
      if (input$confounder.eff.mean == 0 && input$confounder.eff.sd == 0) {
        cat("                (flagged only: confounder effect is 0,",
            "so these taxa are unaffected)\n")
      }
    }
    cat("Zero fraction :", sprintf("%.1f%%", 100 * mean(ct == 0)), "\n\n")
    cat("Library size:\n"); print(summary(colSums(ct)))
    cat("\nCovariate:\n"); print(summary(as.vector(s$covariate)))
  })
  
  output$preview <- renderTable({
    ct <- sim()$otu.tab.sim
    ct[seq_len(min(15, nrow(ct))), seq_len(min(10, ncol(ct))), drop = FALSE]
  }, rownames = TRUE)
  
  output$depth_plot <- renderPlot({
    hist(colSums(sim()$otu.tab.sim), breaks = 30, col = "grey80",
         border = "white", main = "Library size", xlab = "Total counts")
  })
  
  output$shannon_plot <- renderPlot({
    s   <- sim()
    div <- shannon(t(s$otu.tab.sim))
    if (input$covariate.type == "binary") {
      boxplot(div ~ factor(as.vector(s$covariate)),
              col = "grey80", main = "Shannon diversity by group",
              xlab = "Group", ylab = "Shannon")
    } else {
      plot(as.vector(s$covariate), div, pch = 16, col = "grey40",
           main = "Shannon diversity vs covariate",
           xlab = "Covariate", ylab = "Shannon")
    }
  })
  
  output$pca_plot <- renderPlot({
    s   <- sim()
    ct  <- t(s$otu.tab.sim) + 0.5             # samples x taxa, pseudocount
    clr <- log(ct) - rowMeans(log(ct))        # CLR transform
    pc  <- prcomp(clr, center = TRUE)
    ve  <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    grp <- as.vector(s$covariate)
    col <- if (input$covariate.type == "binary") {
      c("#4C72B0", "#DD8452")[factor(grp)]
    } else {
      r <- range(grp)
      grey(0.1 + 0.7 * (grp - r[1]) / diff(r))
    }
    plot(pc$x[, 1], pc$x[, 2], pch = 16, col = col,
         main = "PCA of CLR-transformed counts",
         xlab = sprintf("PC1 (%.1f%%)", ve[1]),
         ylab = sprintf("PC2 (%.1f%%)", ve[2]))
  })
  
  # Reference vs simulated composition, one bar per sample. The reference is
  # subset to the same 1:nOTU rows SimulateMSeqU() used, so the two sides are
  # directly comparable rather than differing because of the nOTU setting.
  output$phylum_plot <- renderPlot({
    s      <- sim()
    ct_sim <- s$otu.tab.sim
    ct_ref <- as.matrix(para()$ref.otu.tab)[seq_len(nrow(ct_sim)), , drop = FALSE]
    
    validate(need(
      !is.null(PHYLUM_LOOKUP),
      paste0("No phylum lookup found. Run save_phylum_lookup.R to write ",
             PARAM_DIR, "/otu2phylum.rds.")
    ))
    
    ct_ref <- ct_ref[, colSums(ct_ref) > 0, drop = FALSE]
    ct_sim <- ct_sim[, colSums(ct_sim) > 0, drop = FALSE]
    
    # A silent miss renders as one big "Unassigned" band, so say it outright.
    hit <- mean(rownames(ct_sim) %in% names(PHYLUM_LOOKUP))
    if (hit < 0.9) {
      showNotification(
        sprintf(paste("Only %.0f%% of taxa matched the phylum lookup. The para",
                      "files and otu2phylum.rds may come from different",
                      "datasets (e.g. V13 vs V35)."), 100 * hit),
        type = "warning", duration = 10)
    }
    
    plot_df <- bind_rows(
      phylum_long(to_relab(ct_ref), "reference"),
      phylum_long(to_relab(ct_sim), "simulated")
    )
    
    top_phyla <- plot_df %>%
      group_by(PHYLUM) %>%
      summarise(m = mean(Abundance), .groups = "drop") %>%
      arrange(desc(m)) %>%
      slice_head(n = 10) %>%
      pull(PHYLUM)
    
    plot_df <- plot_df %>%
      mutate(PHYLUM2 = ifelse(PHYLUM %in% top_phyla, PHYLUM, "Other"))
    
    ggplot(plot_df, aes(x = Sample, y = Abundance, fill = PHYLUM2)) +
      geom_bar(stat = "identity") +
      facet_grid(. ~ source, scales = "free_x", space = "free_x") +
      scale_y_continuous(labels = percent) +
      theme_bw() +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid   = element_blank()
      ) +
      labs(
        title = sprintf("%s, %s, visit %s: phylum-level relative abundance",
                        source_label(input$source), site_label(input$site),
                        input$visit),
        x = "Samples", y = "Relative abundance", fill = "Phylum"
      )
  })
  
  output$dl_rds <- downloadHandler(
    filename = function() {
      sprintf("sim_%s_%s_v%s_seed%s.rds", input$source, input$site,
              input$visit, input$seed)
    },
    content = function(file) saveRDS(sim(), file)
  )
  
  output$dl_csv <- downloadHandler(
    filename = function() {
      sprintf("counts_%s_%s_v%s_seed%s.csv", input$source, input$site,
              input$visit, input$seed)
    },
    content = function(file) write.csv(sim()$otu.tab.sim, file)
  )

  ## ---- survival workflow ---------------------------------

  uploaded_survival_counts <- reactive({
    req(input$surv_count_file)

    path <- input$surv_count_file$datapath
    ext  <- tolower(tools::file_ext(input$surv_count_file$name))
    sep  <- if (ext == "csv") "," else "\t"

    dat <- tryCatch(
      utils::read.table(
        path, header = TRUE, row.names = 1, sep = sep,
        check.names = FALSE, comment.char = "", quote = "\"",
        stringsAsFactors = FALSE
      ),
      error = function(e) e
    )

    read_error <- inherits(dat, "error")
    validate(need(!read_error,
                  if (read_error) paste("Could not read the uploaded file:",
                                        conditionMessage(dat)) else ""))
    validate(need(nrow(dat) > 0 && ncol(dat) > 0,
                  "The uploaded file contains no count columns."))
    validate(need(all(vapply(dat, is.numeric, logical(1))),
                  "Every count column must be numeric."))

    m <- as.matrix(dat)
    if (identical(input$surv_upload_rows, "taxa")) m <- t(m)
    m
  })

  # All downstream survival functions use samples x taxa. SimulateMSeqU()
  # returns taxa x samples, so Step 1 data are transposed exactly once here.
  survival_counts <- reactive({
    if (identical(input$surv_count_source, "generated")) {
      validate(need(input$run > 0,
                    "Run the count simulation in Step 1 before using it here."))
      m <- t(sim()$otu.tab.sim)
    } else {
      m <- uploaded_survival_counts()
    }

    validate(
      need(is.numeric(m), "The survival count matrix must be numeric."),
      need(nrow(m) >= 4, "Use at least 4 samples for survival simulation."),
      need(ncol(m) >= 2, "Use at least 2 taxa for survival simulation."),
      need(all(is.finite(m)), "Counts cannot contain NA, NaN, or Inf."),
      need(all(m >= 0), "Counts cannot be negative.")
    )

    if (is.null(rownames(m))) rownames(m) <- paste0("Sample", seq_len(nrow(m)))
    if (is.null(colnames(m))) colnames(m) <- paste0("Taxon", seq_len(ncol(m)))
    validate(need(!anyDuplicated(colnames(m)),
                  "Taxon names must be unique."))
    m
  })

  output$surv_count_dims <- renderText({
    m <- try(survival_counts(), silent = TRUE)
    if (inherits(m, "try-error")) {
      if (identical(input$surv_count_source, "generated"))
        return("Waiting for a Step 1 count simulation.")
      return("Upload a matrix to continue.")
    }
    sprintf("Available: %d samples x %d taxa", nrow(m), ncol(m))
  })

  observe({
    m <- try(survival_counts(), silent = TRUE)
    if (inherits(m, "try-error")) return()

    taxa <- colnames(m)
    n_timing <- min(10L, length(taxa))
    timing_default <- taxa[seq_len(n_timing)]

    n_occ <- min(5L, length(taxa))
    occ_start <- if (length(taxa) >= n_timing + n_occ) n_timing + 1L else 1L
    occurrence_default <- taxa[seq.int(occ_start,
                                       length.out = min(n_occ,
                                                        length(taxa) - occ_start + 1L))]

    updateSelectizeInput(
      session, "surv_timing_taxa", choices = taxa,
      selected = timing_default, server = TRUE
    )
    updateSelectizeInput(
      session, "surv_occurrence_taxa", choices = taxa,
      selected = occurrence_default, server = TRUE
    )
  })

  survival_run <- eventReactive(input$run_survival, {
    validate(need(requireNamespace("simsurv", quietly = TRUE),
                  paste0("The survival engine requires the 'simsurv' package. ",
                         "Install it with install.packages('simsurv'), then restart the app.")))

    counts <- survival_counts()
    timing_taxa <- input$surv_timing_taxa
    validate(
      need(length(timing_taxa) > 0,
           "Select at least one timing-associated taxon."),
      need(all(timing_taxa %in% colnames(counts)),
           "One or more selected timing taxa are not in the current count matrix.")
    )

    predictor_timing <- list(
      index = timing_taxa,
      betas = rep(input$surv_timing_beta, length(timing_taxa))
    )

    predictor_occurrence <- list()
    occurrence_taxa <- character()
    if (identical(input$surv_design, "two_stage")) {
      occurrence_taxa <- input$surv_occurrence_taxa
      validate(
        need(length(occurrence_taxa) > 0,
             "Select at least one occurrence-associated taxon."),
        need(all(occurrence_taxa %in% colnames(counts)),
             "One or more selected occurrence taxa are not in the current count matrix.")
      )
      predictor_occurrence <- list(
        index = occurrence_taxa,
        betas = rep(input$surv_occurrence_beta, length(occurrence_taxa))
      )
    }

    tte_args <- list(
      model = input$surv_model,
      dist = input$surv_dist,
      params = list(shape = if (input$surv_dist == "weibull")
        input$surv_shape else 1),
      target_median = input$surv_median
    )

    censor_args <- switch(
      input$surv_censor_type,
      none = list(type = "none"),
      administrative = {
        if (identical(input$surv_design, "one_stage")) {
          list(type = "administrative")
        } else {
          list(type = "administrative",
               administrative_time = input$surv_admin_time)
        }
      },
      random = list(
        type = "random",
        random = list(dist = "exponential", rate = input$surv_random_rate)
      ),
      administrative_plus_random = list(
        type = "administrative_plus_random",
        administrative_time = input$surv_admin_time,
        random = list(dist = "exponential", rate = input$surv_random_rate)
      )
    )

    fit <- withProgress(message = "Generating survival outcomes", value = 0.5, {
      if (identical(input$surv_design, "one_stage")) {
        simulate_tte(
          counts = counts,
          design = "one_stage",
          predictor_args = predictor_timing,
          tte_args = tte_args,
          censor_args = censor_args,
          target_event_prop = input$surv_target_event_prop,
          seed = input$surv_seed
        )
      } else {
        simulate_tte(
          counts = counts,
          design = "two_stage",
          predictor_occurrence = predictor_occurrence,
          predictor_timing = predictor_timing,
          tte_args = tte_args,
          censor_args = censor_args,
          target_event_prop = input$surv_target_event_prop,
          seed = input$surv_seed
        )
      }
    })

    list(
      fit = fit,
      counts = counts,
      sample_id = rownames(counts),
      timing_taxa = timing_taxa,
      occurrence_taxa = occurrence_taxa,
      model = input$surv_model,
      dist = input$surv_dist,
      censor_type = input$surv_censor_type,
      seed = input$surv_seed
    )
  })

  survival_table <- reactive({
    run <- survival_run()
    fit <- run$fit
    risk <- if (is.list(fit$eta)) fit$eta$timing else fit$eta

    out <- data.frame(
      sample_id = run$sample_id,
      time = fit$outcome$time,
      status = fit$outcome$status,
      risk_score = risk,
      stringsAsFactors = FALSE
    )

    if (!is.null(fit$tte$pi)) {
      out$susceptibility <- fit$tte$pi
      out$event_eligible <- fit$tte$event_eligible
    }
    out
  })

  output$surv_summary <- renderPrint({
    run <- survival_run()
    fit <- run$fit
    tab <- survival_table()

    cat("Design         :", if (fit$call$design == "one_stage")
      "One-stage" else "Two-stage mixture cure", "\n")
    cat("Model          :", toupper(run$model), "\n")
    cat("Distribution   :", tools::toTitleCase(run$dist), "\n")
    cat("Censoring      :", gsub("_", " ", run$censor_type), "\n")
    cat("Seed           :", run$seed, "\n\n")
    cat("Input counts   :", nrow(run$counts), "samples x",
        ncol(run$counts), "taxa\n")
    cat("Timing taxa    :", length(run$timing_taxa), "\n")
    if (fit$call$design == "two_stage")
      cat("Occurrence taxa:", length(run$occurrence_taxa), "\n")
    cat("Observed events:", sum(tab$status), "of", nrow(tab),
        sprintf("(%.1f%%)\n", 100 * mean(tab$status)))
    cat("Median time    :", round(stats::median(tab$time), 3), "\n")
    if (!is.null(fit$tte$pi)) {
      cat("Mean susceptible probability:",
          round(mean(fit$tte$pi), 3), "\n")
      cat("Realized susceptible subjects:",
          sum(fit$tte$event_eligible), "of", nrow(tab), "\n")
    }
  })

  output$surv_preview <- renderTable({
    utils::head(survival_table(), 10)
  }, digits = 4, rownames = FALSE)

  output$surv_km_plot <- renderPlot({
    tab <- survival_table()
    validate(need(nrow(tab) > 1, "At least two subjects are required."))

    if (stats::sd(tab$risk_score) > 0) {
      risk_group <- factor(
        ifelse(tab$risk_score <= stats::median(tab$risk_score),
               "Lower microbiome risk", "Higher microbiome risk"),
        levels = c("Lower microbiome risk", "Higher microbiome risk")
      )
    } else {
      risk_group <- factor(rep("All subjects", nrow(tab)))
    }

    km <- survival::survfit(
      survival::Surv(time, status) ~ risk_group,
      data = transform(tab, risk_group = risk_group)
    )
    cols <- c("#176b87", "#d77b37")[seq_along(levels(risk_group))]
    plot(km, col = cols, lwd = 2, mark.time = TRUE,
         xlab = "Time", ylab = "Survival probability",
         main = "Kaplan-Meier curves by microbiome risk score")
    if (nlevels(risk_group) > 1L) {
      legend("bottomleft", legend = levels(risk_group), col = cols,
             lwd = 2, bty = "n")
    }
  })

  output$surv_time_plot <- renderPlot({
    tab <- survival_table()
    ggplot(tab, aes(x = time, fill = factor(status))) +
      geom_histogram(position = "identity", alpha = 0.65, bins = 30) +
      scale_fill_manual(values = c("0" = "grey65", "1" = "#176b87"),
                        labels = c("Censored", "Event")) +
      theme_bw() +
      labs(title = "Observed follow-up times", x = "Time", y = "Subjects",
           fill = "Outcome")
  })

  output$surv_risk_plot <- renderPlot({
    tab <- survival_table()
    ggplot(tab, aes(x = risk_score, y = time, color = factor(status))) +
      geom_point(size = 2, alpha = 0.75) +
      geom_smooth(method = "loess", se = FALSE, color = "grey35") +
      scale_color_manual(values = c("0" = "grey55", "1" = "#176b87"),
                         labels = c("Censored", "Event")) +
      theme_bw() +
      labs(title = "Follow-up time versus microbiome risk score",
           x = "Timing linear predictor", y = "Observed time",
           color = "Outcome")
  })

  output$dl_surv_csv <- downloadHandler(
    filename = function() sprintf("survival_outcomes_seed%s.csv", input$surv_seed),
    content = function(file) utils::write.csv(survival_table(), file,
                                               row.names = FALSE)
  )

  output$dl_combined_csv <- downloadHandler(
    filename = function() sprintf("counts_survival_seed%s.csv", input$surv_seed),
    content = function(file) {
      run <- survival_run()
      combined <- cbind(survival_table(),
                        as.data.frame(run$counts, check.names = FALSE))
      utils::write.csv(combined, file, row.names = FALSE)
    }
  )

  output$dl_surv_rds <- downloadHandler(
    filename = function() sprintf("survival_simulation_seed%s.rds", input$surv_seed),
    content = function(file) saveRDS(survival_run(), file)
  )
}

shinyApp(ui, server)
