# Example04: G-CSF PK-PD population FIM vs PopED Design 1 RSE.
utils = system.file("vignette-scripts", "pfim-vignette-utils.R", package = "PFIM")
if (!nzchar(utils)) stop("pfim-vignette-utils.R not found.", call. = FALSE)
source(utils, local = TRUE)
library(PFIM)

paths = pfimVignetteSetupPaths()
plotOptions = list(unitTime = c("hour"), unitOutcomes = c("ng/mL", "10^3/uL"))

# ---- Algebraic helpers (QSS TMDD + myelopoiesis) -----------------------------

.pfimGcsfNt = function(prefix = "") {
  paste0("(", paste0(prefix, c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "NB"),
                     collapse = "+"), ")")
}

.pfimGcsfCp = function() {
  nt = .pfimGcsfNt()
  bb = paste0("(KSI*", nt, "-CENT/VD+KD)")
  paste0("(0.5*(-", bb, "+sqrt((", bb, ")**2+4*KD*CENT/VD)))")
}

.pfimGcsfBaselineICs = function() {
  h10 = "(1+SM1*BAS/(SC1+BAS))"
  h20 = "(1+SM2*BAS/(SC1+BAS))"
  h30 = "(1+SM3*BAS/(SC1+BAS))"
  kinb = paste0("(KMT*NB0/", h10, ")")
  den = paste0("(KBB1*", h30, "+KTT*", h20, ")")
  bm10 = paste0("(", kinb, "*", h10, "/", den, ")")
  bm20 = paste0("(KTT*", h20, "*", bm10, "/", den, ")")
  bm30 = paste0("(KTT*", h20, "*", bm20, "/", den, ")")
  bm40 = paste0("(KTT*", h20, "*", bm30, "/", den, ")")
  bm50 = paste0("(KTT*", h20, "*", bm40, "/", den, ")")
  bm60 = paste0("(KTT*", h20, "*", bm50, "/", den, ")")
  bm70 = paste0("(KTT*", h20, "*", bm60, "/", den, ")")
  bm80 = paste0("(KTT*", h20, "*", bm70, "/", den, ")")
  bm90 = paste0("(KTT*", h20, "*", bm80, "/", den, ")")
  nt0 = paste0("(", bm10, "+", bm20, "+", bm30, "+", bm40, "+", bm50, "+",
                bm60, "+", bm70, "+", bm80, "+", bm90, "+NB0)")
  ac0 = "(BAS*VD)"
  adr0 = paste0("(KSI*", nt0, "*", ac0, "/(KD+BAS))")
  list(
    ABS = "FF*dose_ABS",
    CENT = paste0("(", ac0, "+", adr0, ")"),
    B1 = bm10, B2 = bm20, B3 = bm30, B4 = bm40, B5 = bm50,
    B6 = bm60, B7 = bm70, B8 = bm80, B9 = bm90,
    NB = "NB0"
  )
}

.pfimGcsfModelEquations = function() {
  cp = .pfimGcsfCp()
  nt = .pfimGcsfNt()
  ac = paste0("(", cp, "*VD)")
  adr = paste0("(KSI*", nt, "*", ac, "/(KD+", cp, "))")
  h1 = paste0("(1+SM1*", cp, "/(SC1+", cp, "))")
  h2 = paste0("(1+SM2*", cp, "/(SC1+", cp, "))")
  h3 = paste0("(1+SM3*", cp, "/(SC1+", cp, "))")
  h10 = "(1+SM1*BAS/(SC1+BAS))"
  h20 = "(1+SM2*BAS/(SC1+BAS))"
  h30 = "(1+SM3*BAS/(SC1+BAS))"
  kinb = paste0("(KMT*NB0/", h10, ")")
  den0 = paste0("(KBB1*", h30, "+KTT*", h20, ")")
  bm10 = paste0("(", kinb, "*", h10, "/", den0, ")")
  bm20 = paste0("(KTT*", h20, "*", bm10, "/", den0, ")")
  bm30 = paste0("(KTT*", h20, "*", bm20, "/", den0, ")")
  bm40 = paste0("(KTT*", h20, "*", bm30, "/", den0, ")")
  bm50 = paste0("(KTT*", h20, "*", bm40, "/", den0, ")")
  bm60 = paste0("(KTT*", h20, "*", bm50, "/", den0, ")")
  bm70 = paste0("(KTT*", h20, "*", bm60, "/", den0, ")")
  bm80 = paste0("(KTT*", h20, "*", bm70, "/", den0, ")")
  bm90 = paste0("(KTT*", h20, "*", bm80, "/", den0, ")")
  nt0 = paste0("(", bm10, "+", bm20, "+", bm30, "+", bm40, "+", bm50, "+",
                bm60, "+", bm70, "+", bm80, "+", bm90, "+NB0)")
  ac0 = "(BAS*VD)"
  adr0 = paste0("(KSI*", nt0, "*", ac0, "/(KD+BAS))")
  kin = paste0("(KEL*", ac0, "+KINT*", adr0, ")")
  bsum = "(B1+B2+B3+B4+B5+B6+B7+B8+B9)"
  list(
    Deriv_ABS  = "-KA*ABS",
    Deriv_CENT = paste0("KA*ABS+", kin, "-KEL*", ac, "-KINT*", adr),
    Deriv_B1   = paste0(kinb, "*", h1, "*(B1/", bm10, ")-KBB1*", h3, "*B1-KTT*", h2, "*B1"),
    Deriv_B2   = paste0("KTT*", h2, "*B1-KBB1*", h3, "*B2-KTT*", h2, "*B2"),
    Deriv_B3   = paste0("KTT*", h2, "*B2-KBB1*", h3, "*B3-KTT*", h2, "*B3"),
    Deriv_B4   = paste0("KTT*", h2, "*B3-KBB1*", h3, "*B4-KTT*", h2, "*B4"),
    Deriv_B5   = paste0("KTT*", h2, "*B4-KBB1*", h3, "*B5-KTT*", h2, "*B5"),
    Deriv_B6   = paste0("KTT*", h2, "*B5-KBB1*", h3, "*B6-KTT*", h2, "*B6"),
    Deriv_B7   = paste0("KTT*", h2, "*B6-KBB1*", h3, "*B7-KTT*", h2, "*B7"),
    Deriv_B8   = paste0("KTT*", h2, "*B7-KBB1*", h3, "*B8-KTT*", h2, "*B8"),
    Deriv_B9   = paste0("KTT*", h2, "*B8-KBB1*", h3, "*B9-KTT*", h2, "*B9"),
    Deriv_NB   = paste0("KTT*", h2, "*B9+KBB1*", h3, "*", bsum, "-KMT*NB")
  )
}

# ---- Design matching PopED create.poped.database(...) -----------------------

pk_dense = c(
  0.167, 0.25, 0.333, 0.5, 0.667, 0.75, 0.833, 1, 1.5, 2, 3, 4, 5,
  6, 8, 10, 12, 14, 16, 18, 20, 24
)
anc_dense = c(
  0.333, 0.5, 0.667, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14,
  16, 18, 20, 24
)
pk_times = c(pk_dense, pk_dense + 6 * 24, 7 * 24 + c(4, 24, 48))
anc_times = c(
  anc_dense, (2:6) * 24, anc_dense + 6 * 24, 7 * 24 + c(4, 24, 48, 72)
)

WT = 75
TAU = 24
dose_times = seq(0, 6 * TAU, by = TAU)
dose_ug_kg = c(1, 3, 10)

sigma_pk_prop = sqrt(2.53e-01)
sigma_anc_prop = sqrt(2.27e-02)
sigma_anc_add = sqrt(2.10e+00)

mp = function(name, mu, omega = 0, fixed_mu = FALSE) {
  ModelParameter(
    name = name,
    distribution = LogNormal(mu = mu, omega = omega),
    fixedMu = fixed_mu
  )
}

modelParameters = list(
  mp("FF", 6.26e-01, 0, TRUE),
  mp("KA", 6.42e-01, 0),
  mp("FR", 1.00e+00, 0, TRUE),
  mp("D2", 6.77e+00, 0, TRUE),
  mp("KEL", 1.48e-01, sqrt(3.12e-01)),
  mp("VD", 2.56e+00, sqrt(3.28e-01)),
  mp("KD", 1.27e+00, 0),
  mp("KINT", 1.01e-01, 0),
  mp("KSI", 2.11e-01, sqrt(2.24e-01)),
  mp("KOFF", 0.00e+00, 0, TRUE),
  mp("KMT", 7.23e-02, 0),
  mp("KBB1", 0.00e+00, 0, TRUE),
  mp("KTT", 1.02e-02, 0),
  mp("NB0", 1.65e+00, sqrt(2.98e-01)),
  mp("SC1", 3.21e+00, sqrt(8.03e-01)),
  mp("SM1", 3.43e+01, sqrt(1.28e-02)),
  mp("SM2", 3.23e+01, 0),
  mp("SM3", 0.00e+00, 0, TRUE),
  mp("BAS", 0.0246, 0, TRUE)
)

modelError = list(
  Combined2(
    output = "RespPK",
    sigmaInter = 0,
    sigmaSlope = sigma_pk_prop,
    sigmaInterFixed = TRUE
  ),
  Combined2(
    output = "RespPD",
    sigmaInter = sigma_anc_add,
    sigmaSlope = sigma_anc_prop
  )
)

samplingPK = SamplingTimes(outcome = "RespPK", samplings = pk_times)
samplingPD = SamplingTimes(outcome = "RespPD", samplings = anc_times)
ics = .pfimGcsfBaselineICs()

mk_arm = function(name, dose_ug_per_kg) {
  amt = dose_ug_per_kg * WT
  Arm(
    name = name,
    size = 10,
    administrations = list(
      Administration(
        outcome = "ABS",
        timeDose = dose_times,
        dose = rep(amt, length(dose_times))
      )
    ),
    samplingTimes = list(samplingPK, samplingPD),
    initialConditions = ics
  )
}

design1 = Design(
  name = "gcsf_design1",
  arms = list(
    mk_arm("dose_1ugkg", dose_ug_kg[[1L]]),
    mk_arm("dose_3ugkg", dose_ug_kg[[2L]]),
    mk_arm("dose_10ugkg", dose_ug_kg[[3L]])
  )
)

modelEquations = .pfimGcsfModelEquations()
cp_out = .pfimGcsfCp()
outputs = list(RespPK = cp_out, RespPD = "NB")

# ---- Evaluation -------------------------------------------------------------
# Heavy 11-state ODE FIM: ship / reuse vignettes/data/*.RDS (same pattern as Ex01–03).
filePathEvalPop = file.path( paths$data, "vignette4_evaluation_populationFIM.RDS" )
filePathEvalPopQuad = file.path(
  paths$data, "vignette4_evaluation_populationFIM_quadratic.RDS"
)
# Report only when regenerating the FIM (or PFIM_VIGNETTE_REPORT=true).
needReport = !file.exists( filePathEvalPop ) ||
  identical( Sys.getenv( "PFIM_VIGNETTE_REPORT" ), "true" )

# Linear FD (central differences) matches PopED's stencil. PFIM's default is
# quadratic; that path is optional here because the 11-state grid is expensive.
# Set PFIM_GCSF_FORCE_RUN=true to ignore the RDS (tests). Set
# PFIM_GCSF_FD_BOTH=true (or ship the quadratic RDS) to fill pfim_quad columns.
.runGcsfPopEval = function( linearOnly ) {
  prev_fd = pfim_get_option( "perf.fdLinearOnly" )
  pfim_set_option( perf.fdLinearOnly = linearOnly )
  on.exit( pfim_set_option( perf.fdLinearOnly = prev_fd ), add = TRUE )
  run( Evaluation(
    name = "gcsf_poped_design1",
    modelEquations = modelEquations,
    modelParameters = modelParameters,
    modelError = modelError,
    outputs = outputs,
    designs = list( design1 ),
    fimType = "population",
    odeSolverParameters = list( atol = 1e-8, rtol = 1e-8 )
  ) )
}

.gcsfLoadOrRun = function( rds_path, linearOnly ) {
  runner = function() .runGcsfPopEval( linearOnly )
  if ( identical( Sys.getenv( "PFIM_GCSF_FORCE_RUN" ), "true" ) )
    return( runner() )
  if ( file.exists( rds_path ) )
    return( readRDS( rds_path ) )
  runner()
}

# Absolute relative difference (%) vs the PopED reference. Never negative
# (signed (PFIM − PopED)/PopED can print as "<0" after rounding).
.gcsfRelPct = function( pfim, ref ) {
  rel = 100 * abs( as.numeric( pfim ) - as.numeric( ref ) ) / as.numeric( ref )
  rel[ !is.finite( rel ) ] = NA_real_
  pmax( rel, 0 )
}

.gcsfFmtRse = function( x, digits = 3 ) {
  x = as.numeric( x )
  out = sprintf( paste0( "%.", digits, "f" ), x )
  out[ !is.finite( x ) ] = "—"
  out
}

.gcsfFmtRel = function( x ) {
  x = pmax( 0, abs( as.numeric( x ) ) )
  out = sprintf( "%.3f", x )
  out[ is.finite( x ) & x < 5e-4 ] = "0.000"
  out[ !is.finite( x ) ] = "—"
  out
}

.gcsfKbl = function( df, caption ) {
  kableExtra::kbl(
    df,
    escape = FALSE,
    caption = caption,
    align = c( "l", rep( "r", max( ncol( df ) - 1L, 0L ) ) ),
    row.names = FALSE
  ) |>
    kableExtra::kable_styling(
      bootstrap_options = c( "striped", "hover", "condensed" ),
      full_width = FALSE,
      position = "center",
      font_size = 12
    )
}

evaluationPop = .gcsfLoadOrRun( filePathEvalPop, TRUE )
evaluationPopQuad = NULL
if ( identical( Sys.getenv( "PFIM_GCSF_FD_BOTH" ), "true" ) ||
     file.exists( filePathEvalPopQuad ) )
  evaluationPopQuad = .gcsfLoadOrRun( filePathEvalPopQuad, FALSE )

rse = getRSE( evaluationPop )

goldPath = system.file( "fixtures", "gold-fim.R", package = "PFIM" )
if ( !nzchar( goldPath ) || !file.exists( goldPath ) )
  goldPath = file.path( "inst", "fixtures", "gold-fim.R" )
sys.source( goldPath, envir = environment() )
gcsfGold = .pfimGold$gcsf_design1
poped_rse_mu = gcsfGold$rse_mu
poped_rse_d = gcsfGold$rse_d
# Residual *variances* in the gold table; PFIM ModelError slots = residual *SDs*.
poped_rse_sigma = gcsfGold$rse_sigma_var
poped_var_sigma = gcsfGold$var_sigma

.pfimRseRow = function( rse, prefix, name ) {
  rn = paste0( prefix, name )
  if ( rn %in% rownames( rse ) ) as.numeric( rse[ rn, "RSE" ] ) else NA_real_
}

.pfimValRow = function( rse, prefix, name ) {
  rn = paste0( prefix, name )
  if ( rn %in% rownames( rse ) ) as.numeric( rse[ rn, "parametersValues" ] ) else NA_real_
}

.gcsfCmpTables = function( rse ) {
  cmp_mu = data.frame(
    param = names( poped_rse_mu ),
    pfim = vapply(
      names( poped_rse_mu ),
      function( nm ) .pfimRseRow( rse, "\u03bc_", nm ),
      numeric( 1 )
    ),
    poped = as.numeric( poped_rse_mu ),
    row.names = NULL
  )
  cmp_mu$rel_pct = .gcsfRelPct( cmp_mu$pfim, cmp_mu$poped )

  cmp_d = data.frame(
    param = names( poped_rse_d ),
    pfim = vapply(
      names( poped_rse_d ),
      function( nm ) .pfimRseRow( rse, "\u03c9\u00b2_", nm ),
      numeric( 1 )
    ),
    poped = as.numeric( poped_rse_d ),
    row.names = NULL
  )
  cmp_d$rel_pct = .gcsfRelPct( cmp_d$pfim, cmp_d$poped )

  cmp_sigma = data.frame(
    param = names( poped_rse_sigma ),
    pfim_sd = vapply(
      names( poped_rse_sigma ),
      function( nm ) .pfimValRow( rse, "\u03c3_", nm ),
      numeric( 1 )
    ),
    poped_var = as.numeric( poped_var_sigma ),
    pfim_rse_sd = vapply(
      names( poped_rse_sigma ),
      function( nm ) .pfimRseRow( rse, "\u03c3_", nm ),
      numeric( 1 )
    ),
    poped_rse_var = as.numeric( poped_rse_sigma ),
    row.names = NULL
  )
  # Delta-method bridge: RSE(variance) ≈ 2 * RSE(SD) when estimating σ² = σ_SD².
  cmp_sigma$pfim_rse_var_equiv = 2 * cmp_sigma$pfim_rse_sd
  cmp_sigma$rel_pct = .gcsfRelPct(
    cmp_sigma$pfim_rse_var_equiv, cmp_sigma$poped_rse_var
  )
  list( cmp_mu = cmp_mu, cmp_d = cmp_d, cmp_sigma = cmp_sigma )
}

.gcsfAttachQuad = function( cmp, pfim_quad, poped ) {
  cmp$pfim_quad = pfim_quad
  cmp$rel_pct_quad = .gcsfRelPct( pfim_quad, poped )
  cmp
}

cmps = .gcsfCmpTables( rse )
cmp_mu = cmps$cmp_mu
cmp_d = cmps$cmp_d
cmp_sigma = cmps$cmp_sigma
if ( !is.null( evaluationPopQuad ) ) {
  cmps_q = .gcsfCmpTables( getRSE( evaluationPopQuad ) )
  cmp_mu = .gcsfAttachQuad( cmp_mu, cmps_q$cmp_mu$pfim, cmp_mu$poped )
  cmp_d = .gcsfAttachQuad( cmp_d, cmps_q$cmp_d$pfim, cmp_d$poped )
  cmp_sigma = .gcsfAttachQuad(
    cmp_sigma, cmps_q$cmp_sigma$pfim_rse_var_equiv, cmp_sigma$poped_rse_var
  )
}

# Tests set PFIM_GCSF_CMP_ONLY=true to skip plots / Report() (~40s).
if ( !identical( Sys.getenv( "PFIM_GCSF_CMP_ONLY" ), "true" ) ) {

showOutputEvaluationPop = pfimCapture( {
  show( evaluationPop )
  getFisherMatrix( evaluationPop )
  getSE( evaluationPop )
  getRSE( evaluationPop )
  getDeterminant( evaluationPop )
  getDcriterion( evaluationPop )
} )
writeLines( showOutputEvaluationPop,
            file.path( paths$outputs, "vignette4_evaluation_populationFIM_show.txt" ) )

# ---- Prediction overlay (3 doses, PK | ANC) ---------------------------------
# True ODE on a dense 0..tmax grid (PopED model_num_points). Never interpolate
# FIM sampling times: Design 1 is dense on days 1 and 7 only, so a line through
# the FIM cache skips the five intermediate daily peaks.

.pfimGcsfPopedStylePlot = function( evaluation, n_points = 400L ) {
  model0 = PFIM:::rebuildEvalModel( evaluation, finiteDifference = FALSE )
  design = prop( evaluation, "evaluationDesign" )[[ 1L ]]
  if ( is.null( design ) )
    design = prop( evaluation, "designs" )[[ 1L ]]
  arms = prop( design, "evaluationArms" )
  if ( !length( arms ) )
    arms = prop( design, "arms" )

  dose_ug_kg = c( dose_1ugkg = 1, dose_3ugkg = 3, dose_10ugkg = 10 )
  out_map = c( RespPK = "Model: PK", RespPD = "Model: ANC" )

  samp_by_arm = lapply( arms, function( arm ) {
    st = prop( arm, "samplingTimes" )
    stats::setNames(
      lapply( st, function( s ) prop( s, "samplings" ) ),
      vapply( st, function( s ) prop( s, "outcome" ), character( 1L ) )
    )
  } )

  tmax = max( unlist( samp_by_arm, use.names = FALSE ), 0 )
  dose_t = seq( 0, 6 * 24, by = 24 )
  dense = sort( unique( c(
    seq( 0, tmax, length.out = n_points ),
    unlist( samp_by_arm, use.names = FALSE ),
    dose_t,
    as.vector( outer( dose_t, c( 0.05, 0.167, 0.25, 0.5, 1, 2, 4, 8 ), `+` ) )
  ) ) )
  dense = dense[ dense >= 0 & dense <= tmax + 1e-9 ]

  pieces = Map( function( arm, samp_by ) {
    armName = prop( arm, "name" )
    dose = unname( dose_ug_kg[[ armName ]] )
    group = paste0( "Group ", match( dose, c( 1, 3, 10 ) ) )

    plot_arm = PFIM:::.pfimCloneS7( arm )
    prop( plot_arm, "evaluationModel" ) = list()
    prop( plot_arm, "samplingTimes" ) = list(
      SamplingTimes( outcome = "RespPK", samplings = dense ),
      SamplingTimes( outcome = "RespPD", samplings = dense )
    )

    m = PFIM:::defineModelAdministration( PFIM:::.pfimCloneS7( model0 ), plot_arm )
    y = PFIM:::evaluateModel( m, plot_arm )

    lapply( names( out_map ), function( out ) {
      df = y[[ out ]]
      if ( is.null( df ) || !nrow( df ) || !out %in% names( df ) ) return( NULL )
      samp = samp_by[[ out ]]
      pred_at = vapply( samp, function( tt ) {
        df[[ out ]][ which.min( abs( df$time - tt ) ) ]
      }, numeric( 1L ) )
      list(
        curve = data.frame(
          time = df$time, pred = df[[ out ]],
          Group = group, Model = unname( out_map[[ out ]] ),
          stringsAsFactors = FALSE
        ),
        pts = data.frame(
          time = samp, pred = pred_at,
          Group = group, Model = unname( out_map[[ out ]] ),
          stringsAsFactors = FALSE
        )
      )
    } )
  }, arms, samp_by_arm )

  flat = unlist( pieces, recursive = FALSE )
  flat = flat[ !vapply( flat, is.null, logical( 1L ) ) ]
  curves = do.call( rbind, lapply( flat, `[[`, "curve" ) )
  pts    = do.call( rbind, lapply( flat, `[[`, "pts" ) )
  curves$Model = factor( curves$Model, levels = c( "Model: PK", "Model: ANC" ) )
  curves$Group = factor( curves$Group, levels = paste0( "Group ", 1:3 ) )
  pts$Model = factor( pts$Model, levels = levels( curves$Model ) )
  pts$Group = factor( pts$Group, levels = levels( curves$Group ) )

  cols = c( "Group 1" = "#E41A1C", "Group 2" = "#4DAF4A", "Group 3" = "#377EB8" )
  ggplot2::ggplot( curves, ggplot2::aes( time, pred, color = Group ) ) +
    ggplot2::geom_line( linewidth = 0.7 ) +
    ggplot2::geom_point( data = pts, size = 1.5, alpha = 0.85 ) +
    ggplot2::facet_wrap( ~ Model, scales = "free_y", nrow = 1L ) +
    ggplot2::scale_color_manual( values = cols ) +
    ggplot2::labs(
      x = "Time from first dose (hours)",
      y = "Model Predictions",
      color = "Group"
    ) +
    ggplot2::theme_bw( base_size = 12 ) +
    ggplot2::theme(
      legend.position = "right",
      strip.background = ggplot2::element_rect( fill = "grey92", color = NA ),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plotEval_SE = pfimSafePlot( PFIM::plotSE( evaluationPop ) )
plotEval_RSE = pfimSafePlot( PFIM::plotRSE( evaluationPop ) )
plotOutcomesEvaluationPopedStyle = pfimSafePlot(
  .pfimGcsfPopedStylePlot( evaluationPop )
)
plotsEval = pfimSafePlot( plotEvaluation( evaluationPop, plotOptions ) )
plotOutcomesEvaluationRespPK = if ( !is.null( plotsEval ) )
  plotsEval[[ "gcsf_design1" ]][[ "dose_10ugkg" ]][[ "RespPK" ]] else NULL
plotOutcomesEvaluationRespPD = if ( !is.null( plotsEval ) )
  plotsEval[[ "gcsf_design1" ]][[ "dose_10ugkg" ]][[ "RespPD" ]] else NULL

# HTML Report() rebuild is ~40s here; only regenerate with the FIM cache
# (NOT_CRAN=true + delete the RDS), or set PFIM_VIGNETTE_REPORT=true.
if ( isTRUE( needReport ) ) {
  invisible( pfimSafeReport(
    evaluationPop, paths$reports, "vignette4_evaluation_popFIM_report.html",
    plotOptions
  ) )
}

}  # PFIM_GCSF_CMP_ONLY

invisible( NULL )
