# Frozen PFIM gold numbers.
#
# If a test against this file fails, the implementation changed.
# Do not edit these values to silence a failure. Investigate the code.

.pfimGold <- list(
  gcsf_design1 = list(
    rse_mu = c(
      KA = 3.289140, KEL = 13.518742, VD = 12.604107, KD = 8.265412,
      KINT = 3.822461, KSI = 19.528269, KMT = 4.066511, KTT = 17.017908,
      NB0 = 10.771361, SC1 = 20.851250, SM1 = 8.289988, SM2 = 7.664966
    ),
    rse_d = c(
      NB0 = 26.961113, KEL = 32.644613, VD = 31.369081,
      KSI = 41.941250, SC1 = 29.156946, SM1 = 26.406506
    ),
    rse_sigma_var = c(
      slope_RespPK = 3.900224,
      slope_RespPD = 6.541108,
      inter_RespPD = 7.713365
    ),
    var_sigma = c(
      slope_RespPK = 2.53e-01,
      slope_RespPD = 2.27e-02,
      inter_RespPD = 2.10e+00
    )
  ),
  cas1_noCov_noIOV = list(
    D = 731.067615428391
  ),
  cas2_noCov_withIOV = list(
    D   = 2048.3263557128,
    M11 = 342.236485,
    M88 = 19961.4704,
    M1010 = 36491.07585
  ),
  cas3_withCov_noIOV = list(
    D   = 576.120006672508,
    M11 = 345.593647,
    M55 = 1493.4589259,
    M88 = 17000.86746
  ),
  cas4_withCov_withIOV = list(
    D = 1568.77419074749
  ),
  cas5_withCov_mixed = list(
    D = 1457.36661504908
  ),
  cas6_occasion_iov = list(
    D     = 1837.7481985995,
    M11   = 334.061402,
    M1010 = 37340.564
  )
)
