#' Pneumococcal Conjugate Vaccine Trial
#'
#' The Pneumococcal Conjugate Vaccine Trial was a cluster randomized trial
#' carried out from 1997 to 2000 to assess the safety and efficacy of a
#' seven-valent conjugate pneumococcal vaccine (O'Brien et al., 2003). The
#' study population was Navajo and White Mountain Apache children younger
#' than 2 years, a group with one of the highest documented rates of invasive
#' pneumococcal disease in the world at that time. A total of 38 geographic
#' areas were randomized: 19 areas were offered pneumococcal vaccine and 19
#' were offered a comparator (meningococcal vaccine). These data, made available
#' by Hayes and Moulton (2017) on Harvard Dataverse, include a
#' random subsample of 449 children drawn from the original 8,292 trial
#' participants.
#'
#' @format A data frame with 449 observations and 4 variables:
#' \describe{
#'     \item{randunit}{numeric, distinct for each randomised geographic area}
#'     \item{bpepisodes}{number of bacterial pneumonia episodes}
#'     \item{spnvac}{test vaccine: 0=meningococcal comparator, 1=pneumococcal vaccine}
#'     \item{fakeid}{recoded individual participant id}
#' }
#' @source \url{https://dataverse.harvard.edu/dataverse/crt}
#' @references
#' Hayes, R. J. and Moulton, L. H. (2017). Cluster Randomised Trials 2nd
#' edition. New York: Chapman and Hall/CRC.
#'
#' O'Brien, K. L. et al. (2003). Efficacy and safety of seven-valent conjugate
#' pneumococcal vaccine in American Indian children: group randomised trial.
#' Lancet 362, 355--361.
"pneumovac"

#' The Botswana Combination Prevention Project
#'
#' The Botswana Combination Prevention Project (BCPP) was a pair-matched HIV
#' prevention CRT to test whether a combination treatment and prevention
#' intervention could reduce population-level cumulative HIV incidence over 3
#' years of follow-up. A total of 30 communities were randomized: 15 to the
#' intervention arm (combination treatment and prevention package) and 15 to
#' the control arm (enhanced standard of care). The primary study endpoint was
#' cumulative HIV incidence, measured at scheduled study visits as time to
#' HIV-infection within a cohort of individuals identified as HIV-negative among
#' a 20\% random sample of eligible households at baseline. That is, we have an
#' interval-censored time-to-event outcome for each cohort participant. Since
#' the primary trial results have not yet been published, this is a simulated
#' data set, which was generated to mimic the BCPP by applying an agent-based
#' epidemic model to a dynamic network of simulated sexual partnerships
#' (Goyal, Blitzstein, and De Gruttola, 2013; Wang et al., 2014).
#'
#' @format A data frame with 10465 observations and 5 variables:
#' \describe{
#'     \item{group}{numeric, distinct for each randomised community}
#'     \item{pair.id}{numeric, distinct for each matched pair}
#'     \item{treat}{indicates assignment to the intervention arm}
#'     \item{left}{left bound of interval censored weeks to HIV-infection}
#'     \item{right}{right bound of interval censored weeks to HIV-infection}
#' }
#'
#' @references
#' Goyal, R., Blitzstein, J., and De Gruttola, V. (2013). Simulating Bipartite
#' Networks to Reflect Uncertainty in Local Network Properties. Harvard
#' University Biostatistics Working Paper Series.
#'
#' Wang, R., Goyal, R., Lei, Q., Essex, M., and De~Gruttola, V. (2014). Sample
#' size considerations in the design of cluster randomized trials of combination
#' HIV prevention. Clinical Trials 11, 309–318.
#'
"bcpp"
