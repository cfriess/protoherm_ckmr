

## ---------------------------
##
## Script name: compile_pairs_quickin.R
##
## Purpose of script: return numbers of different kin types from the indivs, 
## pairs and samps output from the FishSims
##
## Author: Eric Anderson
##
##
##
## ---------------------------


#' From the indivs, pairs and samps output from the FishSims, return numbers of different kin types
#' 
#' This was a total rewrite for using the output of `quickin()`.
#' @param samps Either a data frame of the samples output from FishSim, or, if a
#' string, it is assumed to be a path to an rds file holding that output.
#' @param quick_out Either a list that is the output from `quickin()`, or, if a
#' string, it is assumed to be a path to an rds file holding that output.
#' This has to be output from the quickin-for-sex-changers branch of Eric's
#' fork of the fishSim repo.  That is taken care of with renv.
#' @param nll_type binary variable specifying whether we are including HSPs only (0) or POPs as well (1)
#' @return This function returns a list of two tibbles. The `regular` one holds columns:
#'     - c1: birth year (cohort) of the first-born of the pair
#'     - c2: birth year (cohort) of second-born of the pair
#'     - s1: the sampling year of the first-born
#'     - g1: sex of the first-born (vector of 0 if not including POPs)
#'     - PO: number of father-offspring pairs
#'     - MO: number of mother-offspring pairs
#'     - HSD: number of half-sibling pairs with different mtDNA (no shared mother)
#'     - HSS: number of half-sibling pairs with the same mtDNA (shared mother)
#'     - Not: number of pairs that are not any of the above relationships
#'
#' The `extended` one replaces HSD and HSS with:
#'     - HSFF: number of half-sibling pairs with shared parent mother of first and second born
#'     - HSFM: number of half-sibling pairs with shared parent mother of first and father of second
#'     - HSMM: number of half-sibling pairs with shared parent father of first and second born.
#' Note that this returns the number of pairs that are from the same cohort, as
#' well.  Before use for CKMR, you will probably want to filter those out
#' so as to use only cross-cohort pairs.
#' @example
#' compile_pairs(
#'   samps = "data/gag_samps_Ninit320K_20221013.rds",
#'   quick_out = "data/gag_quinres_Ninit320K_20221013.rds",
#'   nll_type = 0
#' )
compile_pairs_quickin <- function(samps, quick_out, nll_type) {
  
  #### Do this for dependencies ####
  require(tidyverse)
  require(fastmatch)
  require(fishSim)
  
  #### deal with the different possible input types ####
  if(is.data.frame(samps)) {
    samps_tib <- as_tibble(samps)
  } else {
    samps_tib <- as_tibble(read_rds(samps))
  }
  
  if(is.list(quick_out)) {
    quick_list <- quick_out
  } else {
    quick_list <- read_rds(quick_out)
  }
  
  
  #### Some functions that will be useful  ####
  
  #' reorder quickin output table so that the first column
  #' is always the first born (or the lexicographically earlier
  #' one if the are from the same cohort)
  #' @param q the data frame from quickin()
  #' @param s the tibble of samples
  pair_id_reord <- function(q, s) {
    q %>%
      as_tibble() %>%
      left_join(s %>% select(Me, BirthY, SampY, Sex) %>% rename(by1 = BirthY, dy1 = SampY, g1 = Sex), by = c("keeps1" = "Me")) %>%
      left_join(s %>% select(Me, BirthY, SampY, Sex) %>% rename(by2 = BirthY, dy2 = SampY, g2 = Sex), by = c("keeps2" = "Me")) %>%
      mutate(
        indiv1 = case_when(
          by1 < by2 ~ keeps1,
          by1 == by2 & keeps1 < keeps2 ~ keeps1,
          TRUE ~ keeps2
        ),
        indiv2 = case_when(
          by1 < by2 ~ keeps2,
          by1 == by2 & keeps1 < keeps2 ~ keeps2,
          TRUE ~ keeps1
        ),
        c1 = case_when(
          by1 < by2 ~ by1,
          by1 == by2 & keeps1 < keeps2 ~ by1,
          TRUE ~ by2
        ),
        c2 = case_when(
          by1 < by2 ~ by2,
          by1 == by2 & keeps1 < keeps2 ~ by2,
          TRUE ~ by1
        ),
        s1 = case_when(
          by1 < by2 ~ dy1,
          by1 == by2 & keeps1 < keeps2 ~ dy1,
          TRUE ~ dy2
        ),
        g1 = case_when(
          by1 < by2 ~ g1,
          by1 == by2 & keeps1 < keeps2 ~ g1,
          TRUE ~ g2
        )
      ) %>%
      select(indiv1, indiv2, c1, c2, s1, g1)
  }
  
  #' re-order the tibble and then add mom and dad on there
  process_quickin <- function(q, s) {
    if(is.null(q)) {
      return(
        tibble()
      )
    }
    pair_id_reord(q, s) %>%
      left_join(s %>% select(Me, Dad, Mum) %>% rename(pa1 = Dad, ma1 = Mum), by = c("indiv1" = "Me")) %>%
      left_join(s %>% select(Me, Dad, Mum) %>% rename(pa2 = Dad, ma2 = Mum), by = c("indiv2" = "Me")) %>%
      select(indiv1:c2, s1, g1, pa1, pa2, ma1, ma2)
  } 
  
  #### Tackle the half-siblings, assigning them to be HSFF or HSMM
  # note that this is mostly to check that we get the same results
  # as the quickin function, as it returns maternal and paternal
  # half-sibs
  hs_tib <- quick_list$HSP %>%
    process_quickin(samps_tib) %>%
    mutate(
      relat = case_when(
        pa1 == pa2 ~ "HSMM",
        ma1 == ma2 ~ "HSFF",
        ma1 == pa2 ~ "HSFM-Unexpected",  # it turns out that the quickin doesn't get these in HSP
        TRUE ~ "Weird! Something is amiss!"
      )
    )
  
  #### Now we get the HSFMs out of there.  I hacked quickin() to
  # return those
  hsfm_tib <- quick_list$qk_HSFM %>%
    process_quickin(samps_tib) %>%
    mutate(relat = "HSFM")
  
  # Now I just have to get the POs and MOs and bind_rows
  # and then add the Total number of pairs from the sample
  # information, and subtract out the true pairs to get the Nots,
  # and bind_rows those on there, then pivot wider.
  
  POs <- quick_list$POP %>% 
    process_quickin(samps_tib) %>%
    mutate(
      relat = case_when(
        indiv1 == ma2 ~ "MO",
        indiv1 == pa2 ~ "PO",
        TRUE ~ "Weird Schtuff with the Parents!"
      )
    )
  
  # put them together and count
  cr_counts <- bind_rows(
    hs_tib,
    hsfm_tib,
    POs
  ) %>%
    count(c1, c2, s1, g1, relat)
  
  # pivot just the related ones:
  relat_counts <- cr_counts %>% 
    pivot_wider(names_from = relat, values_from = n, values_fill = 0L)
  
  tot_relat_counts <- cr_counts %>% 
    group_by(c1, c2, s1, g1) %>%
    summarise(tot_relats = sum(n)) %>%
    ungroup()
  
  # now, count up the number of samples of different cohorts, and use that
  # to get the total number of pairs
  # need to include older individuals' sampling year and sex in order to do POPs
  cc <- samps_tib %>%
    count(BirthY, SampY, Sex) %>%
    rename(c = BirthY, s = SampY, g = Sex) %>%
    mutate(csg = paste(c,s,g, sep = "_"))
  
  # now, get all the possible pair categories and add up
  # compute the number of possible pairs in each
  tot_pairs <- expand_grid(csg1 = cc$csg, csg2 = cc$csg) %>%
    left_join(cc, by = c("csg1" = "csg")) %>%
    rename(n1 = n,
           s1 = s,
           c1 = c,
           g1 = g) %>%
    left_join(cc, by = c("csg2" = "csg")) %>%
    rename(n2 = n,
           s2 = s,
           c2 = c,
           g2 = g) %>%
    filter(c1 <= c2,
           !(c1 == c2 & s1 > s2), # fix duplicates
           !(c1 == c2 & s1 == s2 & g1 == "F" & g2 == "M") # fix duplicates
    ) %>%
    mutate(
      npairs = case_when(
        csg1 == csg2 ~ n1 * (n1 - 1) / 2,
        TRUE ~ n1 * n2 * 1.0
      )
    ) %>%
    filter(npairs > 0) %>%
    group_by(c1, c2, s1, g1) %>%
    summarize(npairs = sum(npairs))
  
  # now, join the total number of relats on there, and from
  # that, compute the number of Not-related pairs
  Nots <- tot_pairs %>%
    left_join(tot_relat_counts, by = c("c1", "c2", "s1", "g1")) %>%
    replace_na(list(tot_relats = 0L)) %>%
    mutate(Not =  npairs - tot_relats) %>%
    select(c1, c2, s1, g1, Not)
  
  # finally, join relat counts onto the Nots and
  # put the columns in the order that we said we would. We have to do some
  # fancy things in case any of columns PO, MO, HSFF, HSFM, HSMM are 0.
  tmp <- Nots %>%
    left_join(relat_counts, by = c("c1", "c2","s1","g1"))
  for(i in setdiff(c("PO", "MO", "HSFF", "HSFM", "HSMM"), names(tmp))) {
    tmp[[i]] <- NA_integer_
  }
  extended <- tmp %>%
    replace_na(list(PO = 0L, MO = 0L, HSMM = 0L, HSFF = 0L, HSFM = 0L)) %>%
    select(c1, c2, s1, g1, PO, MO, HSFF, HSFM, HSMM, Not) %>%
    mutate(g1 = ifelse(g1 == "F", 0, 1)) # recode females as 0 and males as 1
  
  if(!nll_type){
    extended <- extended %>%
      group_by(c1, c2, s1) %>%
      summarize(PO = sum(PO),
                MO = sum(MO),
                HSFF = sum(HSFF),
                HSFM = sum(HSFM),
                HSMM = sum(HSMM),
                Not = sum(Not)) %>%
      mutate(g1 = 9, .after = s1) # dummy variable so we can read the same data into the tpl file
  }
  
  regular <- extended %>%
    mutate(
      HSD = HSFM + HSMM,
      HSS = HSFF,
      .after = MO
    ) %>%
    select(-HSFF, -HSFM, -HSMM)
  
  list(
    regular = regular,
    extended = extended
  )
  
}
