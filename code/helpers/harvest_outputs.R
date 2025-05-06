

## ---------------------------
##
## Script name: harvest_outputs.R
##
## Purpose of script: use the `snake_paths` object obtained by sourcing 
## R/snake_path_list to extract simulation outputs
##
## Author: Eric Anderson
##
##
##
## ---------------------------


require(tidyverse)
#' use the `snake_paths` object obtained by sourcing R/snake_path_list to extract outputs
#' 
#' Retrieves everything into a tibble, possibly with list columns.
#' @param rule  the name of the rule that produced the output
#' @param block either output, log, or benchmark
#' @param key  the snakemake file identifier, like plot_sheet, or rds, or log
#' @param action this is one of: 
#'  - read_rds    --- read contents of an RDS file into a list column
#'  - read_tsv    --- read a tsv file (that must all have the same headers) into a
#'                    list column, but then unnest that column into a bunch of new columns.
#'                    (this is primarily for reading the benchmarks)
#'  - read_subset --- reads the rds in and assumes it is a list, and then just pick out the
#'                    components named in `subs`.
#'  - read_lines  --- read text lines into a vector and store in a list column.  This
#'                    is mostly for reading log files.
#'  - store_path  --- just return a tibble that has all the parameter settings and the path
#'                    to the file.  Useful for collection pdf paths, etc.
#' @examples 
#' # get the key_objs outputs of the estimate rule
#' key_objs <- harvest_outputs(rule = "estimate", block = "output", key = "key_obs", action = "read_rds")
#' 
#' # get a subset of the key_pbjs output.  This is not necessary any longer, and is deprecated, but wanted
#' # to leave it here so we know how to pull subsets from things.
#' key_objs <- harvest_outputs(rule = "estimate", block = "output", key = "key_obs", action = "read_subset", subs = c("abund_res", "para_res", "sim_pop", "ckmr_obs", "ckmr_true"))

#' 
#' # get the rds outputs of the estimate rule
#' est_rds <- harvest_outputs(rule = "estimate", block = "output", key = "rds", action = "read_rds")
#'
#' # get the benchmarks (time and memory) of the make_pedigree rule
#' make_ped_benchmarks <- harvest_outputs(rule = "make_pedigree", block = "benchmark", key = "bm", action = "read_tsv")
#'
#' # get the paths to the plot_sheet output of estimate
#' est_plot_sheet <- harvest_outputs(rule = "estimate", block = "output", key = "plot_sheet", action = "store_path")
#'
#' sample_pop_log <- harvest_outputs(rule = "sample_pop", block = "log", key = "log", action = "read_lines")
#' 
#' sample_pop_benchmarks <- harvest_outputs(rule = "sample_pop", block = "benchmark", key = "bm", action = "read_tsv")
harvest_outputs <- function(
    rule,
    block,
    key,
    action,
    subs = TRUE
) {
  source("R/snake_path_list.R")
  
  # let x be the relevant part of the snake_paths list
  x <- snake_paths[[rule]][[block]][[key]]
  
  files <- Sys.glob(x$glob)
  
  # get the basic tibble structure
  tib <- tibble(
    path = files
  ) %>%
    tidyr::extract(
      path,
      into = x$into,
      regex = x$regex, 
      remove = FALSE,
      convert = TRUE
    ) %>%
    mutate(
      rule = rule,
      block = block,
      key = key,
      .before = path
    )
  
  # here, we make ri_seed and cv_seed %03d again  # TURNS OUT I DIDN'T REALLY WANT TO DO THAT
  #if("ri_seed" %in% names(tib)) {
  #  tib <- tib %>%
  #    mutate(ri_seed =  sprintf("%03d", ri_seed))
  #}
  #if("cv_seed" %in% names(tib)) {
  #  tib <- tib %>%
  #    mutate(cv_seed =  sprintf("%03d", cv_seed))
  #}
  
  # now, do the action
  if (action == "store_path") {
    return(tib)
  } else if(action == "read_rds") {
    tib2 <- tib %>%
      mutate(rds = map(path, read_rds))
  } else if(action == "read_subset") {
    tib2 <- tib %>%
      mutate(rds = map(path, function(x) {y <- read_rds(x); y[subs]} ))
  } else if(action == "read_tsv") {
    tib2 <- tib %>%
      mutate(tsv = map(path, read_tsv, show_col_types = FALSE, progress = FALSE)) %>%
      unnest(tsv)
  } else if(action == "read_lines") {
    tib2 <- tib %>%
      mutate(lines = map(path, read_lines, progress = FALSE))
  } else {
    stop("unrecognized action.  Must be one of \"read_rds\", \"read_tsv\", \"read_lines\", or \"store_path\"")
  }
  tib2
}


# now, here is a function that uses that in order to get the necessary
# things to run the estimate step locally, easily---this is an updated
# version that returns a list.  One that has the sample_pop_out_list,
# and another that has the estimate snakemake object in a column.
# THIS IS DEPRECATED---NOW THAT WE DON'T HAVE A SAMPLE_POP RULE IT
# NO LONGER WORKS!
local_snake_list <- function() {
  
  sample_pop_out_list <- harvest_outputs(
    rule = "sample_pop", 
    block = "output", 
    key = "out_list", 
    action = "read_rds"
  ) %>%
    rename(in_list = rds) %>%
    select(-rule, -block, -key, -path)
  
  estimate_snake_obj <- harvest_outputs(
    rule = "estimate", 
    block = "log", 
    key = "snake_obj", 
    action = "read_rds"
  ) %>% 
    rename(snake_obj = rds) %>%
    select(-rule, -block, -key, -path)
  
  list(
    sample_pop_out_list = sample_pop_out_list,
    estimate_snake_obj = estimate_snake_obj
  )
  
}


# now that we don't have a pop_sample rule, we need to get all the ped.Rdata
# files if we want to re-run the estimate step locally.  However, that is
# 17 Gb worth of data spread across 120 Rdata files.  We don't want to
# have to read and load all those.  So, this function returns a tibble
# that has a list column with the estimate_snake_obj's, and another
# column which is the path to the ped.Rdata.  Then I will use rclone
# with --include to copy all the ped.Rdata's to a directory on google
# drive that maintains the file structure.  I should probably tarball
# them all.  
# NOTE: to run this, you should be able to just do:
# write_local_snake(date = "YYYY-MM_DD")
write_local_snake <- function(loc_snake_dir = "LOCAL_SNAKE", date = Sys.Date()) {
  
  # get the snakemake objects in a list column in a tibble
  message("storing snake objects in a tibble")
  estimate_snake_obj <- harvest_outputs(
    rule = "estimate", 
    block = "log", 
    key = "snake_obj", 
    action = "read_rds"
  ) %>% 
    rename(snake_obj = rds) %>%
    select(-rule, -block, -key, -path)
  
  # now, get a vector of the unique input paths for those ped.Rdata files
  ped_paths <- lapply(estimate_snake_obj$snake_obj, function(x) x@input$in_list) %>%
    unlist() %>%
    unique()
  
  # and write them to a temp file that can be used with tar --include
  tmp <- tempfile()
  cat(ped_paths, file = tmp, sep = "\n")
  
  # make the output directory:
  outdir <-  file.path(loc_snake_dir, date)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outtar <- file.path(outdir, "Ped-Rdata.tar")
  
  message("Compiling pedigree Rdata files into tar file ", outtar)
  # now, tar up the many gigs of ped outdata
  CALL <- paste("tar -cvf", outtar, "-T", tmp, collapse = " ")
  system(CALL)
  
  
  # write out the estimate_snake_obj
  eso_path <- file.path(outdir, "estimate_snake_obj_tibble.rds")
  message("writing snake_object tibble to ", eso_path)
  write_rds(estimate_snake_obj, file = eso_path)
  
  message("Now, do this on the Unix shell:")
  message(" ")
  message("rclone copy  --dry-run  LOCAL_SNAKE  gdrive-rclone:Gag_grouper/LOCAL_SNAKE")
  
  TRUE
}




# this is a simple function.  Once you have whittled a tibble 
# produced by local_snake_list down to a single row by filtering,
# and pass it into this function, along with the corresponding
# sample_pop_out_list tibble, it will return a LOCAL_SNAKE list.
#' @param X the tibble that is like the estimate_snake_obj component
#' of the list returned by local_snake_list, BUT, it has been filtered
#' down to have only one row (the run you are interested in).
#' @param S the sample_pop_out_list portion of the output list from 
#' local_snake_list produced when the tibble from which X comes was made.
make_local_snake <- function(X, S) {
  if(nrow(X) != 1) stop("X must be a tibble with a single row")
  
  Z <- left_join(X, S)
  
  list(
    in_list = Z$in_list[[1]],
    snake_obj = Z$snake_obj[[1]]
  )
}



# get all the plots in a tibble named by eofsr--ped_rep--sample_rep--ri_seed--cv_seed.pdf
#' @param subdir the name of the subdirectory where you want to put them inside of
#' a PLOT_OUTPUTS directory
#' Once this is done, eric gets them to google drive with:
#' cd PLOT_OUTPUTS
#' tar -cvf actual_subdir_name.tar actual_subdir_name
#' gzip actual_subdir_name.tar
#' rm -r actual_subdir_name
#' rclone copy  PLOT_OUTPUTS  gdrive-rclone:Gag_grouper/PLOT_OUTPUTS
rename_plots <- function(X, subdir) {
  # make it less unwieldy
  X2 <- X %>%
    select(eofsr, ped_rep, samp_rep, ri_seed, cv_seed, path)
  
  ddir <- file.path("PLOT_OUTPUTS", subdir)
  dir.create(ddir, recursive = TRUE, showWarnings = FALSE)
  X3 <- X2 %>%
    mutate(call = str_c("cp ", path, " ", ddir, "/", eofsr, "--", ped_rep, "--", samp_rep, "--", ri_seed, "--", cv_seed, ".pdf; "))
  
  # we cycle over these and do them individually, because if you send system
  # a string that is too big, it seems to fail.
  for(i in 1:nrow(X3)) {
    CALL <- X3$call[i]
    system(CALL)
  }
}

