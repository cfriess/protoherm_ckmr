

# simple awk script to parse the snakemake output paths
# in the various rules and return each in a list.  For each
# path we have:
#  - glob:  a string suitable for finding all such files using Sys.glob
#  - regex: a string suitable for parsing the path using tidyr::extract()
#  - into: a string giving the column names for tidyr::extract()
#
# To use it, just run it on the Snakefile like:
#
#   sed 's/^\t*//g; s/^ *//g;' Snakefile | awk -f scripts/snake_paths_to_list.awk > R/snake_path_list.R
#
# It produces R code to make a list called snake_paths.  After this, doing
# source("R/snake_path_list")  
# will load that list into R where it can be used to harvest the outputs

# this is a pretty naive script that assumes that Snakefile is in
# the format that ours is in

BEGIN {
  curr_block = "NULL";
  print "snake_paths = list()"
}


/^rule/ {
  current_rule = $2; 
  gsub(/[: ]/, "", current_rule)
  if(current_rule != "all") {
    print "snake_paths$" current_rule, "= list()"
    print "snake_paths$" current_rule "$output = list()"
    print "snake_paths$" current_rule "$log = list()"
    print "snake_paths$" current_rule "$benchmark = list()"
  }
}

$1~/:/ {curr_block = "NULL"}

$1=="output:" {curr_block = "output";}
$1=="log:" {curr_block = "log"; }
$1=="benchmark:" {curr_block = "benchmark";}

curr_block != "NULL" && (/=/ || /results\/benchmarks/) {
  if(match($0, /=/)) {
    n = split($1, a, /=/)
    name = a[1]
    path = a[2]
  } else if(match($0, /benchmarks/)) {
    name = "bm"
    path = $1
  }
  gsub( /[",]/, "", path) # get quotes and commas out of it
  
  # make the list to hold the regex, into, and glob
  print "snake_paths$" current_rule "$" curr_block "$" name " = list()"
  
  # split path into a vector for handling
  nr = split(path, ar, /[{}]/)
  
  # now, we parse path into regex
  regex = "\""
  for(i=1;i<=nr;i++) {
    if(i%2==1) regex = sprintf("%s%s",regex, ar[i]); 
    else regex = sprintf("%s(.*)", regex);
  }
  regex = sprintf("%s\"", regex)
  # print that regex into the list
  print "snake_paths$" current_rule "$" curr_block "$" name "$regex = " regex
  
  
  # now, we parse path into "into"
  into = "c("
  for(i=1;i<=nr;i++) {
    if(i%2==0) {
      if(i==2) into = sprintf("%s\"%s\"", into, ar[i]); 
      else into = sprintf("%s, \"%s\"", into, ar[i]);
    }
  }
  into = sprintf("%s)", into)
  # print that "into" into the list
  print "snake_paths$" current_rule "$" curr_block "$" name "$into = " into
  
  
  # now, we parse path into glob
  glob = "\""
  for(i=1;i<=nr;i++) {
    if(i%2==1) glob = sprintf("%s%s",glob, ar[i]); 
    else glob = sprintf("%s*", glob);
  }
  glob = sprintf("%s\"", glob)
  # print that regex into the list
  print "snake_paths$" current_rule "$" curr_block "$" name "$glob = " glob
  
}

