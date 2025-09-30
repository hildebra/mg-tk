library(rmarkdown)

# Usage:
# Rscript QC_make_html.R <Rscripts directory> <metagStats.txt filepath> <output filepath>

# Get the command line arguments
argv <- commandArgs(TRUE)

# The working directory needs to be passed in manually
# It should be this R_scripts folder
setwd(argv[1])

# Render the document
rmarkdown::render(
    'QC_html_report.Rmd',
    params=list(statspath = argv[2]), 
    output_file = argv[3]
)
