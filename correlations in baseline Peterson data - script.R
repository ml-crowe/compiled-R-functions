# The below code assumes your dataset is in CSV format, with
# each column corresponding to one variable, and each row
# corresponding to one participant. To avoid issues, it's best to
# remove all columns other than those corresponding to individual
# items from the several psychological measures we're interested in.
# That is, it's best to remove subject ID variable, gender and any other demographics, etc.
#
#
# This script reads that CSV into R as a dataframe, and then uses the
# corr_cross() function from the {lares} package to depict the strongest
# bivariate correlations in the dataset. We are interested in correlations
# higher than perhaps r=.5. If there are tons of those, we can just focus on ones higher than .7.
#
# The code may take a while to run - in my experience this package isn't very optimized.


install.packages('lares', dependencies = TRUE) # Install lares package

library(lares) # load the lares package

items_dataset <- read.csv( # read your CSV file into an R dataframe
  'PASTE YOUR CSV FILE PATH HERE INSIDE THE QUOTES. CHANGE BACKSLASHES \ TO FORWARD SLASHES / IF USING WINDOWS'
)
    # (In Windows, filepaths by default contain backslashes. Change
    #   them to forward slashes for R to read it correctly.
    # Make sure you have read the initial note above; this dataset should
    #   contain only item-level variables, no subject ID, no demographics, etc.)

corrplot1 <- lares::corr_cross(items_dataset, # create ggplot object in R, to use in next steps
                       top = 50,
                       rm.na = TRUE,
                       method = "spearman")
corrplot1 # put the plot itself into Rstudio's 'Plots' pane

getwd() # print out the path to folder where the plot will save as PDF (in the next step)

ggplot2::ggsave('highest-item-correlations.pdf', # this saves the plot as a PDF file -- feel free to change its name
                plot = corrplot1,
                device = 'pdf',
                dpi = 'retina',
                width = 6.5,
                height = 9)
  # If the plotting or the saving of the plot doesn't work, run the 'tidyverse' lines below and
  # return to the 'ggsave' command to attempt the PDF file save again.
  # Also, instead of just specifying the filename you wish to save to, you can also specify an entire path,
  # in case the default working directory -- printed with getwd() -- is not a convenient save location

# The code below will:
# Get full, ranked correlation list into the Rstudio Viewer, ------------------------
# so you can save it as a standalone HTML webpage file
# (to be shared with the group)

# install.packages('tidyverse') # !!! UNCOMMENT this line if you need to first install the tidyverse packages!
install.packages('correlation')
install.packages('gt')
library('tidyverse') # run install.packages('tidyverse') line above if you don't already have it installed
library('correlation')
library('gt')

correlation::correlation(items_dataset, method = 'spearman') %>% # makes a massive dataframe of all correlations
  dplyr::arrange(-rho) %>% # ...and then arranges them by effect size rho,
  filter(p < 0.1) %>% # filters out correlations whose p-value is greater than .1
  select(-S) %>%
  gt::gt() %>%  # ...and then converts that dataframe into a table in Rstudio's Viewer...
  gt::fmt_number(3:7, decimals = 4) #  ...which can then be saved as a standalone webpage (HTML file)...instructions below


# To save the entire correlation table from the Viewer to a standalone HTML
# file that you can share, in the Viewer click 'Export' -> 'Save as Web page'
