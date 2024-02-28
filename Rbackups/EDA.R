# NEITHER OF THE IS PERFECT

# NOT TOO INFORMATIVE - MASKS SUMMARIZE()
# install the dataMaid library
# install.packages("dataMaid")
# load the dataMaid library
library(dataMaid) 
library(cachem)
library(bslib)

# use makeDataReport with HTML as output
makeDataReport(comciq, output = "html", replace = TRUE)


# NOT TOO STABLE
# install the DataExplorer library
# install.packages("DataExplorer")
# load the DataExplorer library
library(DataExplorer) 
# use create_report 
create_report(comciq)


# TAKES FOREVER
# install the SmartEDA library
# install.packages("SmartEDA")
# load the SmartEDA library
library(SmartEDA) 
# use ExpReport
ExpReport(comciq_flevel, op_file = 'SmartEDA_df.html')
