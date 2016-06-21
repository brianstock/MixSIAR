# Brian Stock
# March 10, 2016
# Run all examples

library(MixSIAR)

# find MixSIAR directory
mixsiar.dir <- find.package("MixSIAR")
# folder with script files
paste0(mixsiar.dir,"/example_scripts")

# run Wolves example
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_wolves.R"))
graphics.off()

# run Geese example
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_geese.R"))

# run Lake example
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_lake.R"))

# run Palmyra example
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_palmyra.R"))

# run Killer Whale example
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_killerwhale.R"))

# run Snail example
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_snail.R"))

# run Storm-petrel example
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_stormpetrel.R"))

# run Cladocera example
graphics.off()
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_cladocera.R"))

# run Isopod example
graphics.off()
source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_isopod.R"))

