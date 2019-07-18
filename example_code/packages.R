# ---- First, you need to intall several packages to run the code presented in the lectures ---- # 

# ---- packages ---- #
# Uncomment to proceed...
# midasr #
#install.packages(c("sandwich", "optimx", "quantreg"))
#install.package("midasr")

# GARCH-MIDAS #
#install.packages("mfGARCH")


# Alternativelt, install everything at once, plus adding additional functions: 
#install.packages("devtools")
require("devtools")
install_github("jstriaukas/MIDASLec")
require("MIDASLec")

