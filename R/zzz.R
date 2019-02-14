.onAttach <- function(libname, pkgname) {

    if (!interactive() || stats::runif(1) > 0.1)
        return()
    welcome.message <- paste0(
        " ==================================================================================\n",
        "  _______ ______ _______ _______ __     __         __ __         __          \n",
        " |_     _|      |     __|   _   |  |--.|__|.-----.|  |__|.-----.|  |--.-----.\n",
        "   |   | |   ---|    |  |       |  _  ||  ||  _  ||  |  ||     ||    <|__ --|\n",
        "   |___| |______|_______|___|___|_____||__||_____||__|__||__|__||__|__|_____|\n",
        "     \n",
        "                      _______ _______ _______                        \n" ,
        "                     |     __|   |   |_     _|                \n"   ,
        "                     |    |  |   |   |_|   |_               \n"  ,
        "                     |_______|_______|_______|                 \n"  ,
        "     \n",
        "                           Version:", utils::packageVersion("TCGAbiolinksGUI"), "\n",
        " ==================================================================================\n",
        " Use suppressPackageStartupMessages to eliminate    \n",
        " package startup messages.")
    packageStartupMessage(welcome.message)
}

