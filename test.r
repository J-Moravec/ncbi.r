set_script_dir = function(){
    dir = commandArgs(FALSE) |>
        grep(pattern = "^--file=", value = TRUE) |>
        sub(pattern = "^--file=", replacement = "") |>
        dirname()
    setwd(dir)
    }

# Run from the script's directory when run using Rscript
if(sys.nframe() == 0) set_script_dir()


source("tests/mutr.r")
source("tests/test-mutr.r")
source("tests/helpers.r")
source("ncbi.r")

TEST_INIT()
TEST_FILE("tests/test-ncbi.r")
TEST_PRINT()
