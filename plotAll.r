library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(tibble)
library(ggpmisc)

source("parseInputData.r")

source("plotMycBaspRatios.r")
source("plotMycPhbRatios.r")
source("plotMycAbsolute.r")
source("plotBaspAbsolute.r")
#source("RPKMvsRSEMdemo.r") does not work currently w/ new APIrequest input system