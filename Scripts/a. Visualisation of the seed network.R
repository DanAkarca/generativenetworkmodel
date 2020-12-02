# Akarca D, et al. A generative network model of neurodevelopment.
# Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
# University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
# Email: danyal.akarca@mrc-cbu.cam.ac.uk

# 2nd December 2020.

# This script produces cortical visualisations of the seed network

# clear environment
rm(list=ls())

# initialise (note, downloaded devtools followed by a devtools installation)
# https://research.shalailahaas.com/2020/05/16/ggseg-visualization-with-freesurfer-data-in-r-installation/
# e.g. devtools::install_github("LCBC-UiO/ggseg", build_vignettes = FALSE)

library(ggseg)
library(ggseg3d)
library(ggsegExtra)
library(dplyr)
library(ggplot2)

# common to all subjects
valuesData = data.frame(
  region = c("bankssts","caudal anterior cingulate","caudal middle frontal",
  "cuneus","entorhinal","fusiform",
  "inferior parietal","inferior temporal","isthmus cingulate",
  "lateral occipital","lateral orbitofrontal","lingual",
  "medial orbitofrontal","middle temporal","parahippocampal",
  "paracentral","pars opercularis","pars orbitalis",
  "pars triangularis","pericalcarine","postcentral",
  "posterior cingulate","precentral","precuneus",
  "rostral anterior cingulate","rostral middle frontal","superior frontal",
  "superior parietal","superior temporal","supramarginal",
  "frontal pole","temporal pole","transverse temporal",
  "insula","bankssts","caudal anterior cingulate",
  "caudal middle frontal","cuneus","entorhinal",
  "fusiform","inferior parietal", "inferior temporal",
  "isthmus cingulate","lateral occipital","lateral orbitofrontal",
  "lingual","medial orbitofrontal","middle temporal",
  "parahippocampal","paracentral","pars opercularis",
  "pars orbitalis","pars triangularis","pericalcarine",
  "postcentral","posterior cingulate","precentral",
  "precuneus","rostral anterior cingulate","rostral middle frontal",
  "superior frontal","superior parietal","superior temporal",
  "supramarginal","frontal pole","temporal pole",
  "transverse temporal","insula"),
  p = c(1,0,0,0,2,1,2,0,0,1,1,1,2,0,0,0,0,0,0,1,0,1,1,0,1,2,2,2,0,0,0,0,1,0,1,1,
        1,0,0,2,2,2,0,1,0,0,0,4,0,0,0,0,0,0,2,0,2,1,0,1,2,2,2,1,0,0,0,1),
  hemi=c("right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
         "right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
         "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left",
         "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left"),
  stringasFactors = FALSE)

ggseg(.data=valuesData, mapping=aes(fill=p),position="stacked",colour="black",size=.7)+ theme_void()
