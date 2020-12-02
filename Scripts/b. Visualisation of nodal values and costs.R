# Akarca D, et al. A generative network model of neurodevelopment.
# Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
# University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
# Email: danyal.akarca@mrc-cbu.cam.ac.uk

# 2nd December 2020.

# This script produces cortical visualisations of nodal values and nodal costs, averaged across the sample, generated from the homophily generative models with ggseg.

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

# parameterised values
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
  K = c(9.00319631469687,7.80952205102756,5.62009560895792,6.10141023044186,11.0314434835978,8.88779947003576,11.6050105914385,4.85652847776261,7.14783684153111,6.35589239257115,
        10.1814404818225,7.15538020640347,12.3147138068044,5.55031108650868,5.32073851379553,6.3950386885244,5.81711155689916,5.43309512710916,5.47908153554285,
        8.62292573897021,5.52558422540003,8.7968269163261,8.64679563324738,7.01164227132134,10.4202180801606,9.55644171836597,9.34492194203523,11.3846314827612,6.19869291024068,
        5.05173849338221,4.30280649223875,5.53609371593252,10.666788875843,7.3075702162384,14.2733088682186,10.2691772520641,9.84394604557461,6.61431813150401,5.53213371189131,13.2924283021267,
        13.8740129116333,13.9814740807381,7.15984783673938,11.0950805759163,5.6715624638565,6.22483759754839,6.33152394070051,14.8729461392826,7.41303582353217,
        5.4085895194824,6.57843079729369,5.17891519017658,5.511537772391,6.34663463282856,10.0090341919153,7.09664277140218,10.1937866543313,11.4352797424721,
        7.45787793450807,9.12898968715097,9.59013488487613,11.7432773180615,14.0110439288289,9.87146545235501,4.63977747399204,
        4.66621313020011,8.40038886281483,12.9032617307312),
  hemi=c("right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
         "right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
         "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left",
         "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left"),
  stringasFactors = FALSE)

ggseg(.data=valuesData, mapping=aes(fill=K),position="stacked",colour="black",size=.7)+ theme_void() +
  scale_fill_gradient(low="goldenrod",high="firebrick")

# parameterised costs
costData = data.frame(
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
  D = c(0.00337114474162581,0.00411878482951205,0.0031701175547964,0.00355764061300191,0.00322185109345255,
        0.00344202459248644,0.00282856976279522,0.00337584954272647,0.0039444912679735,0.00273669294189775,
        0.00355383595756841,0.00346329900648464,0.00376206454043992,0.0034778112790836,0.00335770171420816,
        0.00364061778068188,0.00350216466284319,0.00319313055160689,0.00341398401551159,0.00393534097001588,
        0.00330148864527408,0.00394643819800591,0.00338144059666672,0.00373155502330447,0.00424555019493761,
        0.00309069706869843,0.0033882233037397,0.00281770636811846,0.00382414823841649,0.00305651707333514,
        0.0027059830241633,0.00305058589068429,0.00410854176072516,0.0041158272503468,0.00313293775592128,
        0.00417466232611281,0.0030344929770174,0.00397025234225018,0.00323816104308112,0.00340793762587127,
        0.00282854612392152,0.00317717446313397,0.00390634283916468,0.0027606113682659,0.00341057404005369,
        0.00342621647995605,0.003559567449406,0.00332541378215993,0.00359512569597027,0.00357346350732382,
        0.0033561519590417,0.00319057969861316,0.00332231986980347,0.00397012514586261,0.00306139532404386,
        0.00395999797199927,0.00326654452013764,0.00361103870931453,0.0042176034586046,0.0030093041543513,
        0.00312259705951743,0.00285424175964796,0.00352421238844939,0.00289531417472875,0.00287033054061638,
        0.00293413038472717,0.0037422841646931,0.00373582310575959),
  hemi=c("right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
         "right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
         "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left",
         "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left"),
  stringasFactors = FALSE)

ggseg(.data=costData, mapping=aes(fill=D),position="stacked",colour="black",size=.7)+ theme_void() +
  scale_fill_gradient(low="dark blue",high="light blue")
