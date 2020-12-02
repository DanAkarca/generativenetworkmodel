# Akarca D, et al. A generative network model of neurodevelopment.
# Code written by Dr Danyal Akarca, MRC Cognition and Brain Sciences Unit.
# University of Cambridge: https://www.neuroscience.cam.ac.uk/directory/profile.php?da434
# Email: danyal.akarca@mrc-cbu.cam.ac.uk

# 2nd December 2020.

# This script produces visualisations of regional statistical properties of the observed and simulated networks 
# according to the homophily generative mechanisms. It also produces visualisations of the generative errors.

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

# set priors
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
           "transverse temporal","insula")

hemi=c("right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
       "right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right","right",
       "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left",
       "left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left","left")

# measures
observed_degree =  c(4.3148,4.9926,4.7667,2.2407,7.8963,8.3963,8.2037,4.3185,8.7667,7.5296,7.5111,5.8185,8.5815,3.6963,6.6852,6.9370,
                         3.8593,5.8185,5.8741,6.7889,5.6593,10.6556,10.7926,3.6778,8.6741,15.7370,11.6963,9.3074,7.9259,1.6074,0,1.9630,
                         10.5741,3.2889,3.9407,4.2519,4.6593,4.5593,1.9481,7.5148,9.5889,7.7778,4.8593,8.4444,7.9148,8.3556,6.2333,9.7370,
                         3.8148,7.1148,5.8370,4.7407,6.2037,6.3296,6.6407,5.6815,10.3815,11.2074,3.2778,8.5704,15.4000,13.5148,9.6333,8.6741,
                         16,3.5889,1.5704,10.3852) # l. frontal pole is 2.0667 but 16 for caxis, r. frontal pole is 3.8259 but 0 for caxis
simulated_degree =  c(7.3778,6.7778,3.7481,5.0593,9.2444,7.1074,8.6852,4.2741,5.8963,4.7148,8.3519,5.8889,10.9593,4.3333,4.0407,5.2148,
                        5.0667,3.9481,4.8519,7.6444,3.9852,7.6926,6.7444,5.8889,8.8667,7.8926,8.2370,8.4111,5.1259,3.9333,0,4.0222,
                        9.7852,6.5630,11.5444,9.2593,7.4222,5.4148,4.6852,11.4000,11.1778,10.6185,6.1037,8.7926,4.6630,5.1259,4.6741,12.1593,
                        5.4667,4.4630,4.8519,4.1630,4.6148,5.7185,8.1852,6.0074,8.2333,9.4074,7.0037,6.8481,7.6815,8.5926,11.5111,8.0111,
                        16,3.7222,7.3222,11.3333) # l. frontal pole is 3.100 but 16 for caxis, r. frontal pole is 3.1889 but 0 for caxis

observed_clustering =  c(0.6932,0.7161,0.4929,0.5627,0.3823,0.5461,0.5109,0.6496,0.4432,0.2994,0.3525,0.2475,0.5334,0.5278,0.4490,0.5569,
                         0.6455,0.5468,0.4610,0.5284,0.5508,0.3634,0.2855,0.6321,0.3755,0.2173,0.2985,0.4205,0.5655,0.2944,0,0.7156,
                         0.3461,0.3080,0.9248,0.6911,0.7326,0.5359,0.3732,0.4278,0.4772,0.5345,0.5518,0.5057,0.2987,0.3560,0.2379,0.4892,
                         0.5226,0.4320,0.6053,0.5594,0.4836,0.4651,0.5453,0.5527,0.3501,0.2805,0.6710,0.3512,0.2143,0.2814,0.4116,0.5281,
                         1,0.3369,0.4488,0.3603) # l. frontal pole is 0.3759 but 1 for caxis, r. frontal pole is 0.3827 but 0 for caxis
simulated_clustering =  c(0.4327,0.3689,0.3038,0.3753,0.4085,0.4208,0.4262,0.3634,0.3767,0.4519,0.4461,0.4069,0.3871,0.3623,0.3309,0.3499,
                          0.3612,0.3621,0.3195,0.4024,0.3146,0.3848,0.4542,0.3774,0.4430,0.4506,0.3987,0.4135,0.3349,0.2991,0,0.3277,
                          0.3944,0.3819,0.4921,0.4227,0.4900,0.3563,0.3415,0.4715,0.4492,0.5211,0.3512,0.5020,0.3639,0.3571,0.3315,0.4833,
                          0.3141,0.3011,0.3927,0.3265,0.3344,0.3942,0.4684,0.3490,0.4730,0.4321,0.3871,0.5077,0.4789,0.4827,0.4999,0.4996,
                          1,0.3265,0.3888,0.4541) # l. frontal pole is 0.3206 but 1 for caxis, r. frontal pole is 0.3639 but 0 for caxis

observed_betweenness =  c(17.7724,19.0133,54.1358,7.1135,129.6974,80.7821,64.2690,32.2786,134.3674,155.0926,127.6599,158.2244,56.3856,49.6272,66.1882,36.3058,
                          13.2972,35.1626,102.3043,76.5717,35.8217,253.1742,368.6077,13.1741,188.2963,628.1426,467.1063,127.9748,50.9835,6.8054,0,0.0898,
                          184.8693,47.5725,0.9998,16.1920,14.6107,43.0366,7.7418,110.9196,114.8103,61.2970,48.1564,87.3251,172.4425,130.7194,199.2516,84.9265,
                          42.7710,77.7460,25.2505,25.3929,53.8157,90.1260,72.2341,36.9907,265.9208,366.7285,13.2522,205.5081,595.7353,593.3772,164.7598,81.4233,
                          700,40.0867,0.0353,153.1629) # l. frontal pole is 13.3127 but 700 for caxis, r. frontal pole is 33.5628 but 0 for caxis
simulated_betweenness =  c(118.9959,95.4070,46.2830,65.0201,133.2062,109.4777,120.0885,54.4173,81.5552,51.7499,104.9943,76.5535,188.2678,54.4083,54.7818,75.7166,
                           78.9352,48.5457,65.4983,124.8135,49.5667,132.8454,98.6923,75.6499,122.9893,114.2238,146.2879,118.6062,70.7161,54.5116,0,55.1289,
                           155.8386,100.8819,109.3015,155.2325,82.3609,66.4443,56.7300,127.0162,118.9124,92.0023,88.8993,87.3980,69.9615,66.6665,63.7596,121.0312,
                           70.7592,66.4034,58.1338,56.4049,59.2146,71.4007,107.8997,98.9050,102.7690,119.4602,101.1423,77.7579,96.5816,87.8693,112.3620,88.9282,
                           700,46.6299,79.7198,145.8466) # l. frontal pole is 37.2817 but 700 for caxis,  r. frontal pole is 32.5817 but 0 for caxis


observed_length =  c(267.0180,240.1521,281.3447,209.9105,444.2462,499.8768,470.1245,294.4289,407.7684,640.4898,235.0030,415.9026,360.3649,45.4613,314.3159,316.1037,
                     100.0272,160.4112,195.3407,400.0476,128.8929,281.8192,348.2762,129.7990,128.1515,324.3493,421.3419,307.9021,21.0946,75.1370,0,24.4596,
                     3.0491,242.1277,122.9205,162.5748,174.9825,113.8723,55.9308,303.4103,430.0837,294.7733,153.1258,327.2733,293.8612,270.3375,124.9431,217.6821,
                     51.2637,137.1467,173.2305,126.4752,182.8735,86.1857,145.6277,99.3076,165.3896,118.7653,45.3165,82.2582,38.7300,138.4200,87.6117,58.5258,
                     700,6.4931,7.3475,0) # l. frontal pole is 4.5111 but 700 for caxis, r. frontal pole is 40.7145 but 0 for caxis

simulated_length =  c(477.1025,332.0726,221.0429,264.3290,574.6825,431.4743,528.7534,211.9193,270.7494,226.0151,344.3029,229.1427,498.6963,181.9994,188.8531,214.4759,
                      206.5839,157.4263,190.4743,322.2519,164.5161,254.0687,263.9416,201.7401,257.1897,207.6875,209.8235,304.3149,125.9610,119.7818,0,117.1294,
                      177.4102,143.6732,436.2528,245.5328,263.7342,130.8557,144.8667,372.0578,392.1043,290.0662,120.8926,228.8836,103.1050,87.0276,81.3780,245.3663,
                      95.2782,75.5214,88.4364,76.6238,77.3181,65.8270,126.2286,70.4692,108.3326,96.8573,52.1002,85.1197,44.0319,69.2001,50.1995,25.2421,
                      700,10.3571,11.1349,0) # l. frontal pole is 9.7725 but 700 for caxis, r. frontal pole is 111.8208 but 0 for caxis

observed_efficiency = c(0.8391,0.8488,0.6302,0.6046,0.5907,0.7601,0.7271,0.7791,0.6655,0.4503,0.5448,0.3343,0.7485,0.6470,0.6589,0.7674,
                        0.7532,0.7241,0.6178,0.7275,0.7330,0.6147,0.5012,0.7485,0.6036,0.4772,0.5327,0.6521,0.7749,0.3189,0,0.7244,
                        0.6062,0.3577,0.9617,0.8380,0.8532,0.6668,0.4016,0.6213,0.7147,0.7343,0.6975,0.7190,0.4442,0.5636,0.3256,0.7254,
                        0.6391,0.6456,0.7936,0.7084,0.6680,0.6253,0.7269,0.7341,0.5999,0.5113,0.7860,0.5576,0.4703,0.5185,0.6277,0.7530,
                        1,0.3955,0.4503,0.6214) # l. frontal pole is 0.4100 but 1 for caxis, r. frontal pole is ,0.4664 but 0 for caxis

simulated_efficiency = c(0.6056,0.5220,0.3940,0.4912,0.6345,0.5877,0.6472,0.4621,0.5159,0.5619,0.6465,0.5446,0.6285,0.4686,0.4173,0.4717,
                         0.4815,0.4485,0.4284,0.5730,0.4145,0.5572,0.6182,0.5155,0.6427,0.6316,0.5838,0.6273,0.4563,0.3915,0,0.4213,
                         0.6217,0.5351, 0.7171,0.6273,0.6764,0.4812,0.4473,0.6981,0.6808,0.7287,0.4884,0.6919,0.4684,0.4674,0.4358,0.7132,
                         0.4289,0.4009,0.4990,0.4208,0.4377,0.5143,0.6658,0.4802,0.6740,0.6401,0.5440,0.6770,0.6636,0.6848,0.7247,0.6858,
                         1,0.4033,0.5382,0.6816) # l. frontal pole is 0.3889 but 1 for caxis, r. frontal pole is 0.4356 but 0 for caxis

# degree

observed_degreeData = data.frame(region,observed_degree,hemi,stringasFactors = FALSE)
ggseg(.data=observed_degreeData, mapping=aes(fill=observed_degree),position="stacked",colour="black",size=0.7)+ theme_void()

simulated_degreeData = data.frame(region,simulated_degree,hemi,stringasFactors = FALSE)
ggseg(.data=simulated_degreeData, mapping=aes(fill=simulated_degree),position="stacked",colour="black",size=0.7)+ theme_void()

# clustering
observed_clusteringData = data.frame(region,observed_clustering,hemi,stringasFactors = FALSE)
ggseg(.data=observed_clusteringData, mapping=aes(fill=observed_clustering),position="stacked",colour="black",size=0.7)+ theme_void() +
scale_fill_gradient(low="dark orange",high="dark red")

simulated_clusteringData = data.frame(region,simulated_clustering,hemi,stringasFactors = FALSE)
ggseg(.data=simulated_clusteringData, mapping=aes(fill=simulated_clustering),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradient(low="dark orange",high="dark red")

# betweenness
observed_betweennessData = data.frame(region,observed_betweenness,hemi,stringasFactors = FALSE)
ggseg(.data=observed_betweennessData, mapping=aes(fill=observed_betweenness),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradient(low="grey",high="gold")

simulated_betweennessData = data.frame(region,simulated_betweenness,hemi,stringasFactors = FALSE)
ggseg(.data=simulated_betweennessData, mapping=aes(fill=simulated_betweenness),position="stacked",colour="black",size=0.7) + theme_void() +
  scale_fill_gradient(low="grey",high="gold")

# edge length
observed_lengthData = data.frame(region,observed_length,hemi,stringasFactors = FALSE)
ggseg(.data=observed_lengthData, mapping=aes(fill=observed_length),position="stacked",colour="black",size=0.7) + theme_void() +
  scale_fill_gradient(low="white",high="purple")

simulated_lengthData = data.frame(region,simulated_length,hemi,stringasFactors = FALSE)
ggseg(.data=simulated_lengthData, mapping=aes(fill=simulated_length),position="stacked",colour="black",size=0.7) + theme_void() +
  scale_fill_gradient(low="white",high="purple")

# local efficiency
observed_efficiencyData = data.frame(region,observed_efficiency,hemi,stringasFactors = FALSE)
ggseg(.data=observed_efficiencyData, mapping=aes(fill=observed_efficiency),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradient(low="white",high="firebrick")

simulated_efficiencyData = data.frame(region,simulated_efficiency,hemi,stringasFactors = FALSE)
ggseg(.data=simulated_efficiencyData, mapping=aes(fill=simulated_efficiency),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradient(low="white",high="firebrick")

# measure error

# degree error
degree_error = c(3.0630,1.7852,-1.0185,2.8185,1.3481,-1.2889,0.4815,-0.0444,-2.8704,-2.8148,0.8407,0.0704,2.3778,0.6370,-2.6444,-1.7222,
                 1.2074,-1.8704,-1.0222,0.8556,-1.6741,-2.9630,-4.0481,2.2111,0.1926,-7.8444,-3.4593,-0.8963,-2.8000,2.3259,-8,2.0593,
                 -0.7889,3.2741,7.6037,5.0074,2.7630,0.8556,2.7370,3.8852,1.5889,2.8407,1.2444,0.3481,-3.2519,-3.2296,-1.5593,2.4222,
                 1.6519,-2.6519,-0.9852,-0.5778,-1.5889,-0.6111,1.5444,0.3259,-2.1481,-1.8000,3.7259,-1.7222,-7.7185,-4.9222,1.8778,-0.6630,
                 8,0.1333,5.7519,0.9481) # l. frontal pole is 1.0333 but 8 for caxis, r. frontal pole is -0.6370 but -8 for caxis

error_degreeData = data.frame(region,degree_error,hemi,stringasFactors = FALSE)
ggseg(.data=error_degreeData, mapping=aes(fill=degree_error),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradientn(colours = c("light blue","white","dark blue"),na.value="grey")

# clustering error
clustering_error = c(-0.2605,-0.3472,-0.1892,-0.1874,0.0261,-0.1253,-0.0847,-0.2862,-0.0665,0.1525,0.0936,0.1594,-0.1463,-0.1655,-0.1181,-0.2069,
                     -0.2843,-0.1847,-0.1415,-0.1260,-0.2363,0.0213,0.1687,-0.2547,0.0675,0.2333,0.1003,-0.0069,-0.2306,0.0047,-0.5,-0.3878,
                     0.0483,0.0739,-0.4327,-0.2684,-0.2426,-0.1797,-0.0318,0.0437,-0.0280,-0.0134,-0.2006,-0.0037,0.0652,0.0011,0.0936,-0.0059,
                     -0.2085,-0.1309,-0.2126,-0.2329,-0.1492,-0.0708,-0.0768,-0.2037,0.1229,0.1516,-0.2839,0.1564,0.2646,0.2012,0.0883,-0.0284,
                     0.5,-0.0104,-0.0600,0.0938) # l. frontal pole is -0.0554 but 0.5 for caxis, r. frontal pole is -0.0188 but -0.5 for caxis

error_clusteringData = data.frame(region,clustering_error,hemi,stringasFactors = FALSE)
ggseg(.data=error_clusteringData, mapping=aes(fill=clustering_error),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradientn(colours = c("dark orange","white","dark red"),na.value="grey")

# betweenness error
betweenness_error = c(101.2235,76.3937,-7.8528,57.9066,3.5088,28.6956,55.8196,22.1386,-52.8122,-103.3428,-22.6656,-81.6708,131.8823,4.7811,-11.4064,39.4109,
                 65.6380,13.3831,-36.8060,48.2418,13.7450,-120.3288,-269.9154,62.4758,-65.3069,-513.9188,-320.8184,-9.3687,19.7326,47.7062,-520,55.0390,
                 -29.0308,53.3093,108.3016,139.0405,67.7502,23.4077,48.9882,16.0966,4.1021,30.7054,40.7429,0.0729,-102.4810,-64.0529,-135.4920,36.1047,
                 27.9883,-11.3427,32.8834,31.0120,5.3989,-18.7253,35.6656,61.9143,-163.1517,-247.2683,87.8902,-127.7502,-499.1537,-505.5079,-52.3978,7.5049,
                 520,6.5432,79.6845,-7.3163) # l. frontal pole is 23.9690 but 520 for caxis, r. frontal pole is -0.9811 but -520 for caxis

error_betweennessData = data.frame(region,betweenness_error,hemi,stringasFactors = FALSE)
ggseg(.data=error_betweennessData, mapping=aes(fill=betweenness_error),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradientn(colours = c("grey","white","gold"),na.value="grey")

# edge length error
length_error = c(210.0846,91.9205,-60.3018,54.4184,130.4363,-68.4025,58.6289,-82.5096,-137.0190,-414.4747,109.2999,-186.7599,138.3314,136.5381,-125.4628,-101.6278,
                 106.5567,-2.9849,-4.8664,-77.7956,35.6232,-27.7505,-84.3345,71.9411,129.0383,-116.6619,-211.5184,-3.5872,104.8664,44.6447,-420,92.6699,
                 174.3610,-98.4545,313.3323,82.9580,88.7517,16.9834,88.9358,68.6476,-37.9794,-4.7071,-32.2332,-98.3897,-190.7562,-183.3100,-43.5652,27.6841,
                 44.0146,-61.6254,-84.7941,-49.8514,-105.5554,-20.3587,-19.3991,-28.8384,-57.0570,-21.9080,6.7837,2.8616,5.3019,-69.2199,-37.4122,-33.2836,
                 420,3.8640,3.7874,0) # l. frontal pole is 5.2614 but 420 for caxis, r. frontal pole is 71.1064 but -420 for caxis

error_lengthData = data.frame(region,length_error,hemi,stringasFactors = FALSE)
ggseg(.data=error_lengthData, mapping=aes(fill=length_error),position="stacked",colour="black",size=0.7)+ theme_void() +
  scale_fill_gradientn(colours = c("grey","white","purple"),na.value="grey")

# average ranked error
ranked_error = c(61,58,26,48,16,21,13,40,43,64,19,45,55,33,42,47,53,24,15,22,38,27,60,51,18,67,63,
                 2,52,14,5,59,23,41,68,62,56,17,28,32,9,11,29,4,54,44,30,10,37,34,39,31,35,7,
                 12,25,46,50,57,36,66,65,20,3,6,1,49,8)

error_rankData = data.frame(region,ranked_error,hemi,stringasFactors = FALSE)
ggseg(.data=error_rankData, mapping=aes(fill=ranked_error),colour="black",size=0.2)+ theme_void() +
  scale_fill_gradient(low="light green",high="firebrick")

