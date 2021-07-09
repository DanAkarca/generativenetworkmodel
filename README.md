**Code repository for “A generative network model of neurodevelopmental diversity in structural brain organization”**
doi:10.5281/zenodo.4762612

Danyal Akarca [1], Petra E Vértes [2,3], Edward T Bullmore [2,4], the CALM team [1] & Duncan E Astle [1].
1. MRC Cognition and Brain Sciences Unit, University of Cambridge, Cambridge, UK
2. Department of Psychiatry, University of Cambridge, Cambridge, UK
3. The Alan Turing Institute, London, UK
4. Department of Clinical Neurosciences, Wolfson Brain Imaging Centre, University of Cambridge, Cambridge, UK

For any questions regarding the use of this repository, please get in touch at Danyal.akarca@mrc-cbu.cam.ac.uk.

If using the code, please consider citing our paper:

**Akarca, D., Vértes, P.E., Bullmore, E.T. et al. A generative network model of neurodevelopmental diversity in structural brain organization. Nat Commun 12, 4216 (2021). https://doi.org/10.1038/s41467-021-24430-z**

**Requirements**

The following installations are required to use all the attached scripts. However, most will be usable with MATLAB alone. Installation time on a typical computer should take no longer than 60 minutes.
* MATLAB 2019b (installation: https://uk.mathworks.com/help/install/install-products.html) 
* RStudio 1.2.5033 (installation: https://rstudio.com/products/rstudio/download/) 
* Python 3.7.0 (installation: https://www.python.org/downloads/) 
* Brain Connectivity Toolbox, 2019 (installation: https://sites.google.com/site/bctnet/)

**Data availability statement**

The datasets supporting the current study have not been deposited in a public repository because of restrictions imposed by NHS ethical approval, but may be available upon request. As such, all data in this repository and unidentifiable. They are either produced from simulations alone or are random data.

**/Example data**

Relevant data can be found at https://osf.io/h9px4/?view_only=984260dcff444b59819961ece9c724ec

This contains relevant, publicly-available, data that is relevant for the present study and can be used with code within the repository. They contain no identifiable data. They are not necessarily representative of the data used for the current study, but allow for the pipeline to be run.

**/Scripts**

The following MATLAB code was used in our analyses. Each script pertains to a separate part of the workflow. They are annotated within the script. The run time for each script can vary considerably. Most can be run within 3 minutes. Those that loop over multiple subjects can take much longer, depending on hyperparameters. In most cases, the estimated time is provided within the script. They include:

i. Connectome thresholding.m

ii. Computing the seed network.m

iii. Running initial generative models.m

iv. Running the homophily generative model.m

v. Exploring generative model outputs.m

vi. Embedding, errors and predictions.m

vii. Gene expression PLS analysis.m

**The following R and Python code was used to produce some of the visuals. They include:**

a. Visualisation of the seed network.R

b. Visualisation of nodal values and costs.R

c. Visualisation of spatial embedding and errors.R

d. Visualisation of wiring parameters in a radar plot.ipynb

**/Scripts/connectome_pipeline**

This contains all the relevant code that was used for our diffusion tensor imaging analysis, using NiPype. Code has been adapted for a new MRI scanner at our locality, from previously published work. A prior version is found at https://github.com/JoeBathelt/Modularity/. For details regarding the specifics of this pipeline, please get in touch directly.
