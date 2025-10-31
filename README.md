[![DOI](https://zenodo.org/badge/606954488.svg)](https://zenodo.org/badge/latestdoi/606954488)  
[![FAIR checklist badge](https://fairsoftwarechecklist.net/badge.svg)](https://fairsoftwarechecklist.net/v0.2?f=21&a=32102&i=31101&r=113)

# Analysis of ecologically relevant sea ice and ocean variables for the Southern Ocean using a high-resolution model to inform ecosystem studies
This repository contains all material necessary to reproduce the figures and summary statistics presented in the publication *Analysis of ecologically relevant sea ice and ocean variables for the Southern Ocean using a high-resolution model to inform ecosystem studies* by [Denisse Fierro-Arcos](https://github.com/lidefi87), [Stuart Corney](https://www.utas.edu.au/profiles/staff/imas/stuart-corney), [Amelie Meyer](https://www.utas.edu.au/profiles/staff/imas/amelie-meyer), [Hakase Hayashida](https://github.com/hakaseh), [Andrew E. Kiss](https://github.com/aekiss) and [Petra Heil](https://www.antarctica.gov.au/science/meet-our-scientists/dr-petra-heil/). This manuscript has been published in *Progress in Oceanography*: https://doi.org/10.1016/j.pocean.2023.103049.  
  
In this study we examined the suitability of using outputs from the second run of [ACCESS-OM2-01](https://cosima.org.au/index.php/2020/07/29/data-available-0-1-1958-2018-access-om2-iaf-run/), a high-resolution coupled ocean-sea ice model, to answer questions about ecological impacts in the Southern Ocean. These notebooks can be used as a template for testing the suitability of model outputs for ecological applications, as well as quantitative estimates of changes in key environmental variables for the Southern Ocean over the past 50 years. The [Marine Ecosystem Assessment for the Southern Ocean (MEASO) regions](https://sokiaq.atlassian.net/wiki/spaces/MEASO/pages/4348444548/MEASO+Approach+Structure+Format) were used to evaluate and quantify the rate of change in the physical environment of the Southern Ocean. These regions were designed to establish a standard spatial scale for reporting and assessing environmental and ecosystem change in the SO, and to facilitate comparisons across studies and throughout time.  
  
The workflow presented in these notebooks can be adapted to evaluate different physical variables of ecological relevance and to outputs from different ocean models from an ecological perspective, as well as using different regional boundaries to examine change.  

## How to cite
You can access and use the code contained in this repository as described in the licence. If using this code as a basis for your work, remember you must cite its use using the following citation:  
- Fierro-Arcos, D. (2023). Assessing the suitability of ACCESS-OM2-01 outputs for ecological applications (Version 1.0.0) [Computer software]. [https://doi.org/10.5281/zenodo.7700075](https://doi.org/10.5281/zenodo.7700075)  
  
When using the code in a publication, please also include the following citation in addition to the citation above:  
- Fierro-Arcos, D., Corney, S., Meyer, A., Hayashida, H., Kiss, A. E., & Heil, P. (2023). Analysis of ecologically relevant sea ice and ocean variables for the Southern Ocean using a high-resolution model to inform ecosystem studies. Progress in Oceanography, 215. [https://doi.org/10.1016/j.pocean.2023.103049](https://doi.org/10.1016/j.pocean.2023.103049)


## ACCESS-OM2-01 model outputs
ACCESS-OM2-01 is managed by the [Consortium for Ocean-Sea Ice Modelling in Australia (COSIMA)](https://cosima.org.au/) and they have made their outputs available through the Gadi supercomputer, which is managed by the [National Computational Infrastructure (NCI)](https://nci.org.au/). The COSIMA community has developed the [COSIMA Cookbook](https://github.com/COSIMA/cosima-cookbook/wiki) to allow users to search outputs with ease.  
  
To access ACCESS-OM2-01 outputs you will need to [create an NCI account](https://opus.nci.org.au/display/Help/How+to+create+an+NCI+user+account) with access to projects `ik11`, `cj50`, `jk72`, and `hh5`. These projects either host the data or give you access to the COSIMA cookbook, which we use here to query ACCESS databases.  Non-NCI users can download the output data via the [NCI Data Catalogue](https://dx.doi.org/10.25914/608097cb3433f).

## Observational data
We used two datasets to assess the accuracy of ACCESS-OM2-01 in replicating past environmental conditions:  
- Daily sea ice concentrations from the [NASA Goddard-merged Near Real Time NOAA/NSIDC Climate Data Record of Passive Microwave Sea Ice Concentration](https://climatedataguide.ucar.edu/climate-data/sea-ice-concentration-data-nasa-goddard-and-nsidc-based-nasa-team-algorithm) (version 3)
- Global climatological monthly mixed layer depth means from [Sallee and collaborators](https://doi.org/10.5281/zenodo.4073174)
  
The sea ice concentration data are available in Gadi. The mixed layer depth data is available in [Zenodo](https://zenodo.org/record/5776180) and it was stored in a folder called `Observations` in the root directory of this repository.  
  
## MEASO regions
The MEASO regions boundaries were obtained from the [`measoshape`](https://australianantarcticdivision.github.io/measoshapes/) package for `R`. A copy of these boundaries is included in this repository for anyone who is not familiar with `R`. Only one mask used in these notebooks has been included in this repository under the `SupportingData/Masks` folder due to file size restrictions in GitHub. Three masks are used in these notebooks, one matching the ACCESS-OM2-01 grid, one matching the sea ice concentration data from observations, and another one matching the mixed layer depth observations. All masks can be reproduced using the [`0_CreatingMeasoMask.ipynb`](https://github.com/lidefi87/ACCESS-OM2-01_EcologicallyRelevantVariables/blob/main/Scripts/0_CreatingMeasoMask.ipynb) script.  
  
## Requirements to run these notebooks
Given that the ACCESS-OM2-01 outputs are only available through Gadi, the notebooks will only run in there. Ensure you have access to a project in Gadi that has computational allocation to be able to run them without any issues. All notebooks in this repository have been developed using the [`Analysis3-22.10`](https://github.com/coecms/conda-envs/releases/tag/analysis3-22.10) conda environment available in Gadi.

## Citation
This repository can be cited as: Fierro-Arcos, D. (2023). Assessing the suitability of ACCESS-OM2-01 outputs for ecological applications (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.7700075.

