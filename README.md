# Code for the paper "Diet composition uncertainty determines impacts on fisheries following an oil spill", Morzaria-Luna et al. 2018

### Data is available through the Gulf of Mexico Research Initiative Information & Data Cooperative (GRIIDC) repository at https://data.gulfresearchinitiative.org/data/R6.x805.000:0002 [DOI:10.7266/N71G0JSG]

We performed an uncertainty analysis for the Atlantis ecosystem model of the Gulf of Mexico (Atlantis-GOM), under a scenario simulating the effects of the Deepwater Horizon oil spill. We used fish stomach content data to inform parameter distribution for the Atlantis-GOM availability matrix, which represents predator total consumption potential and diet preference. We sampled the fish diet composition distribution and analyzed the variability of functional group biomass and catch predicted by Atlantis-GOM simulations to changes in the availability matrix. We used biomass and catch outputs to fit statistical emulators of thr Atlantis-GOM and then predict biomass and catch given the complete diet parameter space. We used simulated and emulated data to assess changes in recovery time to oil spill effects. 

#### Scripts to set up Azure virtual machines for data analysis and to run the Atlantis Ecosystem Model

1. atlantis_template.md. Code you need to set up your Azure virtual VM to run Atlantis, generalize the image and set it up as a template.

2. atlantis_replicate.md. Code to create a replicate from your template and customize it, including setting up Dropbox and compiling Atlantis.
3. atlantis_server.md. Manage and use your server.

