
##  In this repository are:

  - the CCBH-2022 in csv format (CCBH-2022.csv),
  - the python code to generate the input file for the ASSA-PBN software (GRN_to_ASSA.py),
  - the python code to compute conserved genes in last n time steps (find_basin.py),
  - the text containing detailed instructions and command lines (commands.txt),
  - the binarization code of the RNA-seq data in R language (bin.R),
  - core-genes.csv: the core subnetwork containing 212 genes divided into 166 with both incoming and outgoing edges, and 46 genes with only outgoing edges.
  - Simulation trajectory results and attraction basins (sim.zip). The normalization and binarization results are also in this file. For enhanced folder navigation, it is important to read the README file contained within it.

Algorithms designed by Marcelo Trindade dos Santos (LNCC/MCTI), member of the Computational Modeling of Multidrug-resistant Bacteria research group (Fiocruz).


##  **Information about the RNA-seq _bulk_ experiments:**
 
For Antibiotic Growth Curves:
  - _Pseudomonas aeruginosa_ was incubated for 10 hours (growth curve) in Mueller-Hinton medium. After this period, antimicrobials were added:
  - Imipenem = 128 mg/l
  - Polymyxin B = 1 mg/ml
  - Control (without addition of antibiotic)
    
  - The cultures were further incubated for 1 hour, and then 1 ml was collected from each tube to perform RNA extraction (500 μl from each tube + 2 volumes of RNA protect).

  - For biofilm studies, the same inoculum used for the growth curves was used to inoculate a 12-well plate. After 24 hours of incubation at 37°C, the planktonic cells were transferred to tubes for OD600 reading.
  - The biofilm adhered to the plate was suspended in 3 ml of PBS (same volume as used at the beginning of the experiment). OD600 reading was also performed on this suspension.

  - 1 ml of the biofilm suspension was collected for RNA extraction. RNAs were extracted using the RNeasy Minikit from Qiagen.

  - For experiments with different carbon sources, M9 medium was used, and sugars were added at 40 mM concentration. Pseudomonas was inoculated, and a growth curve was generated. As they reached mid-log phase, aliquots were withdrawn for RNA extraction.

  - The planktonic samples were grown in tubes (16x100mm) containing 3 ml of Mueller Hinton II medium (same medium used for biofilm and other tests). The tubes were incubated at 37°C with agitation (150 rpm). After 11 hours of incubation, they reached a OD600 of C9 = 0.75 and C10 = 0.78. 500 μl were collected from each tube for RNA extraction.

  - For the biofilms, a 12-well plate was used. 3 ml of medium was added to each well. The plate was incubated for 24 hours at 37°C without agitation to allow biofilm formation. After incubation, suspended cells were discarded, and the wells were washed with PBS to remove non-adherent cells. Then, the biofilm was disrupted (3 ml of PBS was added, and cells were loosened using a pipette). The volume from each well was transferred to a sterile glass tube for OD600 reading. C11 = 1.27 and C12 = 1.14. 500 μl were collected from each well for RNA extraction. The biofilm is three-dimensional, but this is only observable under microscopy.

  - To the collected aliquots for RNA extraction, RNAprotect (to protect RNA from degradation) was added, and extraction was performed using the Qiagen Rneasy Minikit.
