Dataset Attribution and Usage
-----------------------------
 
* Dataset title: Linking microenvironment modification to species interactions and demography in an alpine plant community

* Dataset description: Demography and spatial dataset from an alpine plant community in southwestern Colorado. Includes survival, growth, and reproduction data from eight census years.

* Persistent Identifier: https://doi.org/10.6078/D1342B

* Dataset contributors: Courtenay A. Ray, Rozália E. Kapás, Øystein H. Opedal, Benjamin W. Blonder

* Journal: Oikos

* Collection methods: Please see associated publication at https://doi.org/10.1111/oik.09235

* License: This work is licensed under a CC0 1.0 Universal (CC0 1.0) Public Domain Dedication license. We are very interested in collaborating with these data. Please, let us know if you are too!

* Suggested dataset citation: Ray, Courtenay A.; Kapás, Rozália E.; Opedal, Øystein H.; Blonder, Benjamin W. (2022), Data from: Linking microenvironment modification to species interactions and demography in an alpine plant community, Dryad, Dataset, https://doi.org/10.6078/D1342B

* Funding: 
  * Rocky Mountain Biological Laboratory Graduate Fellowship, Snyder Endowment
  * Rocky Mountain Biological Laboratory Graduate Fellowship, Langenheim Endowment
  * Colorado Mountain Club Foundation, Award: Kurt Gerstle Fellowship
  * Arizona State University, Award: Research and Training Initiatives Graduate Student Support
  * Peder Sather Grant
  * Norges Forskningsråd KLIMAFORSK, Award: 250233
  * Natural Environment Research Council, Award: NE/M019160/1
  * Norwegian Research Council Center of Excellence, Award: 223257
  * KLIMAFORSK
  
---

Description of Files

Details for: Density_Overlap_Demog_Data_2014_2019.csv
---------------------------------------

*Description: Contains demographic and precipitation data from 2014 to 2019 as well as estimates of aboveground spatial overlaps

* Size: 907 KB

* Variables:
  * Year: Year of data collection (2014-2019)	
  * Plot: 2 x 2 m census plot (1-50)	
  * Tag: Individual numeric identifier	
  * X.global..cm.: X position (cm) of plant within entire field site with respect to the uppermost left corner (e.g. Top left of plot 1)
  * Y.global..cm.: Y position (cm) of plant within entire field site with respect to the uppermost left corner (e.g. Top left of plot 1)	
  * Taxon: Genus and species of focal plant
  * Length..cm.: Length of plant across longest axis (cm)
  * Length..cm..prev: Length of plant across longest axis (cm) in prior year
  * Growth: Difference in length (cm) between current and prior year	
  * Height..cm.: Maximum height of plant (cm)
  * X..Capitulescences: Number of capitulescences on focal plant	
  * Is.seedling: Whether plant is a new recruit from seed (0/1)	
  * Is.recruit: Whether the plant is a new recruit from seed or vegetative reproduction (0/1)
  * Survival: Plants have a survival==0 in the first year of two consecutive years of size==0 or the last year present in the data before 'ghosting'	
  * ghost.dead: Ghost dead plants (==1) are plants that are alive (sizes>0) that have survival scores of 0 because they 'ghost'. Ghosting is when a plant leaves a data set without two years of being recorded dead in the census. The biggest source of ghosting is due to merging and dying individuals in multi-individual patches. The only truly dead plants in the demography data set have size==0 and survival==0. 
  * true.dead: 1==plant is truly dead, i.e., plants that are (size==0) and (survival==0)
  * avg.precip.mm: Average daily precipitation during the growing season (mm)
  * Inter.Overlap.Area:	Percentage of the focal individual that is overlapped by interspecific taxa based on aboveground spatial extents 
  * Intra.Overlap.Area: Percentage of the focal individual that is overlapped by intraspecific taxa based on aboveground spatial extents 


Details for: Belowground_Density_Overlap_Demog_Data_2014_2019.csv
---------------------------------------

*Description: Contains demographic and precipitation data from 2014 to 2019 as well as estimates of belowground spatial overlaps for each individual

* Size: 930 KB

* Variables:
  * Year: Year of data collection (2014-2019)	
  * Plot: 2 x 2 m census plot (1-50)
  * Tag: Individual numeric identifier	
  * X.global..cm.: X position (cm) of plant within entire field site with respect to the uppermost left corner (e.g. Top left of plot 1)
  * Y.global..cm.: Y position (cm) of plant within entire field site with respect to the uppermost left corner (e.g. Top left of plot 1)	
  * Taxon: Genus and species of focal plant
  * Length..cm.: Length of plant across longest axis (cm)
  * Length..cm..prev: Length of plant across longest axis (cm) in prior year
  * Growth: Difference in length (cm) between current and prior year	
  * Height..cm.: Maximum height of plant (cm)
  * X..Capitulescences: Number of capitulescences on focal plant	
  * Is.seedling: Whether plant is a new recruit from seed (0/1)	
  * Is.recruit: Whether the plant is a new recruit from seed or vegetative reproduction (0/1)
  * Survival: Plants have a survival==0 in the first year of two consecutive years of size==0 or the last year present in the data before 'ghosting'	
  * ghost.dead: Ghost dead plants (==1) are plants that are alive (sizes>0) that have survival scores of 0 because they 'ghost'. Ghosting is when a plant leaves a data set without two years of being recorded dead in the census. The biggest source of ghosting is due to merging and dying individuals in multi-individual patches. The only truly dead plants in the demography data set have size==0 and survival==0. 
  * true.dead: 1=plant is truly dead, i.e., plants that are (size==0) and (survival==0)
  * avg.precip.mm: Average daily precipitation during the growing season (mm)
  * Inter.Overlap.Area:	Percentage of the focal individual that is overlapped by interspecific taxa based on belowground spatial extents 
  * Intra.Overlap.Area: Percentage of the focal individual that is overlapped by intraspecific taxa based on belowground spatial extents 

Details for: demography_2014-2021_manual.edits.csv
---------------------------------------

*Description: Contains the full demography data from the census plot from 2014 to 2021.

* Size: 774 KB 

* Variables:
  * Year: Year of data collection (2014-2021)
  * Plot: 2 x 2 m census plot (1-50)
  * Tag: Individual numeric identifier	
  * Taxon: Genus and species of focal plant	
  * X..cm.: X position (cm) of plant within plot with respect to the lower left corner
  * Y..cm.: Y position (cm) of plant within plot with respect to the lower left corner	
  * Tag.offset.x: X position (cm) of tag with respect to the focal plant
  * Tag.offset.y: Y position (cm) of tag with respect to the focal plant	
  * Length..cm.: Length of plant across longest axis (cm)
  * Height..cm.: Maximum height of plant (cm)
  * X..Capitulescences: Number of capitulescences on focal plant
  * Is.seedling: Whether plant is a new recruit from seed (0/1)
  * Died.this.census.final: Whether the plant died in current year (0/1)
  * Length..cm..prev: Length of plant across longest axis (cm) in prior year
  * Growth: Difference in length (cm) between current and prior year
  * Survival: Plants have a survival==0 in the first year of two consecutive years of size==0 or the last year present in the data before 'ghosting'
  * Is.recruit: Whether the plant is a new recruit from seed or vegetative reproduction (0/1)

Was data derived from another source?
If yes, list source(s): No

- - -
END OF README
