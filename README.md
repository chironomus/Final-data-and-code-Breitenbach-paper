# Final-data-and-code-Breitenbach-paper
This is code and data for the paper " V Baranov, J Jourdan, F Pilotto, R Wagner, P Haase, 2020. 
Complex and nonlinear climate‐driven changes in freshwater insect communities over 42 years
Conservation Biology, 2020 https://doi.org/10.1111/cobi.13477
ABSTRACT
The ongoing biodiversity crisis becomes evident in the widely observed decline in abundance and diversity of species, profound changes in 
community structure and shifts in the species' phenology. Insects are among the most affected groups, with documented decreases 
in abundance up to 76% in the last 25–30 years in some terrestrial ecosystems. Identifying the underlying drivers is a major 
obstacle as most ecosystems are affected by multiple stressors simultaneously and in situ measurements of environmental 
variables are often missing. In our study, we investigated a headwater stream belonging to the most common stream type in 
Germany located in a nature reserve with no major anthropogenic impact except climate change. We used the most comprehensive
quantitative long‐term dataset on aquatic insects available that includes weekly measurements of species‐level insect abundance, 
daily water temperature and stream discharge as well as measurements of additional physico‐chemical variables for a 
42‐year period (1969‐2010). Overall, water temperature increased by 1.88°C and discharge patterns changed significantly. 
These changes were accompanied by an 81.6% decline in insect abundance, but an increase in richness (+8.5%), Shannon diversity (+22.7%), 
evenness (+22.4%) and interannual turnover (+34%). Moreover, the community's trophic structure and 
phenology changed: the duration of emergence increased by 15.2 days while the peak of emergence moved 13.4 days earlier. 
Additionally, we observed short‐term fluctuations (<5 years) in almost all metrics as well as complex and non‐linear responses 
of the community towards climate change that would have been missed by simply using snapshot data or shorter time series. 
Our results indicate that climate change has already altered biotic communities severely even in protected areas, 
where no other interacting stressors (pollution, habitat fragmentation etc.) are present. 
This is a striking example of the scientific value of comprehensive long‐term data in capturing the complex responses of
communities towards climate change.

Article impact statement: In a nature reserve, climate change can lead to a significant decline in insect abundance and severe 
restructuring of communities.

Dataset is a result of the long-term research project on the aquatic biota abundance and diversity ran by the Max Plank Institute of 
Limnology research station, at Schlitz, Hessen, Germany. Majority of the project data were made available by the Max Planck Society 
at their institutional repository https://pure.mpg.de/pubman/faces/ViewItemFullPage.jsp?itemId=item_969687.
Our dataset, hovewer is complemented by the data collected by the station in the last year of the project (2005-2010).
Data are optimised for running the R -code attached here. 

Attached data and code were used to produce all the main figure and all the statistical analyses in the paper. 
Data were processed as follows -original diversity and abundance data available from the Max Planck Society 
https://pure.mpg.de/pubman/faces/ViewItemFullPage.jsp?itemId=item_969687 are containinf sex-split data, 
with abundance separated for the males and females. We have merged the abundance of both sexes, 
with only one abundance value for every species available therfore. 
For the processing of the phenological data we had calculated duration for which adult specimens of every species flying (duration)
and central tendency of the emegence- i.e. the median date of the emergence period. More details are provided in the code.

R script contains additional explanations on how the data were pre-processed and should be run. 
Data from year 2006 were omitted from all analyses as sampling was not done over the entire year.
Please download scripts "code phenology of 31 spp_1410" & "community code fig1-3_JJ_VB" and the data. 
Code contains explanations regarding each dataset. Put all files in the same directory, run "community code fig1-3_JJ_VB" first.
