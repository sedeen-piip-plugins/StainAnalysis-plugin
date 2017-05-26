<h3 align="center">Stain Analysis Plugin</h3>
This plugin calculates the percentage area of a region of interest which has been stained using an IHC biomarker. The color deconvolution method described by Ruifrok <sup>[1]</sup> is used to separate the color channels. The stain of interest is selected and a threshold image is then obtained by applying a linear thresholding algorithm on the deconvolved image. Finally, the percentage of pixels that have been stained is calculated using the threshold image. All the processing is done in the highest resolution.

## User Manual
1.	Open the WSI image and define the region of interest (ROI). 
2.	Load the “StainAnalysis” plugin from the pulldown list of Algorithms (Fig 1.)

<div align="right">
  <img src="https://github.com/sedeen-piip-plugins/StainAnalysis-plugin/blob/master/Images/StainAnalysis_1_1.png"/>
</div>

<div align="right">
<strong>Fig1.</strong> 
  <font size="2" color="black">Select the Stain Analysis Plugin.</font>
</div>

3.	Select the stain using the Selected Stain pulldown menu (Fig 2.). The plugin provides a number of "built in" stain vectors:   Hematoxylin and Eosin (H&E), Hematoxylin and DAB (H DAB), and Hematoxylin, Eosin and DAB (H&E DAB) [1] .

<div align="right">
  <img src="https://github.com/sedeen-piip-plugins/StainAnalysis-plugin/blob/master/Images/StainAnalysis_1_2.png"/>
</div>

<div align="right">
<strong>Fig2.</strong> 
  <font size="2" color="black">Select stain option.</font>
</div>

4.	Users can also determine their own vectors to achieve an accurate stain separation. This option is available in two ways: “From ROI” or “Load From File”. If the “From ROI” option is set, users also need to choose three ROIs to compute the stain vectors. Select small ROIs areas which are intensely stained with only one of the stain, without empty background. If the staining method uses only 2 colors instead of 3, for the 3rd selection just select a small ROI from the background. The plug-in will compute the stain vectors based on three regions [1] 
The computed stain vectors will be saved in default directory (the Sedeen folder in the same directory as the image): "Path to the image/sedeen/" as “StainsFile.csv” file for future use. The stain vectors can be provided by “.csv” file by selecting “Load From File” option. An example of the “.csv” file will be provided in the plug-in directory. This file needs to be copied in default directory: “Path to the image/sedeen/”.
5.	The stain to be used for the analysis should be selected from the Display pull down menu.
6.	If desired, modify the Threshold value which is in the range 0.0 to 50.0
7.	The Processing ROI option allows the % staining to be calculated over a user defined ROI. (Fig 3). If an ROI is not selected then the calculation will be done over the displayed region.

<div align="right">
  <img src="https://github.com/sedeen-piip-plugins/StainAnalysis-plugin/blob/master/Images/StainAnalysis_1_3.png"/>
</div>

<div align="right">
<strong>Fig3.</strong> 
  <font size="2" color="black">Select the Processing ROI.</font>
</div>

8.	Clicking on the Run button will execute the algorithm. The calculated % area is shown in the results panel. A color overlay shows which pixels are included in the % stain calculation. Use the Show Result checkbox to toggle the overlay on and off. The overlay is generated at the resolution of the displayed image but the calculation is carried out at full image resolution.

<div align="right">
  <img src="https://github.com/sedeen-piip-plugins/StainAnalysis-plugin/blob/master/Images/StainAnalysis_1_4.png"/>
</div>

<div align="right">
<strong>Fig4.</strong> 
  <font size="2" color="black">Select the Processing ROI.</font>
</div>

<sup>[1]</sup> A. C. Ruifrok and D. A. Johnston, “Quantification of histochemical staining by color deconvolution,” Anal. Quant. Cytol. Histol., vol. 23, no. 4, pp. 291–299, 2001.

## Authors
**Azadeh Yazanpanah**

Stain Analysis Plugin has been developed by Martel lab at Sunnybrook Research Institute (SRI), University of Toronto.
[Funding provided by NIH.](https://itcr.nci.nih.gov/funded-project/pathology-image-informatics-platform-visualization-analysis-and-management)

## Copyright & License
