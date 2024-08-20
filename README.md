# Nixtamalization Moisture Content Prediction in a Maize Inbred Diversity Panel
Michael Burns <br> University of Minnesota - Applied Plant Sciences

## Synopsis:
Maize is a critical crop in the United States, primarily used for livestock feed and ethanol production rather than direct 
human consumption. Consequently, research into food-grade maize varieties has been limited. Food-grade maize is typically 
selected from existing No. 2 yellow dent germplasm. This narrow selection limits the ability to enhance specific traits 
crucial for food products, such as moisture content during nixtamalization, a traditional cooking process used to prepare masa. 
Nixtamalization involves cooking maize kernels in an alkaline solution to remove the pericarp and soften the grain for grinding. 
The moisture content of the resulting masa affects its taste, texture oil content, and more making it an important trait for 
food-grade maize.

The nixtamalization process impacts moisture content through interactions with various kernel macromolecules, including fiber, 
protein, and starch. These components interact differently with the alkaline solution, affecting how much moisture the grain 
absorbs. The macromolecular composition of the pericarp, germ, and endosperm play significant roles in determining nixtamalization 
moisture content. Multiple methods have been proposed for assessing nixtamalization moisture content, however they tend to focus
on the manufacturer and require anywhere from 100g to 1,000lbs of seed. This is not a viable option in early maize breeding pipelines
where large numbers of small-quantity samples are present. To develop a breeding platform, a higher throughput process is needed.

To address this, we utilized Near-Infrared (NIR) spectroscopy combined with machine learning to develop a high-throughput method 
for predicting moisture content during nixtamalization. By training a machine learning model on NIR spectra and moisture content 
data we were able to train an accurate (Spearman R = 0.85) on only 30g of grain per sample. This model was then applied to a larger 
dataset, including 501 genotypes across five environments, allowing for a comprehensive analysis of trait variation and genetic factors. 
The study's approach offers a rapid and scalable solution for evaluating nixtamalization moisture content, providing valuable 
insights for improving food-grade maize germplasm and informing breeding strategies for masa-based products.

## Publication:
To learn more about the utility and biological insights of this work, please see:

Burns et al. (2021). Predicting moisture content during maize nixtamalization using machine learning with NIR spectroscopy.
Theoretical and Applied Genetics. 134(11):3743-3757.
doi: 10.1007/s00122-021-03926-8.

## Repository Files:
- Manuscript_Code
  - All code used in the creation and analysis of the machine learning model to predict moisture content after nixtamalization.
- Data
  - All relevant data used to create the machine learning model and created from the analysis of the machine learning model.
 
