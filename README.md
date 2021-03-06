# Regularized_S-map_Python

We develop this python script using the Scikit-learn’s ElasticNetCV function. The regularized S-map is used to analyze the varying bacterial interaction in community. Please cite Yu, Z., Gan, Z., Huang, H., Zhu, Y. & Meng, F. Regularized S-Map Reveals Varying Bacterial Interactions. Appl. Environ. Microbiol. 86, e01615-01620, doi:10.1128/aem.01615-20 (2020) if this method is useful for your research.

Detailed info about the regularized S-map can be found in the previous studies.

# Requirement

Linux version 5.4.0-45-generic (buildd@lcy01-amd64-014) (gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1 18.04)) #49~18.04.2-Ubuntu SMP Wed Aug 26 16:29:02 UTC 2020

Python==3.6.5

# Instruction

## How to prepare the input time-series data?

Please refer to the csv files in the Input folder. Three demo files are provided. In each input file, the columns represent OTU numbers, the rows represent time points, and the values represent the relative abundance of an OTU at a time point.

Commonly, we inferer the Jacobian matrices using the continuous time series. However, if you don't have enough time points but have the replicate time series from independet reactors. You can combine the replicate OTU tables as the input. However, you should change the uncontinuous in line 127 and specify the uncontinuous states list in the input in line 58 of Regularized_S-map_new.py in the meantime.

## How to create the Jacobian matrix?

After preparing your input time-series data, run

***python Regularized_S-map_new.py -I ../Input/intersect_OTU_0.csv -O ../Output/0***

The Regularized_S-map_new.py is more user-friendly than the Regularized_S-map.py. Just ignore the Regularized_S-map.py file.

## How to evulate the Jacobian matrix inference?

The Jacobian matrices and quality data will be generated in the Output folder.In the Output folder, there are two subfolder named coef and fit_result, respectively. The coef folder contains the Jacobian matrices and the fit_result folder contains the quality data. In order to evaluate the inference quality, run the Inference_quality_control.R in Function folder.

You can change the input path (line 5) and directory (i.e. sample ID; line 6) in the Inference_quality_control.R if you save your results or name your data differently.

Two plots (RMSE and RMSE/STD boxplots) will be generated for each sample and can help you to determine the best theta to use. Commonly speaking,  the smaller RMSE the better quality and the RMSE/STD of a good inference should be lower than 1.

The Jacobian element file is named like 0_0.1_coefs.csv, where the first digital number (0) repressnts the target OTU and the second digital number (0.1) represents the theta used for inference. This csv table provides the info about the effects of other OTUs (column) on target OTU at each time point (row).

## The final Jacobian elements after the quality control?

After the inference quality control, a csv table named XX_best_coefs.csv will be generated. This table can provided the best inference of the effects of other OTUs (column) on target OTU at each time point (row). The row names looks like 0-1, where the first digital number represents the target OTUs, and the second digital number represents the time point.

You can use these Jacobian elements to calculate the volume contraction rate  andharmony level according to the Cenci S et al. and Zhong Y et al. (See references for further info).

## What are the good values for some parameters used in the inference?

Several parameters are used for the inference, including the l_grid, iteration, CV, test_len, dominate_threshold, zero_frequency_threshold, and theta.

**l_grid** (line 123): This represents the lambda λ gradient. It controls the penalty in elastic net regression.

**iteration** (line 124): This represents the iteration times for the elastic net regression.

**CV** (line 125): This repressnts the K-fold cross validation. Commonly, 10 fold cross validation is used.

**test_len** (line 126): This represents the number of test set. We only randomly picked 5 observations as the test set. You can use as many as you want, but make sure enough observations are used for the training.

**dominate_threshold** (line 128): OTUs whose abundances are more than this threshold will be included in the inference. This parameter will be used when the select in line 130 equals True.

**zero_frequency_threshold** (line 129): OTUs whose absence frequency are less than this threshold will be included in the inference. This parameter will be used when the select in line 130 equals True.

**theta** (line 152): The theta list used for the inference. You can specify the thetas (θ) list you want to try. Only the states closer to the target state will be used for regression when θ is large. Try different theta and find out how it affects the inference quality.

# References

Zou H, Hastie T. 2005. Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67:301-320.

Pedregosa et al., 2011, Scikit-learn: Machine Learning in Python,JMLR 12, pp. 2825-2830.

Deyle ER, May RM, Munch SB, Sugihara G. 2016. Tracking and forecasting ecosystem interactions in real time. Proceedings Biological Sciences 283:20152258.

Cenci S SS. 2019. Non-parametric estimation of the structural stability of non-equilibrium community dynamics. Nature ecology & evolution 3:912.

Yu, Z., Gan, Z., Huang, H., Zhu, Y. & Meng, F. Regularized S-Map Reveals Varying Bacterial Interactions. Appl. Environ. Microbiol. 86, e01615-01620, doi:10.1128/aem.01615-20 (2020). **(This is my paper! Would love to hear your opinions on it!)**
