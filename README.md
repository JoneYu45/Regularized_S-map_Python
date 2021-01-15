# Regularized_S-map_Python

We develop this python script using the Scikit-learn’s ElasticNetCV function. The regularized S-map is used to analyze the varying bacterial interaction in community. Please cite Yu, Z., Gan, Z., Huang, H., Zhu, Y. & Meng, F. Regularized S-Map Reveals Varying Bacterial Interactions. Appl. Environ. Microbiol. 86, e01615-01620, doi:10.1128/aem.01615-20 (2020) if this method is useful for your research.

Detailed info about the regularized S-map can be found in the previous studies.

# Requirement

Linux version 5.4.0-45-generic (buildd@lcy01-amd64-014) (gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04)) #49~18.04.2-Ubuntu SMP Wed Aug 26 16:29:02 UTC 2020

Python==3.6.5

# Instruction

**How to prepare the input time-series data? **

Please refer to the csv files in the Input folder. Three demo files are provided. In each input file, the columns represent OTU numbers, the rows represent time points, and the values represent the relative abundance of an OTU at a time point.

**How to create the Jacobian matrix? **

(Old version) After preparing your input time-series data, enter the Function folder and change the input directory in the Regularized_S-map.py (Line 88) accordingly.
Then run 

***python Regularized_S-map.py ***

(New version) After preparing your input time-series data, run

***python Regularized_S-map_new.py -I ../Input/intersect_OTU_0.csv -O ../Output/0 ***

**How to evulate the Jacobian matrix inference?**

The Jacobian matrices and quality data will be generated in the Output folder.In the Output folder, there are two subfolder named coef and fit_result, respectively. The coef folder contains the Jacobian matrices and the fit_result folder contains the quality data. In order to evaluate the inference quality, run the Inference_quality_control.R in Function folder.

You can change the input path (line 5) and directory (6) in the Inference_quality_control.R if you save your results or name your data differently.

Two plots (RMSE and RMSE/STD boxplot) will be generated for each sample and can help you to determine the best theta to use. Commonly speaking,  the smaller RMSE the better quality and the RMSE/STD of a good inference should be lower than 1.

The Jacobian element file is named like 0_0.1_coefs.csv, where the first digital number (0) repressnts the target OTU and the second digital number (0.1) represents the theta used for inference. This csv table provides the info about the effects of other OTUs (column) on target OTU at each time point (row).

# References

Zou H, Hastie T. 2005. Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67:301-320.

Pedregosa et al., 2011, Scikit-learn: Machine Learning in Python,JMLR 12, pp. 2825-2830.

Deyle ER, May RM, Munch SB, Sugihara G. 2016. Tracking and forecasting ecosystem interactions in real time. Proceedings Biological Sciences 283:20152258.

Cenci S SS. 2019. Non-parametric estimation of the structural stability of non-equilibrium community dynamics. Nature ecology & evolution 3:912.

Yu, Z., Gan, Z., Huang, H., Zhu, Y. & Meng, F. Regularized S-Map Reveals Varying Bacterial Interactions. Appl. Environ. Microbiol. 86, e01615-01620, doi:10.1128/aem.01615-20 (2020). **(This is my paper! Would love to hear your opinions on it!)**
