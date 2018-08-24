# Error-Analysis
This repository contains Matlab functions for the evaluation of closeness of two datasets. The datasets can be inputted as either time series or stack of images.
Given one or multiple target time series and a same number of corresponding reference time series, use the "comp_TS.m" function can compute statistics and error metrics of and between the target and the reference time seies.
Given a stack of 2-D images and a number of reference time series on locations within the coordinate of the 2-D images, use the "RS2TS.m" function first to extract target the time series of the locations of interest; then use the "comp_TS.m" function to compute statistics and error metrics of and between the target and the reference time seies.
Given a stack of target 2-D images and also a stack of reference images, use the "comp_RS.m" function can compute statistics and error metrics of and between the target and the reference images.
