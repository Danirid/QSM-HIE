#!/usr/bin/env python
# -*- coding: utf-8 -*-

####################################
# Filename: Plot.py
# Author: Daniel Ridani 
# Date of creation: August 6th, 2021
# Python version: 3.8.3
# 
# Description : This script plot a boxplot/violin plot for susceptibility values across different part of the brain
# 
#  
# Inputs
# ----------
# directory : 
#     Path to csv excel file that contain the data
# 
# outpath : 
#      PNG image of a boxplot/violin plot
#
####################################
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

df = pd.read_csv('ismrm_data2.csv')

df["Method"] = df["Background field removal"]+ '_' + df["Dipole inversion"]

plt.figure()

g = sns.catplot(x="Method", y="Value", palette="tab10",
          data=df, kind="box",col="Region", col_wrap=1, aspect=4) # for a violin plot type kind="violin" instead of kind="box"
sns.set(font_scale = 1.5)
g.set(ylim=(0,0.07))
sns.despine(offset=0, trim=True)
plt.grid(axis = 'y',alpha=0.2, color="black")
plt.savefig("CV_results", dpi=600)



