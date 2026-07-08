# Charting upper cervical spinal cord normative modeling and its clinical applications
# 1. Introduction
This repository provides scripts for constructing upper cervical spinal cord normative models using Generalized Additive Models for Location, Scale and Shape (GAMLSS) and for applying these models to diverse clinical tasks, including disease classification, clinical association, and prognosis prediction.

The repository includes:

Example datasets for normative model fitting, calibration, and validation.
Public normative reference packages for upper cervical spinal cord centile/Z-score calculation.
Scripts for each stage of the analytical workflow.
Example outputs for reference and validation.
# 2. Datasets
Example datasets are provided solely for demonstration and code testing purposes. They include small synthetic samples that mimic the format of the original datasets. These examples are not suitable for scientific analyses.

The full multi-center datasets used in published studies cannot be shared publicly due to data-use restrictions. Researchers interested in collaboration or access to the trained normative models should contact the authors directly via email.
# 3. Models
The public upper cervical spinal cord normative references are provided as release packages rather than raw .rds model objects. Please see the following two compressed packages:

• Upper_cervical_spinal_cord_Normative_PublicRelease.zip: public normative reference based on the original upper cervical spinl cord models.
• Upper_cervical_spinal_cord_Normative_PublicRelease_Age8_85_ICV.zip: public normative reference based on the age 8-85 ICV-adjusted upper cervical spinl cord models.
Each package contains dense age- and sex-specific centile and Z-score reference tables, a curve atlas, an example input file, and an R helper script for calculating centiles and Z-scores in new data.
> **Note:** The complete internal RDS model files are not included in this repository because they may contain fitted-model internals and subject-level information. The public release packages preserve practical normative scoring ability while excluding individual-level data, site labels, identifiers, and explicit sample-size fields.
# 4. Scripts
## 4.1 Required R packages and installation
Check-and-install-packages.R Automatically checks and installs all necessary R packages required for modeling, application, and visualization.
## 4.2 Normative curve estimation and peak age determination
Normative-model-fit.R Fits normative lifespan models using GAMLSS based on upper cervical spinal cord measures, estimates peak ages, and generates normative trajectories.
## 4.3 Bootstrap analysis of normative curve
## 4.4 ICV adjustment and model comparison
## 4.5 Model calibration using new datasets
## 4.6 Applying normative models to disease dataset
## 4.7 Applying normative models to individual-level data
## 4.8 Statistical analyses of deviation scores across diseases
## 4.9 Clinical downstream tasks
# 5. License
MIT License

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
