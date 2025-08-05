Mental and Physical Health Multimorbidity: Evidence of Neural, Environmental, and Genetic Signatures of a 'Disease' Factor
================================================================================================================================================

This repository contains the scripts and resources used in our study. The project is structured into multiple directories to reflect the different types of analyses described in the manuscript.

* * *

**Repository Structure**
------------------------

### 1. **CFA**

* **Description**: Contains R scripts for confirmatory factor analysis (CFA), including:
  * Unifactor model
  * Bifactor model
  * Corfactor model
* **Tools Used**:
  * `Lavaan` (4.4.3)
  * `Psych` (4.3.3)
  * `SemTools` (0.5-6)

### 2. **ESEM**

* **Description**: Contains MPLUS scripts for Exploratory Structural Equation Modeling (ESEM). Includes the implementation of ESEM-Bifactor models.
* **Tools Used**:
  * Mplus7
  * `BifactorIndicesCalculator` (0.2.2)

### 3. **Mediation Effects Model**

* **Description**: Scripts for analyzing:
  * Extreme deviations
  * Correlations with the 'disease' factor (`d-factor`)
  * Mediation effects
* **Tools Used**:
  * `Pyprocessmacro` (1.0.12)
  * `Statsmodels` (0.14.0)
  * `Pandas` (2.0.3)

### 4. **Normative Model**

* **Description**: Contains Python scripts to develop the normative model used in this study.  
  This implementation is based on the repository from [Predictive Clinical Neuroscience](https://github.com/predictive-clinical-neuroscience/NM_educational_OHBM24).
* **Tools Used**:
  * `Pcntoolkit` (0.28)
  * `Scikit-learn` (1.3.0)
  * `Numpy` (1.25.0)
  * `Pandas` (2.0.3)

* * *

**Software and Tools**
----------------------

### **Python Packages**

* `Numpy` (1.25.0)
* `Pandas` (2.0.3)
* `Scikit-learn` (1.3.0)
* `Matplotlib` (3.7.2)
* `Scipy` (1.11.1)
* `Statsmodels` (0.14.0)
* `Pcntoolkit` (0.28)
* `Seaborn` (0.12.2)
* `Pyprocessmacro` (1.0.12)
* `Gwaslab` (3.4.48)

### **Matlab Toolboxes**

* BrainNet (1.7)
* NIfTI (1.27)

### **R Packages**

* `Lavaan` (4.4.3)
* `Psych` (4.3.3)
* `SemTools` (0.5-6)
* `BifactorIndicesCalculator` (0.2.2)

### **Other Software**

* Mplus7
* fastGWA

* * *

**Usage Instructions**
----------------------

**Important**: The final version of the cleaned and validated scripts will be uploaded prior to the manuscript's publication.

* * *

**Citation**
------------

If you use any scripts or methods from this repository, please cite the associated paper:  
**Mental and Physical Health Multimorbidity: Evidence of Neural, Environmental, and Genetic Signatures of a 'Disease' Factor**
