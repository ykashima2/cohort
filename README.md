

# 1. Required packages

```{r }
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot)
library(patchwork)
```


# 2. Input file

input data doenload from below.

```
wget https://kero.hgc.jp/cgi-bin/download/Cohort_Multiome/Analysis_input_dataKCohrot.tar.gz
Analysis_input_dataKCohrot.tar.gz
tar zxvf Analysis_input_dataKCohrot.tar.gz
mv  Analysis_input_dataKCohrot data
mkdir output
```

# 3. Analysis step

## 3-1. Figure 2  Corellation analysis

### Corellation analysis

```{bash}
Rscript script/correlation_Fig2.r
```

### Cross varidation like cheke

```{bash}
Rscript script/correlation_Fig2_CrossVaridation_sex.r
```

### Femake &  Male Corellation

```{bash}
Rscript script/FemakeMale_DiffDetection.r
```

## 3-2. Fiure 3 age assotiation

```{bash}
Rscript script/LoessSmoothing_Age_wrapper.r
```



