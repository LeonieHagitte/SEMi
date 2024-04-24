---
title: "Filter T1 observations and predictors"
output: html_document
date: "2023-11-22"
author: "A. Brandmaier"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Filter only T1 waves from longitudinal data

```{r}
df4 <- read_delim("df4.csv", delim = ";", 
    escape_double = FALSE, col_types = cols(...1 = col_skip(), 
        mergeid = col_skip(), yrbirth = col_skip()), 
    trim_ws = TRUE)

df_ <- df4 %>% dplyr::select(ends_with("_1") | starts_with("single") | starts_with("gender"))

df_ <- df_  %>% rename_with(~sub("_1$", "", .), everything())

df4 <- df_
```