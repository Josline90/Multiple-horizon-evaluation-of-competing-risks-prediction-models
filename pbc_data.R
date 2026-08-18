# Data preparation
#--------------------------------libraries
pacman::p_load(dplyr, data.table, tidyr, janitor)

#-----Load PBC (liver disease) dataset: CR data; 0-censored, 1-transplantation, 2-death
data(pbc, package="survival")
glimpse(pbc); str(pbc)

cols <- c("chol", "trig", "platelet", "copper")
pbc[cols] <- lapply(pbc[cols], as.numeric)

setDT(pbc); dim(pbc) 

colSums(is.na(pbc)) # Any missingness
pbc_312 <- head(pbc, 312); dim(pbc_312) # participants assigned treatment (N = 312)
missing_prop <- pbc_312 %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing_Count") %>%
  mutate(
    Missing_Percent = round(Missing_Count / nrow(pbc_312) * 100, 2),
    Missing = paste0(Missing_Count, " (", Missing_Percent, "%)")
  ) %>%
  select(Variable, Missing)
missing_prop

pbc_276 <- na.omit(pbc_312); dim(pbc_276); names(pbc_276) #---remove missing

pbc_276$stage <- factor(
  pbc_276$stage, levels = c(1, 2, 3, 4), ordered = TRUE
)
pbc_276$edema <- factor(
  pbc_276$edema, levels = c(0, 0.5, 1), ordered = TRUE
)

#------ some summaries
tbl <- pbc_276
#--categorical
tbl %>%
  tabyl(status, trt) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting() %>%
  adorn_ns()
#---continuous
tbl %>%
  group_by(status) %>%
  summarise(
    median_age = median(age, na.rm = TRUE),
    Q1_age = quantile(age, 0.25, na.rm = TRUE),
    Q3_age = quantile(age, 0.75, na.rm = TRUE)
  )
table(tbl$status)

#---some recoding
pbc_276[, `:=`(trt = +(trt == 1), sex = +(sex == "m"))]
#pbc_276[, time:=time/365.25] #--- time in years
pbc_276[, time:=time/30.44] # in months

# Evaluation time points
pr <- c(.25, .5, .75) 
q_pbc <- quantile(pbc_276$time[pbc_276$status == 2], probs = pr, na.rm = TRUE, names = TRUE)
q_pbc # months 22, 39, 75-max 150

dat <- pbc_276 
dat <- dat[, !names(dat) %in% c("id", "ascites"), with = FALSE]; dim(dat)

# Convert categorical variables to factors.
factor_vars <- c("trt", "sex", "hepato", "spiders")
factor_vars <- intersect(factor_vars, names(dat))
dat[, (factor_vars) := lapply(.SD, factor), .SDcols = factor_vars]

 # ! ordered factor variables problematic with RSF !
 ordered_cols <- names(dat)[sapply(dat, is.ordered)]
 if (length(ordered_cols) > 0) {
   dat[, (ordered_cols) := lapply(.SD, function(x) factor(as.character(x), ordered = FALSE)),
       .SDcols = ordered_cols]
 }

str(dat)


