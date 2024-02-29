library(tidyverse)
library(magrittr)
library(haven)
library(grf)
library(car)
library(cowplot)
library(stargazer)
source("helpers.R")

covariates <- c('hhh_gender', 'program_knowledge', 'hhh_age', 'monthly_spending', 'any_drops', 'maths_results', 'hhh_literacy', 'age', 'gender', 'hhid', 'schoolid', 'group', 'group_cond', 'benef')

X_vars <- c('hhh_gender', 'program_knowledge', 'hhh_age', 'monthly_spending', 'hhh_literacy', 'age', 'gender', 'benef')


maths <- read_dta("LCT Morroco/cct_aser_an.dta")
baseline <- read_dta("LCT Morroco/cct_baseline_an.dta")
endline <- read_dta("LCT Morroco/cct_endline_an.dta")
school_baseline <- read_dta("LCT Morroco/cct_preliminary_survey_an.dta")
# dropout <- read_dta('LCT Morroco/cct_school_visits_an.dta')
admin_data <- read_dta('LCT Morroco/cct_tayssir_admin_data_an.dta')
knowledge <- read_dta("LCT Morroco/cct_knowledge_households_year1_an.dta")

# Identify individual maths respondents
maths_link <- endline %>%
  filter(!is.na(hhid_endline) & hhid_endline != "") %>%
  labelled::look_for_and_select("Child's ID") %>%
  select(starts_with("d")) %>%
  bind_cols(endline["hhid_endline"]) %>%
  gather(key = "Variable", value = "Name", -hhid_endline) %>%
  mutate(child_num = str_sub(Variable, -1, -1)) %>%
  left_join(maths, by = "hhid_endline") %>%
  filter(Name == id_enf_test)

# Link to baseline data
# Get the child ids in baseline
selected <- baseline %>%
  labelled::look_for_and_select("Member ID") %>%
  select(starts_with("d")) %>%
  bind_cols(baseline["hhid"]) %>%
  gather(key = "Variable", value = "Name", -hhid) %>%
  mutate(child_num = str_sub(Variable, -1, -1)) %>%
  left_join(maths_link, by = c("hhid", "Name")) %>%
  filter(Name == id_enf_test) %>%
  select(hhid, child_num.x, Name) %>%
  mutate(child_num.x = as.numeric(child_num.x)) %>%
  group_by(hhid) %>%
  slice(1) %>%
  ungroup()

names(selected) <- c("hhid", "child_num", "member_id")
baseline1 <- left_join(selected, baseline, by = "hhid")

# Get section a self
baseline <- baseline1 %>%
  get_vars_by_regex("^a", "_[:digit:]+", "member_id") %>%
  get_vars_by_regex("^c[[:digit:]]", "[:digit:]+", "child_num") %>%
  get_vars_by_regex("^d[[:digit:]]", "_[:digit:]+", "child_num") %>%
  get_vars_by_regex("^e[[:digit:]]", "_[:digit:]+", "child_num")


joined <- baseline %>%
  left_join(maths, by = "hhid", suffix = c("", ".y")) %>%
  select_at(
    vars(-ends_with(".y"))
  )

baseline %<>%
  mutate(schoolunitid = factor(schoolunitid),
         group_control = as.numeric(group != 0),
         group = car::recode(group, "0 = NA;"),
         group_cond = car::recode(group, "1=1; else=2")
  )

week_to_month <- function (x) {
  (x / 7) * 30
}

year_to_month <- function (x) {
  x / 12
}


expenditure <- baseline %>%
  select(starts_with('h1_'),
         starts_with('h3_'),
         starts_with('h4_')
  ) %>%
  select_if(is.numeric) %>%
  mutate_at(
    vars(starts_with('h1_')),
    week_to_month
  ) %>%
  mutate_at(
    vars(starts_with('h4_')),
    week_to_month
  )

expenditure[expenditure < 0] <- 0

baseline$monthly_spending <- rowSums(expenditure, na.rm = TRUE) %>% log()

num_kids <- baseline %>%
  select(hhid, matches("^d1_[0-9]{1,2}")) %>%
  gather(-hhid, key = "question", value = "value") %>%
  mutate(hhid = hhid %>% as_factor() %>% as.character()) %>%
  group_by(hhid) %>%
  dplyr::summarise(est_num_kids = sum(value == 1, na.rm = T))


assets <- baseline %>%
  select(b10_5, b10_6, b10_20, b10_26,  b10_21, b10_28, b4_1, b5_1, b9_5, monthly_spending) %>%
  mutate(b4_1 = car::recode(b4_1, "1:3 = 1; 4:5 = 2; -99 = NA"),
         b5_1 = car::recode(b5_1, "1:3 = 1; 5 = 1; 4 = 2; 6 = 2; -99 = NA"),
         b9_5 = car::recode(b9_5, "1 = 1; 2 = 2; -99:-1 = NA")
  )

assets[assets==-99]<-NA

maths_results <- rowSums(maths[c('t6', 't5_1', 't4_1', 't3_1')])

knowledge1 <- knowledge %>%
  transmute(
    hhid = hhid,
    program_knowledge = ksm_c4_1 == 6
  )

joined1 <- baseline[c('hhid', 'group', 'group_cond', 'group_control', 'monthly_spending', 'schoolid')] %>%
  left_join(num_kids, by = "hhid") %>%
  left_join(knowledge1, by = "hhid") %>%
  filter(!duplicated(hhid))

joined2 <- baseline %>%
  left_join(num_kids, by = "hhid") %>%
  left_join(knowledge1, by = "hhid") %>%
  filter(!duplicated(hhid))


maths %<>%
  transmute(maths_results = as.numeric(t6 == 1) +
           (t5_2 / 2) +
           (t4_2 / 5) +
           (t3_2 / 5),
           hhid = hhid
         )

joined1 %<>%
  left_join(maths, by = "hhid")

joined2 %<>%
  left_join(maths, by = "hhid")

# Add endline ID to be able to join to maths data
joined1 %<>%
  dplyr::inner_join(endline[c('hhid', 'hhid_endline')], by = 'hhid', suffix = c("", ".y")) %>%
  left_join(admin_data[c('hhid_endline', 'stud_id_tayssir', 'datenaiseleve', 'genre_el')], by = 'hhid_endline', suffix = c("", ".y")) %>%
  select_at(
    vars(-ends_with(".y"))
  )

joined2 %<>%
  dplyr::inner_join(endline[c('hhid', 'hhid_endline')], by = 'hhid', suffix = c("", ".y")) %>%
  left_join(admin_data[c('hhid_endline', 'stud_id_tayssir', 'datenaiseleve', 'genre_el')], by = 'hhid_endline', suffix = c("", ".y")) %>%
  select_at(
    vars(-ends_with(".y"))
  )


joined_num2 <- joined2 %>%
  mutate(benef = as_factor(benef)) %>%
  mutate_at(vars(-hhid, -schoolunitid), as.numeric)

joined_num2[joined_num2 < 0] <- NA

joined_num2 <- joined_num2[complete.cases(joined_num2[c("maths_results")]),]



joined_num2 %>% write_rds('morocco2.rds')

