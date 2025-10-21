###
# stutils.R - Subtye Data Pre-processing Utility functions
# 08/21/2025 - V1.0 - Data Preprocessing and Exploration
###

getInputData <- function(SITE_REGION_XREF,TICGEN_CLIN_DATA) {
  ###
  #' Get Input Data
  ###
  siteregionxref <- read_excel(SITE_REGION_XREF, sheet = 1) %>%
    dplyr::rename(SiteID = Site) %>%
    type_convert()
  
  ### Read from TIC Genetics file and create tcg.
  tcg <- read_excel(TICGEN_CLIN_DATA, sheet = 1) %>%
    type_convert()
  
  # get region
  tcg <- tcg %>%
    inner_join(siteregionxref) %>%
    mutate(Region = as.factor(Region)) %>%
    select(-Country)
  return(tcg)
}
### End Get Input Data

###
# combineDiseaseDiag - create uniform diagnosis for all diseases
###
combineDiseaseDiag <- function(tcg) {
  ### Disease Diagnosis Breakdowns
  tcx <- tcg %>%
    # Combine mutually exclusive diagnoses into a factor
    imap_dfc(~ if (.y %in% c(
      "Tourette DisorderID",
      "Tic Disorder-NOSID",
      "Transient Tic DisorderID",
      "Provisional Tic DisorderID",
      "Chronic Tic Disorder-Combined subtypeID",
      "Chronic Tic Disorder-Vocal subtypeID",
      "Chronic Tic Disorder-Motor subtypeID"
    )) {
      if_else(.x == 1, str_remove(.y, "ID$"), "", missing = "")
    } else{
      .x
    }) %>%
    unite(
      "TicDiag",
      c(
        "Tourette DisorderID",
        "Tic Disorder-NOSID",
        "Transient Tic DisorderID",
        "Provisional Tic DisorderID",
        "Chronic Tic Disorder-Combined subtypeID",
        "Chronic Tic Disorder-Vocal subtypeID",
        "Chronic Tic Disorder-Motor subtypeID"
      ),
      sep = ""
    ) %>%
    mutate(TicDiag = na_if(TicDiag, "") %>%
             as.factor()) %>%
    
    # Combine mutually-exclusive ADHD columns
    imap_dfc(~ if (.y %in% c(
      "ADHD, Combined TypeID",
      "ADHD, Predominantly Inattentive TypeID",
      "ADHD, Predomniantly Hyperactive-Impulsive TypeID",
      "Subclinical ADHDID"
    )) {
      if_else(.x == 1, str_remove_all(.y, "^ADHD, &ID$"), "", missing = "")
    } else{
      .x
    }) %>%
    unite(
      ADHDDiag,
      c(
        "ADHD, Combined TypeID",
        "ADHD, Predominantly Inattentive TypeID",
        "ADHD, Predomniantly Hyperactive-Impulsive TypeID",
        "Subclinical ADHDID"
      ),
      sep = ""
    ) %>%
    mutate(
      ADHDDiag = na_if(ADHDDiag, "") %>%
        replace_na("No ADHD") %>%
        as.factor() %>% relevel("No ADHD")
    ) %>%
    # Make levels of age columns unique
    mutate(across(ends_with("AgeOnset"), ~ cut(., breaks = c(0, 10, Inf)))) %>%
    imap_dfc(~ if (str_detect(.y, "AgeOnset$")) {
      .x <- as.character(.x) %>%
        replace_na(paste(.y, "NA", sep = "_")) %>%
        as.factor()
      levels(.x) <- paste(.y, c("Early", "Late", "Never"), sep = "_")
      .x
    } else{
      .x
    }) %>%
    mutate(OCDiag = OcdDiagnois, TrichDiag = TrichID)
  return(tcx)
}
### end combineDiseaseDiag

###
#' reportSplits by TIC Type / Diagnosis
###
reportSplits <- function(df) {
  print("##########################################################")
  print(paste0("Splits by Tic Type - Data Set: ", deparse(substitute(df))))
  print("##########################################################")
  print("NOTE: Diagnosis Override might switch from one type to another; Need Notes")
  deparse(substitute(df))
  # TD
  print("TIC by type")
  table(df$TicDiag)  %>%
    as.data.frame() %>% pander()
  
  print("Tic Type Override - split by type - Need Notes")
  df %>% filter(TicDxOverride == TRUE) %>%
    select(TicDiag) %>% table() %>%
    as.data.frame() %>% pander()
  
  # OCD
  print("OCD Diagnosis")
  table(df$`OCDiag`) %>%
    as.data.frame() %>% pander()
  
  print("OCD Diagnosis Override - split by type - Need Notes")
  df %>% filter(OCDxOverride == TRUE) %>%
    select(OCDiag) %>% table() %>%
    as.data.frame() %>% pander()
  
  # ADHD
  print("ADHD Diagnosis")
  table(df$ADHDDiag) %>%
    as.data.frame() %>% pander()
  
  print("ADHD Diagnosis Override - split by type")
  df %>% filter(ADHDDxOverride == TRUE) %>%
    select(ADHDDiag) %>% table() %>%
    as.data.frame() %>% pander()
  
  # Multi-Birth
  print("Multi-Birth by Type")
  table(df$MultiBirthDescr) %>%
    as.data.frame() %>% pander()
  
  print("TD Diagnosis by type - Multi Birth")
  crosstab <- table(df$MultiBirthDescr, df$TicDiag) %>%
    as.data.frame() %>% pander()
  
  print("Multi-Birth by Zygosity")
  table(df$MultiBirthTypeDescr) %>%
    as.data.frame() %>% pander()
  
  print("TD Diagnosis by type - Zygosity")
  crosstab <- table(df$MultiBirthTypeDescr, df$TicDiag) %>%
    as.data.frame() %>% pander()
  
  # AgeOnset Splits
  print("TD Diagnosis by type - AgeOnset")
  crosstab <- table(df$TicAgeOnset, df$TicDiag) %>%
    as.data.frame() %>% pander()
  
  print("OCD Diagnosis by type - AgeOnset")
  crosstab <- table(df$OCAgeOnset, df$OCDiag) %>%
    as.data.frame() %>% pander()
  
  print("Trich Diagnosis by type No/Yes/Unratable - AgeOnset")
  crosstab <- table(df$TrichAgeOnset, df$TrichDiag) %>%
    as.data.frame() %>% pander()
}
### End get splits
#################################################################
###
# Subjects with Tourettes Diagnosis - One per Family
###
getOneChildPerFamWithTD <- function(tcx) {
  tct <- tcx %>%
    filter(TicDiag == "Tourette Disorder") %>%
    filter(substr(SubjID_Display, 15, 18) > 4999) %>%
    filter(substr(SubjID_Display, 15, 18) < 6000) %>%
    arrange(SubjID_Display) %>%
    group_by(substr(SubjID_Display, 1, 13)) %>%
    slice(1) %>% ungroup()
  tct <- tct %>% select(-c("substr(SubjID_Display, 1, 13)"))
  nrow(tct)
  table(substr(tct$SubjID_Display, 15, 18))
  return(tct)
}
reportDiagOnsetMultiSplits <- function(df) {
  print("##########################################################")
  print(paste0("Diagnostics, Age Of Onset and Multi Birth splits - Data Set: ", deparse(substitute(df))))
  print("##########################################################")
  # OCD
  print("OCD Diagnosis")
  table(df$OCDiag) %>%
    as.data.frame() %>% pander()
  # ADHD
  print("ADHD Diagnosis")
  table(df$ADHDDiag) %>%
    as.data.frame() %>% pander()
  # Trich
  print("Trich Diagnosis")
  table(df$TrichDiag) %>%
    as.data.frame() %>% pander()
  
  # Diag override
  print('TIC Dx Override')
  table(df$TicDxOverride) %>%
    as.data.frame() %>% pander()
  print('OCD Dx Override')
  table(df$OCDxOverride, df$OCDiag) %>%
    as.data.frame() %>% pander()
  print('ADHD Dx Override')
  table(df$ADHDDxOverride, df$ADHDDiag) %>%
    as.data.frame() %>% pander()
  
  # Age of Onset
  print("TIC age of Onset")
  table(df$TicAgeOnset) %>%
    as.data.frame() %>% pander()
  print("OCD age of Onset")
  table(df$OCAgeOnset, df$OCDiag) %>%
    as.data.frame() %>% pander()
  print("Trich age of Onset")
  table(df$TrichAgeOnset) %>%
    as.data.frame() %>% pander()
  
  # Multi-Birth
  print('Multi-Birth / Zygosity')
  table(df$MultiBirthDescr, df$MultiBirthTypeDescr) %>%
    as.data.frame() %>% pander()
}
#################################################################
getNoOverrides <- function(df) {
  ###
  # Subjects with Tourettes Diagnosis - One per Family - No Dx Overrides
  ###
  tcno <- tcx %>%
    filter(TicDxOverride == FALSE &
             OCDxOverride == FALSE &
             ADHDDxOverride == FALSE) %>%
    filter(TicDiag == "Tourette Disorder") %>%
    filter(substr(SubjID_Display, 15, 18) > 4999) %>%
    filter(substr(SubjID_Display, 15, 18) < 6000) %>%
    arrange(SubjID_Display) %>%
    group_by(substr(SubjID_Display, 1, 13)) %>%
    slice(1) %>% ungroup()
  tcno <- tcno %>% select(-c("substr(SubjID_Display, 1, 13)"))
  nrow(tcno)
  table(substr(tcno$SubjID_Display, 15, 18))
  return(tcno)
}
##############################################################
###
# Flag count and filtering
###
getOneChildPerFamWithTDnoFlags <- function(tcno) {
  tcno %>% filter(
    `FLAG Psychosis` |
      `FLAG Other Neurological condition` |
      `FLAG Congenital anomalies` |
      `FLAG Genetic Syndrome/ Chromosomal abn`
  ) %>% nrow()
  tcno %>% filter(
    `FLAG Psychosis` |
      # `FLAG Other Neurological condition` |
      `FLAG Congenital anomalies` |
      `FLAG Genetic Syndrome/ Chromosomal abn` |
      `FLAG Other significant MEDICAL history`
  ) %>% nrow()
  tcno %>% filter(
    `FLAG Psychosis` |
      `FLAG Other Neurological condition` |
      `FLAG Congenital anomalies` |
      `FLAG Genetic Syndrome/ Chromosomal abn` |
      `FLAG Other significant PSYCHIATRIC history`
  ) %>% nrow()
  tcno %>% filter(
    `FLAG Psychosis` |
      `FLAG Other Neurological condition` |
      `FLAG Congenital anomalies` |
      `FLAG Genetic Syndrome/ Chromosomal abn` |
      `FLAG Other significant MEDICAL history` |
      `FLAG Other significant PSYCHIATRIC history`
  ) %>% nrow()
  tcno %>% filter(
    `FLAG Psychosis` |
      `FLAG Other Neurological condition` |
      `FLAG Congenital anomalies` |
      `FLAG Genetic Syndrome/ Chromosomal abn` |
      `FLAG Other significant MEDICAL history` |
      `FLAG Other significant PSYCHIATRIC history` |
      `FLAG Atypical presentation`
  ) %>% nrow()
  
  ###
  # tcnof - TS, No Dx Override, No Flags
  ###
  tcnof <- tcno %>% filter(
    !(
      `FLAG Psychosis` |
        `FLAG Other Neurological condition` |
        `FLAG Congenital anomalies` |
        `FLAG Genetic Syndrome/ Chromosomal abn` |
        `FLAG Other significant MEDICAL history` |
        `FLAG Other significant PSYCHIATRIC history` |
        `FLAG Atypical presentation`
    )
  )
  
  ###
  # Remove Unable to Rate
  ###
  tcnof <- tcnof %>%
    filter(OCDiag != "Unable to rate") %>%
    filter(TrichDiag != 9)
  
  # tcnof - stats
  print(paste0("Number of cases in tcnof: ", nrow(tcnof)))
  print(paste0("Number of parameters in tcnof: ", ncol(tcnof)))
  print(paste0("Number of probands in tcnof: ", sum(tcnof$ProbandID)))
  
  return(tcnof)
}
##############################################################
# Remove unwanted columns
# Excluding irrelevant metadata; filter out diagnosis IDs and retain descriptions, internal IDs and dates
###
genDataForClustering <- function(tcnof) { 
  tcnof <- tcnof %>%
    select(
      -c(
        GrantNum,
        DxSmrySubmitted,
        SubjectType,
        DxSmryID,
        NihGUID,
        SiteID,
        SubsiteID,
        FamilyID,
        SubjID,
        FatherID,
        MotherID,
        RUCDRNum,
        CreateDate,
        CreatedBy,
        LastUpdDate
      )
    ) %>%
    select(
      -c(
        ProbandID,
        GenderID,
        OcdDiagnoisID,
        MultiBirthID,
        MultiBirthTypeID,
        YearOfBirth,
        EvalDate,
        AgeAtEval,
        OtherGeneticStudy,
        OtherGeneticStudyDescr
      )
    ) %>%
    select(
      -c(
        EthnicityID,
        EthnicityDescr,
        RaceWhite,
        RaceBlack,
        RaceAmIndian,
        RaceAsian,
        RacePacific,
        RaceUnknown,
        RaceOther
      )
    ) %>%
    select(
      -c(
        `FLAG Atypical presentation`,
        `FLAG Psychosis`,
        `FLAG Other Neurological condition`,
        `FLAG Congenital anomalies`,
        `FLAG Genetic Syndrome/ Chromosomal abn`,
        `FLAG Other significant PSYCHIATRIC history`,
        `FLAG Other significant MEDICAL history`,
        FlagOtherDescr
      )
    ) %>%
    select(-c(TicDxOverride, OCDxOverride, ADHDDxOverride))
  # If diagnosis is true and age of onset is never, use Age @ Evaluation as proxy
  tcnof <- tcnof %>% mutate(TicAgeOnset = if_else(
    TicAgeOnset != "TicAgeOnset_Never",
    TicAgeOnset,
    if_else(Evaluated < 10, "TicAgeOnset_Early", "TicAgeOnset_Late")
  ))
  
  tcnof <- tcnof %>%
    mutate(OCAgeOnset = if_else(
      (
        OCDiag != "No OC disorder/symptoms" &
          OCAgeOnset == "OCAgeOnset_Never"
      ),
      if_else(Evaluated < 10, "OCAgeOnset_Early", "OCAgeOnset_Late"),
      OCAgeOnset
    ))
  
  tcnof <- tcnof %>%
    mutate(TrichAgeOnset = if_else(
      (TrichDiag == 1 & TrichAgeOnset == "TrichAgeOnset_Never"),
      if_else(Evaluated < 10, "TrichAgeOnset_Early", "TrichAgeOnset_Late"),
      TrichAgeOnset
    ))
  
  # Remove parameters not needed for subtyping.
  TD <- tcnof %>% select(-c(Evaluated, MultiBirthTypeDescr, TicDiag, ADHDAgeOnset))
  # Rename Fields
  TD <- TD %>%
    rename(
      Subject = SubjID_Display,
      Sex = GenderDescr,
      MultiBirth = MultiBirthDescr,
      TrichDiag = TrichDiag,
      ASDDiag = `FLAG PDD or Autism Spectrum`
    )
  # Fix levels to meaningful values for reporting based on feedback from Prof.
  TD <- TD %>%
    mutate(
      MultiBirth = fct_collapse(
        as.factor(MultiBirth),
        MultiBirth = c("Twin", "Multiple (>twin)"),
        SingleBirth = "No"
      ),
      Sex = as.factor(Sex),
      TicSxCur = if_else(TicSxCur == 1, TRUE, FALSE),
      OCSxCur = if_else(OCSxCur == 1, TRUE, FALSE),
      ADHDSxCur = if_else(ADHDSxCur == 1, TRUE, FALSE),
      TrichSxCur = if_else(TrichSxCur == 1, TRUE, FALSE),
      TicAgeOnset = as.factor(TicAgeOnset),
      OCAgeOnset = as.factor(OCAgeOnset),
      TrichAgeOnset = as.factor(TrichAgeOnset),
      OCDiag = as.factor(OCDiag),
      ADHDDiag = as.factor(ADHDDiag),
      ASDDiag = as.logical(ASDDiag),
      TrichDiag = as.logical(TrichDiag)
    ) %>%
    select(
      Subject,
      Sex,
      MultiBirth,
      TicAgeOnset,
      TicSxCur,
      ADHDDiag,
      ADHDSxCur,
      OCDiag,
      OCAgeOnset,
      OCSxCur,
      TrichDiag,
      TrichAgeOnset,
      TrichSxCur,
      ASDDiag,
      Region
    )
  return(TD)
}
##############################################################
