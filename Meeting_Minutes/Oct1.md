# October 1, 2024

## Meeting Agenda
- Variables of interest
  - Ryan paper: Gender, medication, smoking status, biopsy location, BMI (relevant to inflammation or condition)
  - Halfvarson paper: Year diagnosed, UC_extent, BMI, medication, sex, calprotein
 
- Potential research questions
  - Ryan paper:
    - How is smoking status linked to inflammation regardless of condition (CD/UC)?
    - How does the biopsy location affect the outcome of both conditions?
    - How does medication affect inflammation?

  - Halfvarson paper:
    - *Disease state (calprotein?) and date diagnosed*

- Determine sample size from above

- Finalize research question and aims
- Write out tentative workflow correlating to each research aim
  - QIIME2 quality check and R analysis
- Discuss action items, delagation of tasks, timeline, and potential other meeting time

## Meeting Notes
- Smoking and biopsy location integrate well into one another; examine community differences in smoking status, enriching inflammatory microbes
  - Compare sex: male and female
  - Check if the presence of inflammation is associated with different microbes
- Ryan dataset will be of better use
- Cannot combine with the Halfvarson dataset due to no similar variables
- Discuss confounding variables - examine literature
  - Sex
 
- Look at metadata statistics to narrow down the research question

  - Aim 1: Correlation between smoking, inflammation, and condition (compare groups based on sex)
    - DC microbiome
  - Aim 2: Community and diversity differences (richness indices, PCA plot); combine the conditional statments to compare and group into categories
    - Create new metadata file
  - Aim 3: Only pursue the significant results for further analysis with biopsy location
    - Create a diversity metrics matrix to examine significance between groups at different locations
      - Separate richness for each item in the matrix
      - Make one for each sex
  - Aim 4: Check for inflammatory microbes that are associated with each condition

## Action Items
- Make metadata for categories first
- Go through QC - create new manifest
- Create figures and split writing for the proposal
