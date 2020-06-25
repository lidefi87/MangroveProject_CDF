Mangrove fish assemblages from the Galapagos Islands
================
Denisse Fierro Arcos
2020-06-11

## Data analysis: Mangrove fish assemblages from the Galapagos Islands

Date of creation: 2018-04-30  
Date of latest update: 2020-03-03  
Version \# 3 Prepared for the Sharks Ecology Project of the Charles
Darwin Foundation (CDF).  
Related to work published in xxxx titled “First archipelago-wide
characterisation and analysis of spatial distribution patterns of
mangrove fish assemblages in the Galapagos, a UNESCO World Heritage
Site”.

### Uploading relevant libraries

``` r
library(tidyverse)
library(vegan)
library(statsr)
library(patchwork)
library(ggpubr)
library(cowplot)
library(ggrepel)
```

### Calling Kernel Density Function developed by Langlois et al, 2012

``` r
source(file = "Modified_Langlois_etal_KDE.R")
```

### Uploading general data to be used in the analysis

``` r
#Fish database: includes information information such as max length, and trophic level
FishData <- read.csv("Data/FishDB.csv")

#List of non-fish species sampled by both sampling methods
NonFish = c("Cardisoma crassum", "Chelonia mydas", "Chelonia  mydas", "Cheloniidae sp", "Phalacrocorax harrisi", "Spheniscus mendiculus", "Scyllarides astori", "Zalophus wollebaeki")
```

### Underwater Visual Census (UVC) data

``` r
#Count data
UVC <- read.csv("Data/UVC_Data_AllPts.csv")

#UVC environmental data
UVCEnviro <- read.csv("Data/UVCEnviro_AllPts.csv")
```

#### UVC Data transformation

``` r
#Filtering out non-fish species from UVC database and creating tibble table from UVC data frame containing counts
UVC <- UVC %>% tbl_df() %>% filter(!Species %in% NonFish) %>% 
  mutate(Species = factor(Species))

#Filtering out non-mangrove areas and creating a tibble table from UVC environmental data
UVCEnviro <- UVCEnviro %>% filter(Habitat == "Mangrove") %>% tbl_df()

#Creating new tibble table containing counts data per site and per species
CountsUVC <- UVC %>% 
  right_join(., select(UVCEnviro, SiteCode), by = "SiteCode") %>%
  mutate(SiteName = factor(paste(SiteCode, "U", sep = "_")), 
         Species = factor(Species)) %>% 
  group_by(Bioregion, Island, SiteName, Species) %>% 
  summarise_at(vars(Counts), funs(N = sum(., na.rm = T))) %>% 
  complete(Species, nesting(Bioregion, SiteName))
head(CountsUVC, n = 5)
```

    ## # A tibble: 5 x 5
    ## # Groups:   Bioregion, Island, SiteName [1]
    ##   Island   Species                 Bioregion SiteName     N
    ##   <chr>    <fct>                   <chr>     <fct>    <int>
    ## 1 Floreana Abudefduf concolor      CSE       BAR1_1_U    NA
    ## 2 Floreana Abudefduf troschelii    CSE       BAR1_1_U    63
    ## 3 Floreana Aetobatus ocellatus     CSE       BAR1_1_U    NA
    ## 4 Floreana Anisotremus interruptus CSE       BAR1_1_U    NA
    ## 5 Floreana Apogon atradorsatus     CSE       BAR1_1_U    NA

``` r
#Checking species are correctly named as per FishData and attaching Family and Genus data
CountsUVC <- FishData %>% select(Family, Genus, ScientificName, ValidName, TrophicCat) %>% 
  rename("Species"="ScientificName") %>% 
  left_join(CountsUVC, ., by = "Species") %>% 
  select(-Species) %>% 
  rename("Species"="ValidName") %>% 
  mutate(Species = factor(Species)) %>% 
  ungroup()
```

### Baited Remote Underwater Video (BRUV) data

``` r
#MaxN data
BRUV <- read.csv("Data/BRUVSMaxN.csv")
#Environmental data
BRUVEnviro <- read.csv("Data/BRUVEnviro.csv")
```

#### BRUV Data transformation

``` r
#Filtering out non-fish species from BRUV database, removing MaxN for different life stages/sex (MaxN column will conserve actual MaxN observed) and creating tibble table from BRUV data frame containing counts
BRUV <- BRUV %>% tbl_df() %>%
  filter(!SpeciesName %in% NonFish) %>% 
  filter(.,!duplicated(select(.,SiteName, Family, Genus, SpeciesName)))%>% 
  mutate(SpeciesName = factor(SpeciesName))

#Filtering out non-mangrove areas and creating a tibble table from BRUV environmental data
BRUVEnviro <- BRUVEnviro %>% filter(Habitat == "Mangrove") %>% tbl_df()

#Creating new tibble table containing counts data per site and per species
CountsBRUV <- BRUV %>% select(-c(Species, Stage, MaxN.by.stage)) %>% 
  right_join(., select(BRUVEnviro, SiteName, Bioregion, Island), by = "SiteName") %>%
  mutate(SiteName = factor(paste(SiteCode, Replicate, "B", sep = "_")), SpeciesName = factor(SpeciesName)) %>% 
  group_by(Bioregion, Island, SiteName, SpeciesName) %>% 
  summarise_at(vars(MaxN), funs(N = sum(., na.rm = T))) %>% 
  rename("Species"="SpeciesName") %>% 
  complete(Species, nesting(Bioregion, SiteName)) 
head(CountsUVC, n = 5)
```

    ## # A tibble: 5 x 8
    ##   Island  Bioregion SiteName     N Family     Genus    Species        TrophicCat
    ##   <chr>   <chr>     <fct>    <int> <chr>      <chr>    <fct>          <chr>     
    ## 1 Florea~ CSE       BAR1_1_U    NA Pomacentr~ Abudefd~ Abudefduf con~ Omnivore  
    ## 2 Florea~ CSE       BAR1_1_U    63 Pomacentr~ Abudefd~ Abudefduf tro~ Carnivore 
    ## 3 Florea~ CSE       BAR1_1_U    NA Myliobati~ Aetobat~ Aetobatus oce~ Carnivore 
    ## 4 Florea~ CSE       BAR1_1_U    NA Haemulidae Anisotr~ Anisotremus i~ Carnivore 
    ## 5 Florea~ CSE       BAR1_1_U    NA Apogonidae Apogon   Apogon atrado~ Carnivore

``` r
#Checking species are correctly named as per FishData and attaching Family and Genus data
CountsBRUV <- FishData %>% select(Family, Genus, ScientificName, ValidName, TrophicCat) %>% 
  rename("Species"="ScientificName") %>% 
  left_join(CountsBRUV, ., by = "Species") %>% 
  select(-Species) %>% 
  rename("Species"="ValidName") %>% 
  mutate(Species = factor(Species)) %>% 
  ungroup()
```

### Fish counts statistics both methods

``` r
#Creating dataset containing information for both methods
CountsALL <- rbind(CountsBRUV %>% mutate(Method = "BRUV"), 
                   CountsUVC %>% mutate(Method = "UVC"))
```

#### Families

``` r
BasicStats <- function(df, TaxLvl){
  if(TaxLvl == "Family"){
    df <- df %>% 
      filter(!grepl("^Silvery*", get(TaxLvl)))
  } else if(TaxLvl == "Species"){
    df <- df %>% 
      filter(!grepl("* sp$", get(TaxLvl)))}
  #Both methods - all data
  df %>% 
    distinct_at(vars(TaxLvl)) %>% 
    nrow() %>% 
    paste0(TaxLvl, ": ", .) %>% 
    print()
  #By Method and Bioregion
  df %>% 
    filter(!is.na(N)) %>% 
    group_by(Bioregion, Method) %>%
    distinct_at(vars(TaxLvl)) %>% 
    summarise_at(vars(TaxLvl), funs(N = n())) %>% 
    print()
  #By Method
  df %>% 
    filter(!is.na(N)) %>% 
    group_by(Method) %>%
    distinct_at(vars(TaxLvl)) %>% 
    summarise_at(vars(TaxLvl), funs(N = n())) %>% 
    print()
  #By bioregion
  df %>% 
    filter(!is.na(N)) %>% 
    group_by(Bioregion) %>%
    distinct_at(vars(TaxLvl)) %>% 
    summarise_at(vars(TaxLvl), funs(N = n())) %>% 
    print()}

#Unique to each method
UniqueMethod <- function(df, TaxLvl){
  if(TaxLvl == "Family"){
    df <- df %>% 
      filter(!grepl("^Silvery*", get(TaxLvl)))
  } else if(TaxLvl == "Species"){
    df <- df %>% 
      filter(!grepl("* sp$", get(TaxLvl)))}
  #UVC
  print("UVC")
  df %>% 
  filter(Method == "UVC" & !is.na(N)) %>% 
  distinct_at(vars(TaxLvl)) %>% 
  anti_join(distinct_at(df %>% filter(Method == "BRUV" & !is.na(N)) %>% 
                          distinct_at(vars(TaxLvl)), vars(TaxLvl))) %>% 
    print()
  #UVC per bioregion
  print("UVC per bioregion")
  df %>% 
    filter(Method == "UVC" & !is.na(N)) %>% 
    group_by(Bioregion) %>% 
    distinct_at(vars(TaxLvl)) %>% 
    anti_join(distinct_at(df %>% filter(Method == "BRUV" & !is.na(N)) %>% group_by(Bioregion) %>% 
                            distinct_at(vars(TaxLvl)), vars(TaxLvl))) %>% 
    print()
  #BRUV
  print("BRUV")
  df %>% 
    filter(Method == "BRUV" & !is.na(N)) %>% 
    distinct_at(vars(TaxLvl)) %>% 
    anti_join(distinct_at(df %>% filter(Method == "UVC" & !is.na(N)) %>% 
                            distinct_at(vars(TaxLvl)), vars(TaxLvl))) %>% 
    print()
  #BRUV per bioregion
  print("BRUV per bioregion")
  df %>% 
    filter(Method == "BRUV" & !is.na(N)) %>% 
    group_by(Bioregion) %>% 
    distinct_at(vars(TaxLvl)) %>% 
    anti_join(distinct_at(df %>% filter(Method == "UVC" & !is.na(N)) %>% group_by(Bioregion) %>% 
                            distinct_at(vars(TaxLvl)), vars(TaxLvl))) %>% 
  print()}

BasicStats(CountsALL, "Family")
```

    ## Note: Using an external vector in selections is ambiguous.
    ## i Use `all_of(TaxLvl)` instead of `TaxLvl` to silence this message.
    ## i See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

    ## [1] "Family: 36"
    ## # A tibble: 4 x 3
    ## # Groups:   Bioregion [2]
    ##   Bioregion Method     N
    ##   <chr>     <chr>  <int>
    ## 1 CSE       BRUV      29
    ## 2 CSE       UVC       29
    ## 3 Western   BRUV      22
    ## 4 Western   UVC       19
    ## # A tibble: 2 x 2
    ##   Method     N
    ##   <chr>  <int>
    ## 1 BRUV      33
    ## 2 UVC       30
    ## # A tibble: 2 x 2
    ##   Bioregion     N
    ##   <chr>     <int>
    ## 1 CSE          34
    ## 2 Western      24

``` r
UniqueMethod(CountsALL, "Family")
```

    ## [1] "UVC"

    ## Joining, by = "Family"

    ## # A tibble: 3 x 1
    ##   Family     
    ##   <chr>      
    ## 1 Cirrhitidae
    ## 2 Apogonidae 
    ## 3 Sciaenidae 
    ## [1] "UVC per bioregion"

    ## Joining, by = c("Family", "Bioregion")

    ## # A tibble: 7 x 2
    ## # Groups:   Bioregion [2]
    ##   Family        Bioregion
    ##   <chr>         <chr>    
    ## 1 Cirrhitidae   CSE      
    ## 2 Labrisomidae  CSE      
    ## 3 Apogonidae    CSE      
    ## 4 Sciaenidae    CSE      
    ## 5 Kyphosidae    CSE      
    ## 6 Pomacanthidae Western  
    ## 7 Balistidae    Western  
    ## [1] "BRUV"

    ## Joining, by = "Family"

    ## # A tibble: 6 x 1
    ##   Family       
    ##   <chr>        
    ## 1 Scombridae   
    ## 2 Hemiramphidae
    ## 3 Chanidae     
    ## 4 Gobiidae     
    ## 5 Centropomidae
    ## 6 Mullidae     
    ## [1] "BRUV per bioregion"

    ## Joining, by = c("Family", "Bioregion")

    ## # A tibble: 10 x 2
    ## # Groups:   Bioregion [2]
    ##    Family        Bioregion
    ##    <chr>         <chr>    
    ##  1 Scombridae    CSE      
    ##  2 Hemiramphidae CSE      
    ##  3 Diodontidae   CSE      
    ##  4 Chanidae      CSE      
    ##  5 Gobiidae      CSE      
    ##  6 Centropomidae Western  
    ##  7 Carangidae    Western  
    ##  8 Scombridae    Western  
    ##  9 Mullidae      Western  
    ## 10 Labrisomidae  Western

#### Genus

``` r
BasicStats(CountsALL, "Genus")
```

    ## [1] "Genus: 67"
    ## # A tibble: 4 x 3
    ## # Groups:   Bioregion [2]
    ##   Bioregion Method     N
    ##   <chr>     <chr>  <int>
    ## 1 CSE       BRUV      48
    ## 2 CSE       UVC       50
    ## 3 Western   BRUV      35
    ## 4 Western   UVC       31
    ## # A tibble: 2 x 2
    ##   Method     N
    ##   <chr>  <int>
    ## 1 BRUV      54
    ## 2 UVC       55
    ## # A tibble: 2 x 2
    ##   Bioregion     N
    ##   <chr>     <int>
    ## 1 CSE          62
    ## 2 Western      45

#### Species

``` r
BasicStats(CountsALL, "Species")
```

    ## [1] "Species: 92"
    ## # A tibble: 4 x 3
    ## # Groups:   Bioregion [2]
    ##   Bioregion Method     N
    ##   <chr>     <chr>  <int>
    ## 1 CSE       BRUV      61
    ## 2 CSE       UVC       66
    ## 3 Western   BRUV      42
    ## 4 Western   UVC       41
    ## # A tibble: 2 x 2
    ##   Method     N
    ##   <chr>  <int>
    ## 1 BRUV      69
    ## 2 UVC       72
    ## # A tibble: 2 x 2
    ##   Bioregion     N
    ##   <chr>     <int>
    ## 1 CSE          85
    ## 2 Western      60

``` r
UniqueMethod(CountsALL, "Species")
```

    ## [1] "UVC"

    ## Joining, by = "Species"

    ## # A tibble: 23 x 1
    ##    Species                    
    ##    <fct>                      
    ##  1 Cirrhitus rivulatus        
    ##  2 Xenichthys agassizii       
    ##  3 Aetobatus ocellatus        
    ##  4 Malacoctenus tetranemus    
    ##  5 Arothron hispidus          
    ##  6 Calamus brachysomus        
    ##  7 Lutjanus jordani           
    ##  8 Apogon atradorsatus        
    ##  9 Odontoscion eurymesops     
    ## 10 Ophioblennius steindachneri
    ## # ... with 13 more rows
    ## [1] "UVC per bioregion"

    ## Joining, by = c("Species", "Bioregion")

    ## # A tibble: 42 x 2
    ## # Groups:   Bioregion [2]
    ##    Species                 Bioregion
    ##    <fct>                   <chr>    
    ##  1 Lutjanus aratus         CSE      
    ##  2 Cirrhitus rivulatus     CSE      
    ##  3 Xenichthys agassizii    CSE      
    ##  4 Aetobatus ocellatus     CSE      
    ##  5 Malacoctenus tetranemus CSE      
    ##  6 Arothron hispidus       CSE      
    ##  7 Calamus brachysomus     CSE      
    ##  8 Lutjanus jordani        CSE      
    ##  9 Apogon atradorsatus     CSE      
    ## 10 Odontoscion eurymesops  CSE      
    ## # ... with 32 more rows
    ## [1] "BRUV"

    ## Joining, by = "Species"

    ## # A tibble: 20 x 1
    ##    Species                   
    ##    <fct>                     
    ##  1 Aetobatus laticeps        
    ##  2 Eucinostomus currani      
    ##  3 Arothron meleagris        
    ##  4 Thalassoma grammaticum    
    ##  5 Scomberomorus sierra      
    ##  6 Haemulon steindachneri    
    ##  7 Orthopristis lethopristis 
    ##  8 Orthopristis cantharinus  
    ##  9 Selar crumenophthalmus    
    ## 10 Sphoeroides angusticeps   
    ## 11 Gymnothorax dovii         
    ## 12 Xenomugil thoburni        
    ## 13 Naucrates ductor          
    ## 14 Chanos chanos             
    ## 15 Muraena clepsydra         
    ## 16 Triaenodon obesus         
    ## 17 Centropomus viridis       
    ## 18 Kyphosus elegans          
    ## 19 Pseudupeneus grandisquamis
    ## 20 Labrisomus dendriticus    
    ## [1] "BRUV per bioregion"

    ## Joining, by = c("Species", "Bioregion")

    ## # A tibble: 38 x 2
    ## # Groups:   Bioregion [2]
    ##    Species                   Bioregion
    ##    <fct>                     <chr>    
    ##  1 Aetobatus laticeps        CSE      
    ##  2 Eucinostomus currani      CSE      
    ##  3 Arothron meleagris        CSE      
    ##  4 Thalassoma grammaticum    CSE      
    ##  5 Scomberomorus sierra      CSE      
    ##  6 Haemulon steindachneri    CSE      
    ##  7 Orthopristis lethopristis CSE      
    ##  8 Orthopristis cantharinus  CSE      
    ##  9 Serranus psittacinus      CSE      
    ## 10 Selar crumenophthalmus    CSE      
    ## # ... with 28 more rows

#### Individuals

``` r
#Individual Fish
CountsALL %>% summarise(Tot = sum(N, na.rm = T))
```

    ## # A tibble: 1 x 1
    ##     Tot
    ##   <int>
    ## 1 35029

``` r
#Individuals per method and per bioregion
CountsALL %>% group_by(Method, Bioregion) %>% 
  summarise(Tot = sum(N, na.rm = T))
```

    ## # A tibble: 4 x 3
    ## # Groups:   Method [2]
    ##   Method Bioregion   Tot
    ##   <chr>  <chr>     <int>
    ## 1 BRUV   CSE        9307
    ## 2 BRUV   Western    2420
    ## 3 UVC    CSE       19786
    ## 4 UVC    Western    3516

#### Number of fish classified as Silvery Fish

``` r
#Fish that could not be identified to species level (Silvery Fish)
CountsALL %>% filter(grepl("*silvery*", Species)) %>% summarise(Tot = sum(N, na.rm = T)) %>%
  pull() #pull extracts the number from summary
```

    ## [1] 2343

### Biodiversity measures

``` r
#Constructing matrix for multivariate analysis
FM_ALL <- CountsALL %>% 
  ungroup() %>% 
  select(SiteName, Species, N) %>% 
  pivot_wider(names_from = Species, values_from = N) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "SiteName") %>% as.matrix() %>% replace_na(., 0) %>% 
  .[order(rownames(.)),]
head(FM_ALL[,1:3])
```

    ##          Abudefduf concolor Abudefduf troschelii Aetobatus laticeps
    ## ALB_1_B                   0                    0                  0
    ## ALB_2_B                   0                    0                  1
    ## ALB_3_B                   0                    2                  0
    ## ALB1_1_U                  0                    0                  0
    ## ALB1_2_U                  0                    0                  0
    ## ALB1_3_U                  0                    4                  0

``` r
#Creating database of biodiversity measures
BioAll <- tibble(SiteName = row.names(FM_ALL),
                 SpRich = specnumber(FM_ALL),
                 N = rowSums(FM_ALL),
                 Shannon = round(diversity(FM_ALL), 3),
                 Pielou = round(Shannon/log(N), 3)) %>% 
  left_join(select(CountsALL, Bioregion, SiteName, Method) %>% unique(), ., by = "SiteName")
```

    ## Warning: Column `SiteName` joining factor and character vector, coercing into
    ## character vector

#### Plotting mean abundance

``` r
#Number of individuals
p1 <- ggplot(BioAll, aes(Method, log(N+1)))+geom_boxplot()+
  labs(y = "Abundance \n (log(x+1))")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p1
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/mean%20abundance-1.png)<!-- -->

``` r
p2 <- ggplot(BioAll, aes(Bioregion, log(N+1)))+geom_boxplot()+facet_grid(.~Method)+
  labs(y = "Abundance \n (log(x+1))")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p2
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/mean%20abundance-2.png)<!-- -->

``` r
p3 <- ggplot(BioAll, aes(Bioregion, log(N+1)))+geom_boxplot()+
  labs(y = "")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p3
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/mean%20abundance-3.png)<!-- -->

``` r
MeanAbund <- ggarrange(ggarrange(p1, p3, ncol = 2, labels = c("A", "B"),
                                 font.label = list(size = 12)),
            p2, nrow = 2, labels = c("", "C"), font.label = list(size = 12))
MeanAbund
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/mean%20abundance-4.png)<!-- -->

``` r
rm(p1, p2, p3)
ggsave("Figures/LogAbundance_AllPts.tiff", MeanAbund, device = "tiff", dpi = 500)
```

    ## Saving 7 x 5 in image

``` r
#Testing for differences between methods and bioregions - Univariate PERMANOVA
adonis(log(N+1) ~ Bioregion*Method, data = BioAll, method = "euc", permutations = 9999)
```

    ## 
    ## Call:
    ## adonis(formula = log(N + 1) ~ Bioregion * Method, data = BioAll,      permutations = 9999, method = "euc") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
    ## Bioregion          1    12.262 12.2616  9.7262 0.04184 0.0024 **
    ## Method             1     3.316  3.3161  2.6304 0.01131 0.1093   
    ## Bioregion:Method   1     2.685  2.6848  2.1296 0.00916 0.1488   
    ## Residuals        218   274.830  1.2607         0.93769          
    ## Total            221   293.092                 1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Plotting species richness

``` r
#Species richness
p1 <- ggplot(BioAll, aes(Method, SpRich))+geom_boxplot()+
  labs(y = "Species richness")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p1
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/species%20richness%20plot-1.png)<!-- -->

``` r
p2 <- ggplot(BioAll, aes(Bioregion, SpRich))+geom_boxplot()+
  facet_grid(~Method)+labs(y = "Species richness")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p2
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/species%20richness%20plot-2.png)<!-- -->

``` r
p3 <- ggplot(BioAll, aes(Bioregion, SpRich))+geom_boxplot()+
  labs(y = "")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p3
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/species%20richness%20plot-3.png)<!-- -->

``` r
SpRich <- ggarrange(ggarrange(p1, p3, ncol = 2, labels = c("A", "B"), 
                              font.label = list(size = 12), vjust = 1),
                    p2, nrow = 2, labels = c("", "C"), 
                    font.label = list(size = 12), vjust = 1)
SpRich
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/species%20richness%20plot-4.png)<!-- -->

``` r
rm(p1, p2, p3)
ggsave("Figures/SpeciesRichness_AllPts.tiff", SpRich, device = "tiff",
       dpi = 500)
```

    ## Saving 7 x 5 in image

``` r
#Testing for differences between methods and bioregions - Univariate PERMANOVA
adonis(SpRich ~ Bioregion*Method, data = BioAll, method = "euc", permutations = 9999)
```

    ## 
    ## Call:
    ## adonis(formula = SpRich ~ Bioregion * Method, data = BioAll,      permutations = 9999, method = "euc") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Bioregion          1      73.7   73.73   5.716 0.01763 0.0179 *  
    ## Method             1    1243.0 1243.04  96.380 0.29726 0.0001 ***
    ## Bioregion:Method   1      53.3   53.28   4.131 0.01274 0.0445 *  
    ## Residuals        218    2811.6   12.90         0.67237           
    ## Total            221    4181.7                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Plotting Shannon diversity index

``` r
#Species richness
p1 <- ggplot(BioAll, aes(Method, Shannon))+geom_boxplot()+
  labs(y = "Shannon index")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p1
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Shannon%20diversity-1.png)<!-- -->

``` r
p2 <- ggplot(BioAll, aes(Bioregion, Shannon))+geom_boxplot()+
  facet_grid(~Method)+labs(y = "Shannon index")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p2
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Shannon%20diversity-2.png)<!-- -->

``` r
p3 <- ggplot(BioAll, aes(Bioregion, Shannon))+geom_boxplot()+
  labs(y = "")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p3
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Shannon%20diversity-3.png)<!-- -->

``` r
ShannonInd <- ggarrange(ggarrange(p1, p3, ncol = 2, labels = c("A", "B"), 
                              font.label = list(size = 12), vjust = 1.25),
                    p2, nrow = 2, labels = c("", "C"), 
                    font.label = list(size = 12), vjust = 1)
ShannonInd
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Shannon%20diversity-4.png)<!-- -->

``` r
rm(p1, p2, p3)
ggsave("Figures/ShannonIndex_AllPts.tiff", ShannonInd, device = "tiff", dpi = 500)
```

    ## Saving 7 x 5 in image

``` r
#Testing for differences between methods and bioregions - Univariate PERMANOVA
adonis(Shannon ~ Bioregion*Method, data = BioAll, method = "euc", permutations = 9999)
```

    ## 
    ## Call:
    ## adonis(formula = Shannon ~ Bioregion * Method, data = BioAll,      permutations = 9999, method = "euc") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Bioregion          1     0.523  0.5229  2.1256 0.00896 0.1479    
    ## Method             1     4.131  4.1312 16.7940 0.07076 0.0004 ***
    ## Bioregion:Method   1     0.100  0.1003  0.4078 0.00172 0.5160    
    ## Residuals        218    53.626  0.2460         0.91856           
    ## Total            221    58.380                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Plotting Pielou evenness index

``` r
#Species evenness
p1 <- ggplot(BioAll, aes(Method, Pielou))+geom_boxplot()+
  labs(y = "Pielou's evenness")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p1
```

    ## Warning: Removed 1 rows containing non-finite values (stat_boxplot).

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Pielou%20evenness-1.png)<!-- -->

``` r
p2 <- ggplot(BioAll, aes(Bioregion, Pielou))+geom_boxplot()+
  facet_grid(~Method)+labs(y = "Pielou's evenness")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p2
```

    ## Warning: Removed 1 rows containing non-finite values (stat_boxplot).

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Pielou%20evenness-2.png)<!-- -->

``` r
p3 <- ggplot(BioAll, aes(Bioregion, Pielou))+geom_boxplot()+
  labs(y = "")+theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 11), 
        axis.text = element_text(family = "sans", size = 12))
p3
```

    ## Warning: Removed 1 rows containing non-finite values (stat_boxplot).

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Pielou%20evenness-3.png)<!-- -->

``` r
Pielou <- ggarrange(ggarrange(p1, p3, ncol = 2, labels = c("A", "B"), 
                              font.label = list(size = 12), vjust = 1,
                              hjust = -0.3),
                    p2, nrow = 2, labels = c("", "C"), 
                    font.label = list(size = 12), vjust = 1)
```

    ## Warning: Removed 1 rows containing non-finite values (stat_boxplot).
    
    ## Warning: Removed 1 rows containing non-finite values (stat_boxplot).
    
    ## Warning: Removed 1 rows containing non-finite values (stat_boxplot).

``` r
Pielou
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Pielou%20evenness-4.png)<!-- -->

``` r
rm(p1, p2, p3)
ggsave("Figures/Pielou_AllPts.tiff", Pielou, device = "tiff", dpi = 500)
```

    ## Saving 7 x 5 in image

``` r
#Testing for differences between methods and bioregions - Univariate PERMANOVA
adonis(Pielou ~ Bioregion*Method, data = BioAll %>% drop_na(), method = "euc", permutations = 9999)
```

    ## 
    ## Call:
    ## adonis(formula = Pielou ~ Bioregion * Method, data = BioAll %>%      drop_na(), permutations = 9999, method = "euc") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## Bioregion          1    0.0102 0.010151  0.4601 0.00197 0.4924    
    ## Method             1    0.2635 0.263531 11.9450 0.05112 0.0009 ***
    ## Bioregion:Method   1    0.0940 0.094043  4.2627 0.01824 0.0389 *  
    ## Residuals        217    4.7875 0.022062         0.92867           
    ## Total            220    5.1552                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Plotting fish community composition - Families

``` r
#By method
Counts <- CountsALL %>%
  drop_na(N) %>% 
  select(Method, Bioregion, Family, N) %>% 
  group_by(Method) %>% 
  mutate(TotalAb = sum(N, na.rm = T), Family = as.character(Family)) %>% 
  group_by(Method, Family) %>% 
  summarise(Sum = sum(N, na.rm = T), 
            Mean = mean(N, na.rm = T),
            SE = plotrix::std.error(N, na.rm = T),
            TotalAb = mean(TotalAb, na.rm = T), 
            Prop = Sum/TotalAb*100) %>% 
  mutate(Family = replace(Family, Prop < 5, "Other Families")) %>% 
  group_by(Method, Family) %>% 
  summarise(Prop = sum(Prop, na.rm = T), N = n()) %>% 
  ungroup() %>% 
  mutate(Method = factor(Method),
         Family = factor(Family, ordered = T),
         Ord = case_when(Family == "Other Families" ~ 1,
                         TRUE ~ Prop)) %>%
  mutate(Family = fct_reorder(Family, Ord, .desc = F))
```

``` r
#No longer needed as "Other families" annotations have been removed
# #Labels
# annot_text <- Counts %>% select(-Prop) %>% 
#   filter(grepl("Other Families", Family)) %>% 
#   unite(Label, Family, N, sep = ": ")

# Graphs - Bar plots
#Factor(1) included under x to be able to use Method under facets
p1 <- Counts %>% ggplot(aes(x = factor(1), y = Prop, 
                            #Reordering Families based on order and method prior to using it as
                            #categories under fill
                            fill = tidytext::reorder_within(Family, Ord, Method)))+
  #Color "white" refers to the colour of the line separating each block within the bar
  geom_bar(position = "stack", color = "black", stat = "identity")+
  facet_grid(facets = Method ~ .)+
  #Using a color scale that is colour blind friendly
  scale_fill_manual(values = c("#DDDDDD", "#DDDDDD", "#332288", "#88CCEE", "#999933", "#44AA99",
                               "#117733", "#DDCC77", "#999933", "#cc6677", "#882255", "#cc6677"))+
  #Flipping bars so they are horizontal
  coord_flip()+
  #Changing labels in axes
  labs(x = "Methods", y = "")+
  #Changing theme to meet magazine requirements
  theme_bw()+theme(axis.title = element_text(family = "sans", size = 12),
                   axis.text.y = element_blank(),
                   axis.text = element_text(family = "sans", size = 12),
                   #Removing legend
                   legend.position = "none",
                   #Removing y axis ticks
                   axis.ticks.y = element_blank())
#Removing in figure annotations including the actual number of other families sampled as per
#reviewer's comments.
# +
#   annotate(geom = "text", x = 2.5, y = 75, 
#            label = annot_text %>% filter(grepl("BRUV", Method)) %>% 
#              select(Label) %>% paste(.))+
#   annotate(geom = "text", x = 1.5, y = 75, 
#            label = annot_text %>% filter(grepl("UVC", Method)) %>% 
#              select(Label) %>% paste(.))
p1
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/graph%20fish%20community%20by%20method-1.png)<!-- -->

``` r
rm(Counts)
```

``` r
##By method and bioregion
CountsBioM <- CountsALL %>% 
  drop_na(N) %>% 
  select(Method, Bioregion, Family, N) %>% 
  group_by(Method, Bioregion) %>% 
  mutate(TotalAb = sum(N, na.rm = T), Family = as.character(Family)) %>% 
  group_by(Method, Bioregion, Family) %>% 
  summarise(Sum = sum(N, na.rm = T), 
            Mean = mean(N, na.rm = T),
            SE = plotrix::std.error(N, na.rm = T),
            TotalAb = mean(TotalAb, na.rm = T), 
            Prop = Sum/TotalAb*100) %>% 
  mutate(Family = replace(Family, Prop < 5, "Other Families")) %>% 
  group_by(Method, Bioregion, Family) %>% 
  summarise(Prop = sum(Prop, na.rm = T), N = n()) %>% 
  ungroup() %>% 
  mutate(Family = factor(Family, ordered = T), 
         Ord = case_when(Family == "Other Families" ~ 1,
                         TRUE ~ Prop)) %>% 
  mutate(Family = fct_reorder(Family, Ord, .desc = F))
```

``` r
#Graph - Bar plots
p2 <- CountsBioM %>% 
  #Creating a new column to reorder families using both methods and bioregions
  mutate(Levels = paste0(Method, Bioregion)) %>% 
  ggplot(aes(x = "", y = Prop, fill = tidytext::reorder_within(Family, Ord, Levels)))+
  geom_bar(position = "stack", color = "black", stat = "identity")+
  facet_grid(Method~Bioregion)+
  scale_fill_manual(values = c("#dddddd", "#dddddd", "#dddddd", "#dddddd", "#332288", "#aa4499",
                               "#332288", "#cc6677", "#882255", "#999933", "#44aa99", "#88ccee",
                               "#44aa99", "#999933", "#ddcc77", "#999933", "#117733", "#ddcc77",
                               "#882255", "#cc6677", "#cc6677", "#999933", "#882255"))+
  coord_flip()+
  labs(x = "", y = "")+
  theme_bw()+theme(axis.text.x = element_text(family = "sans", size = 12), 
                   axis.ticks.y = element_blank(),
                   legend.position = "none")
p2
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/graph%20fish%20community%20by%20method%20&%20bioregion-1.png)<!-- -->

``` r
#Extracting legend from dummy graph including all families shared across graphs to be used in the
#final Family Composition figure
#Creating graph with correct colour palette and giving desired format to the legend
g <- CountsBioM %>% ggplot(aes(x = "", y = Prop, fill = Family))+
  geom_bar(position = "stack", stat = "identity", color = "black")+
  scale_fill_manual(values = c("#dddddd", "#332288", "#aa4499", "#88ccee", "#44aa99", "#999933",
                               "#117733", "#ddcc77", "#882255", "#cc6677"))+
  #Legend to be divided in two rows with title on top of categories
  guides(fill = guide_legend(nrow = 2, title.position = "top"))+
  #Changing font and size of legend
  theme(legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12))

#Extracting legend to be used in final composite figure
p2leg <- get_legend(g)

rm(CountsBioM, g)
```

``` r
CountsBio <- CountsALL %>% 
  drop_na(N) %>% 
  select(Bioregion, Family, N) %>% 
  group_by(Bioregion) %>% 
  mutate(TotalAb = sum(N, na.rm = T), Family = as.character(Family)) %>% 
  group_by(Bioregion, Family) %>% 
  summarise(Sum = sum(N, na.rm = T), 
            Mean = mean(N, na.rm = T),
            SE = plotrix::std.error(N, na.rm = T),
            TotalAb = mean(TotalAb, na.rm = T), 
            Prop = Sum/TotalAb*100) %>% 
  mutate(Family = replace(Family, Prop < 5, "Other Families")) %>% 
  group_by(Bioregion, Family) %>% 
  summarise(Prop = sum(Prop, na.rm = T), N = n()) %>% 
  ungroup() %>% 
  mutate(Family = factor(Family, ordered = T), 
         Ord = case_when(Family == "Other Families" ~ 1,
                         TRUE ~ Prop)) %>% 
  mutate(Family = fct_reorder(Family, Ord, .desc = F))

#Graph - Bar plots
p3 <- CountsBio %>% ggplot(aes(x = factor(1), y = Prop, 
                               fill = tidytext::reorder_within(Family, Ord, Bioregion)))+
  geom_bar(position = "stack", color = "black", stat = "identity")+
  facet_grid(facets = Bioregion~.)+
  scale_fill_manual(values = c("#dddddd", "#dddddd", "#ddcc77", "#332288", "#44aa99", "#ddcc77",
                               "#44aa99", "#999933", "#117733", "#882255", "#999933", "#cc6677",
                               "#882255"))+
  coord_flip()+
  labs(x = "Bioregion", y = "% total abundance")+
  theme_bw()+theme(axis.title = element_text(family = "sans", size = 12),
                   axis.text.x = element_text(family = "sans", size = 12),
                   axis.text.y = element_blank(),
                   legend.position = "none",
                   axis.ticks.y = element_blank())
p3
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/fish%20community%20by%20bioregion-1.png)<!-- -->

``` r
#Merging the two plots into one panel
CommComp <- ggarrange(plot_grid(p2leg, plot_grid(p1, p2+theme(legend.position = "none"), 
                                     align = "hv", axis = "lb", labels = c("A", "B"), 
                                     hjust = -1, rel_widths = c(0.8, 1)), 
                    ncol = 1, rel_heights = c(0.4, 1, 0.7),
          p3, labels = c("", "", "C"), nrow = 3, hjust = -1, vjust = 1))
CommComp
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/fish%20community%20by%20bioregion-2.png)<!-- -->

``` r
ggsave("Figures/CommunityFamCompAllPts.tiff", CommComp, device = "tiff", dpi = 500, 
       height = 18.23, width = 20, units = "cm")
rm(p1, p2, p3, CountsBio, p2leg)
```

### Comparing trophic levels between methods and across regions

``` r
#Creating dataframe with summary of trophic data - Blank category appears for individuals not identified to species level
Troph <- CountsALL %>% 
  group_by(Method, SiteName) %>% 
  mutate(TotalAb = sum(N, na.rm = T)) %>% 
  group_by(Method, SiteName, TrophicCat) %>%
  mutate(TotTrop = sum(N, na.rm = T), 
         PropTrop = TotTrop / TotalAb) %>% 
  group_by(Method, Bioregion, SiteName, TrophicCat) %>% 
  summarise_at(vars(PropTrop), funs(mean(., na.rm = T))) 

#Summary and graph at a method level
p1 <- Troph %>% group_by(Method, TrophicCat) %>% 
  summarise(MeanTrop = mean(PropTrop, na.rm = T), 
            SE_Trop = plotrix::std.error(PropTrop, na.rm = T)) %>% 
  filter(!TrophicCat == "") %>% 
  ggplot(aes(TrophicCat, MeanTrop, fill = Method))+
  geom_bar(stat = "identity", position = position_dodge(), color = "black")+
  scale_fill_manual(values = c("#E1E1E1", "white"))+
  geom_errorbar(aes(ymin = MeanTrop, ymax = MeanTrop+SE_Trop), width = 0.2,
                position = position_dodge(.9))+
  theme_bw()+scale_y_continuous(expand = expansion(mult = c(0,0.2)))+
  labs(x = "", y = "")+
  theme(axis.title = element_text(family = "sans", size = 12), 
                   axis.text = element_text(family = "sans", size = 12),
        legend.position = "top")+
  annotate(geom = "text", x = 1:4, y = c(0.1, 1, 0.1, 0.375), label = "*", size = 10)
p1
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Trophic%20levels-1.png)<!-- -->

``` r
#Summary and graph at a method and bioregion level
p2 <- Troph %>% group_by(Method, Bioregion, TrophicCat) %>%
  summarise(MeanTrop = mean(PropTrop, na.rm = T),
            SE_Trop = plotrix::std.error(PropTrop, na.rm = T)) %>%
  filter(!TrophicCat == "") %>%
  ggplot(aes(TrophicCat, MeanTrop, fill = Bioregion))+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  scale_fill_manual(values = c("#808080", "#000000"))+
  geom_errorbar(aes(ymin = MeanTrop, ymax = MeanTrop+SE_Trop), width = 0.2,
                position = position_dodge(.9))+
  facet_wrap(~Method)+theme_bw()+
  labs(x = "", y = "")+
  theme(axis.title = element_text(family = "sans", size = 12),
                   axis.text = element_text(family = "sans", size = 12),
        legend.position = "top")
p2
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Trophic%20levels-2.png)<!-- -->

``` r
#Summary and graph across bioregions
p3 <- Troph %>% group_by(Bioregion, TrophicCat) %>%
  summarise(MeanTrop = mean(PropTrop, na.rm = T),
            SE_Trop = plotrix::std.error(PropTrop, na.rm = T)) %>%
  filter(!TrophicCat == "") %>%
  ggplot(aes(TrophicCat, MeanTrop, fill = Bioregion))+
  geom_bar(stat = "identity", position = position_dodge(), color = "black")+
  scale_fill_manual(values = c("#808080", "#000000"))+
  geom_errorbar(aes(ymin = MeanTrop, ymax = MeanTrop+SE_Trop), width = 0.2,
                position = position_dodge(.9))+
  theme_bw()+scale_y_continuous(expand = expansion(mult = c(0,0.2)))+
  labs(x = "", y = "")+
  theme(axis.title = element_text(family = "sans", size = 12),
                   axis.text = element_text(family = "sans", size = 12),
        legend.position = "top")+
  annotate(geom = "text", x = 3, y = 0.125, label = "*", size = 10)
p3
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Trophic%20levels-3.png)<!-- -->

``` r
TrophCat <- ggarrange(p1, ggarrange(p2, p3, labels = c("B", "C"), nrow = 2, common.legend = T),
                      labels = c("A", "", ""), nrow = 2, heights = c(1, 2))
TrophCat <- annotate_figure(TrophCat, left = text_grob("Mean relative abundance", rot = 90, hjust = 0.5, size = 12))
TrophCat
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Trophic%20levels-4.png)<!-- -->

``` r
rm(p1, p2, p3)
ggsave("Figures/TrophicProportions.tiff", TrophCat, device = "tiff", dpi = 500)
```

    ## Saving 7 x 5 in image

``` r
#Testing for differences between methods and bioregions - Univariate PERMANOVA
#Carnivores
adonis(PropTrop ~ Bioregion*Method, data = Troph %>% filter(TrophicCat == "Carnivore"), 
       method = "euc", permutations = 9999, na.rm = T)
```

    ## 
    ## Call:
    ## adonis(formula = PropTrop ~ Bioregion * Method, data = Troph %>%      filter(TrophicCat == "Carnivore"), permutations = 9999, method = "euc",      na.rm = T) 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Bioregion          1    0.1638  0.1638   3.043 0.00956 0.0813 .  
    ## Method             1    5.2257  5.2257  97.094 0.30501 0.0001 ***
    ## Bioregion:Method   1    0.0101  0.0101   0.188 0.00059 0.6660    
    ## Residuals        218   11.7330  0.0538         0.68483           
    ## Total            221   17.1325                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Herbivores
adonis(PropTrop ~ Bioregion*Method, data = Troph %>% filter(TrophicCat == "Herbivore"), 
       method = "euc", permutations = 9999, na.rm = T)
```

    ## 
    ## Call:
    ## adonis(formula = PropTrop ~ Bioregion * Method, data = Troph %>%      filter(TrophicCat == "Herbivore"), permutations = 9999, method = "euc",      na.rm = T) 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
    ## Bioregion          1   0.06959 0.069594  9.7854 0.04199 0.0018 **
    ## Method             1   0.03433 0.034332  4.8274 0.02072 0.0262 * 
    ## Bioregion:Method   1   0.00295 0.002947  0.4144 0.00178 0.5228   
    ## Residuals        218   1.55041 0.007112         0.93551          
    ## Total            221   1.65728                  1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Apex
adonis(PropTrop ~ Bioregion*Method, data = Troph %>% filter(TrophicCat == "Apex"), 
       method = "euc", permutations = 9999, na.rm = T)
```

    ## 
    ## Call:
    ## adonis(formula = PropTrop ~ Bioregion * Method, data = Troph %>%      filter(TrophicCat == "Apex"), permutations = 9999, method = "euc",      na.rm = T) 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## Bioregion          1   0.00680 0.006800  1.4177 0.00617 0.2372    
    ## Method             1   0.03822 0.038225  7.9694 0.03469 0.0001 ***
    ## Bioregion:Method   1   0.01137 0.011365  2.3695 0.01031 0.1289    
    ## Residuals        218   1.04562 0.004796         0.94883           
    ## Total            221   1.10201                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Omnivores
adonis(PropTrop ~ Bioregion*Method, data = Troph %>% filter(TrophicCat == "Omnivore"), 
       method = "euc", permutations = 9999, na.rm = T)
```

    ## 
    ## Call:
    ## adonis(formula = PropTrop ~ Bioregion * Method, data = Troph %>%      filter(TrophicCat == "Omnivore"), permutations = 9999, method = "euc",      na.rm = T) 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Bioregion          1    0.0689 0.06886  1.4807 0.00598 0.2256    
    ## Method             1    1.2431 1.24306 26.7315 0.10792 0.0001 ***
    ## Bioregion:Method   1    0.0686 0.06864  1.4761 0.00596 0.2241    
    ## Residuals        218   10.1374 0.04650         0.88014           
    ## Total            221   11.5180                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rm(Troph)
```

### Comparing fish assemblages between methods and bioregions

``` r
#Fish count data
#Comparing species composition between methods
FishComp <- FM_ALL %>% as.data.frame() %>%
  mutate(SiteName = rownames(FM_ALL)) %>% 
  left_join(select(BioAll, SiteName, Bioregion, Method), by = "SiteName") %>% 
  group_by(Method) %>% 
  summarise_if(is.numeric, mean) %>% 
  pivot_longer(-c(Method), names_to = "Species", values_to = "MeanProp") %>% 
  pivot_wider(names_from = Method, values_from = MeanProp) %>% 
  left_join(select(FishData, ValidName, TrophicCat) %>% distinct(), by = c("Species" = "ValidName")) %>%
  mutate(TrophicCat = as.character(TrophicCat)) %>% 
  mutate(TrophicCat = replace(TrophicCat, TrophicCat == "", "Unallocated")) %>% 
  ggplot(aes(x = BRUV, y = UVC))+geom_point(aes(shape = TrophicCat), size = 2)+
  geom_text_repel(aes(label = ifelse(BRUV > 10 | UVC > 10, Species, "")),
                  nudge_x = 1.25,
                  segment.size = 0.5,
                  size = 4)+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.2)))+
  geom_abline(slope = 1, linetype = "dashed")+
  theme_bw()+
  labs(x = "", y = "")+
  theme(axis.title = element_text(family = "sans", size = 12), 
                   axis.text = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12, face = "bold"),
        legend.position = "top")+
  scale_shape_discrete(name = "Trophic Category")+
  guides(shape = guide_legend(title.position = "top", title.hjust = 0.5))
FishComp
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/fish%20assemblages-1.png)<!-- -->

``` r
#Comparing species composition across bioregions - Only showing names of herbivores responsible for at least 5% of
#total composition
FishCompBio <- FM_ALL %>% as.data.frame() %>%
  mutate(SiteName = rownames(FM_ALL)) %>% 
  left_join(select(BioAll, SiteName, Bioregion, Method), by = "SiteName") %>% 
  group_by(Method, Bioregion) %>% 
  summarise_if(is.numeric, mean) %>% 
  pivot_longer(-c(Method, Bioregion), names_to = "Species", values_to = "MeanProp") %>% 
  pivot_wider(names_from = Method, values_from = MeanProp) %>% 
  left_join(select(FishData, ValidName, TrophicCat) %>% distinct(), by = c("Species" = "ValidName")) %>%
  mutate(TrophicCat = as.character(TrophicCat)) %>% 
  mutate(TrophicCat = replace(TrophicCat, TrophicCat == "", "Unallocated")) %>% 
  ggplot(aes(x = BRUV, y = UVC))+geom_point(aes(shape = TrophicCat), size = 2)+
  geom_text_repel(aes(label = ifelse(TrophicCat == "Herbivore" & (BRUV > 2.5 | UVC > 2.5), Species, "")),
                  nudge_x = 1.25,
                  segment.size = 0.5,
                  size = 4)+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.2)))+
  geom_abline(slope = 1, linetype = "dashed")+
  theme_bw()+facet_grid(Bioregion~.)+
  labs(x = "", y = "")+
  theme(axis.title = element_text(family = "sans", size = 12), 
                   axis.text = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12, face = "bold"),
        legend.position = "top")+
  scale_shape_discrete(name = "Trophic Category")+
  guides(shape = guide_legend(title.position = "top", title.hjust = 0.5))
FishCompBio
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/fish%20assemblages-2.png)<!-- -->

``` r
#Putting two graphs together
TrophicComp <- ggarrange(FishComp, FishCompBio, ncol = 2, common.legend = T, widths = c(1, 0.7))
TrophicComp <- annotate_figure(TrophicComp, left = text_grob("% UVC Abundance", vjust = 1.25,
                                                          rot = 90, hjust = 0.5, size = 12))
TrophicComp <- annotate_figure(TrophicComp, bottom = text_grob("% BRUV Abundance", hjust = 0.5, 
                                                               vjust = -0.5, size = 12))
TrophicComp
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/fish%20assemblages-3.png)<!-- -->

``` r
ggsave("Figures/MethodDetectionDifferences.tiff", TrophicComp, device = "tiff", dpi = 500)
```

    ## Saving 7 x 5 in image

``` r
rm(FishComp, FishCompBio)

#Simplfying the fish counts database: Calculation mean abundance values per site sampled
FM_Mean <- FM_ALL %>% as.data.frame() %>%
  mutate(SiteName = rownames(FM_ALL), 
         SiteCode = paste0(str_extract(SiteName, "[A-Z]{3,4}"), str_extract(SiteName, "[_][A-Z]{1}"))) %>%
  mutate(SiteCode = case_when(grepl("*_U$" , SiteCode) ~ SiteCode,
                              grepl("*_B$" , SiteCode) ~ SiteCode,
                              TRUE ~ paste0(SiteCode, "_B"))) %>% 
  group_by(SiteCode) %>% summarise_if(is.numeric, mean) %>% 
  remove_rownames() %>% column_to_rownames(var = "SiteCode") %>% 
  as.matrix()

#Creating new matrix with adjusted abundance values - Values as percentages
Norm_FM <- prop.table(FM_Mean, 1)*100

#Fourth root transformation applied to give more weight to rarer species
Trans_FM <- sqrt(sqrt(Norm_FM))

#Calculating Bray Curtis dissimilarity
BC_FM <- vegdist(Trans_FM)

#PCO
pco_FM <- wcmdscale(BC_FM, eig = T)
plot(pco_FM)
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/fish%20assemblages-4.png)<!-- -->

### Environmental variables

``` r
#Environmental variables
#Pooling environmental data together - Getting mean values per site sampled for further analysis
MeanEnviro <- plyr::rbind.fill(BRUVEnviro %>% 
                                group_by(SiteCode) %>% 
                                summarise_if(is.numeric, mean) %>% 
                                select(-c(Depth, Replicate)) %>%
                                left_join(select(BRUVEnviro, SiteCode, Bioregion) %>% 
                                            distinct(), by = "SiteCode") %>% 
                                mutate(SiteCode = paste0(SiteCode, "_B"), Method = "BRUV"), 
                              UVCEnviro %>% 
                                mutate(SiteCode = str_extract(SiteCode, "[A-Z]{3,4}")) %>% 
                                group_by(SiteCode) %>% 
                                summarise_if(is.numeric, mean) %>%  
                                select(-c(Depth)) %>% 
                                left_join(., select(UVCEnviro, SiteCode, Bioregion) %>% 
                                            mutate(SiteCode = str_extract(SiteCode, "[A-Z]{3,4}")) %>% 
                                            distinct(), by = "SiteCode") %>% 
                                mutate(SiteCode = paste0(SiteCode, "_U"), Method = "UVC")) %>%
  arrange(SiteCode)


AllEnviro <- plyr::rbind.fill(BRUVEnviro %>% unite(SiteCode, SiteCode, Replicate) %>% 
                                mutate(SiteCode = paste0(SiteCode, "_B"), Method = "BRUV") %>% 
                                select(-c(SiteName, Time.of.day, Facing.mangroves., Tide, 
                                          Substrate.type, Island, Depth, Habitat)) %>%
                                select(-starts_with("Yr")),
                              UVCEnviro %>% 
                                mutate(SiteCode = paste0(SiteCode, "_U"), Method = "UVC") %>% 
                                select(-c(Depth, SiteName, Habitat, Island)) %>% 
                                select(-starts_with("Yr"))) %>% 
  arrange(SiteCode) %>% write.csv("Outputs/RelativeFishAbundanceAllPts.csv")


#Visualising enviromental variables - Draftsman plots. To check if transformation is needed
MeanEnviro %>% select_if(is.numeric) %>% 
  mutate_at(vars(matches("area")), log) %>% GGally::ggpairs()
```

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning in (function (data, mapping, alignPercent = 0.6, method = "pearson", :
    ## Removed 22 rows containing missing values

    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).
    
    ## Warning: Removed 22 rows containing missing values (geom_point).

    ## Warning: Removed 22 rows containing non-finite values (stat_density).

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Enviro%20variables-1.png)<!-- -->

``` r
#The following variables variables were removed from further analysis to avoid multicollinearity:
#All variables removed had correlations >= 0.90
#Yearly Temperature and yearly chlorophyll as they showed high correlation with dry and wet season measurements
library(corrplot)
```

    ## corrplot 0.84 loaded

``` r
Env_Cor <- MeanEnviro %>% select_if(is.numeric) %>% select(-starts_with("Yr")) %>% 
  mutate_at(vars(matches("area")), log) %>% cor(use = "pairwise.complete.obs")
Env_PCor <- MeanEnviro %>% select_if(is.numeric) %>% select(-starts_with("Yr")) %>% 
  mutate_at(vars(matches("area")), log) %>% cor.mtest(use = "pairwise.complete.obs")
corrplot.mixed(Env_Cor, tl.pos = "lt", p.mat = Env_PCor$p, sig.level = 0.05)
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/Enviro%20variables-2.png)<!-- -->

``` r
#Normalising environmental data
Norm_Env <- MeanEnviro %>% select(-starts_with("Yr")) %>% 
  mutate_at(vars(matches("area")), log) %>% 
  mutate_if(is.numeric, list(~scale(.) %>% as.vector)) %>% 
  remove_rownames() %>% column_to_rownames(var = "SiteCode") %>% 
  select_if(is.numeric) %>% as.matrix() 

#Calculating Euclidean dissimilarity matrix
Euc_Env <- vegdist(Norm_Env, "euc", na.rm = T)
```

### Comparing fish lengths

#### Extracting data

``` r
#Uploading BRUVS lengths
B_Lengths <- read.csv("Data/BRUVSLengths.csv") %>% 
  tbl_df() %>% filter(!SpeciesName %in% NonFish) %>% 
  select(Length..mm., SiteCode, SpeciesName, Number) %>% 
  mutate(Length..mm. = Length..mm./10) %>% 
  rename("Species" = "SpeciesName", "Length_cm" = "Length..mm.") %>% 
  right_join(., select(BRUVEnviro, SiteCode, Bioregion) %>% distinct(), by = "SiteCode") %>% 
  filter(Length_cm > 0)
B_Lengths <- FishData %>% select(Family, Genus, ScientificName, ValidName, TrophicCat) %>% 
  rename("Species"="ScientificName") %>% 
  left_join(B_Lengths, ., by = "Species") %>% 
  select(-Species) %>% 
  rename("Species"="ValidName") %>% 
  mutate(Species = factor(Species)) %>% 
  ungroup()

#Creating function to extract length data for taxonomic groups of interest from BRUVS and UVC databases
ExtractLengths <- function(Sp){
  df <- rbind(UVC %>% filter(grepl(Sp, Species)) %>% 
  select(Bioregion, EstimatedSize, Counts, Species) %>% 
  rename("Length_cm" = "EstimatedSize") %>% 
  uncount(Counts) %>% 
  mutate(Method = "UVC"),
  B_Lengths %>% filter(grepl(Sp, Species)) %>% 
  select(Bioregion, Length_cm, Species) %>% 
  mutate(Method = "BRUV")) %>% 
    mutate(Bioregion = factor(Bioregion), Method = factor(Method))
  return(df)}

#Extracting length data 
Molfax <- ExtractLengths("Mycteroperca olfax")
Mugil <- ExtractLengths("Mugil")
Lutjanus <- ExtractLengths("Lutjanus")
```

#### KDE calculations

``` r
#Mycteroperca olfax
Molfax_figM <- kde.compare(length = Molfax$Length_cm,
            group = Molfax$Method,
            align = "no",
            nboot = 500,
            xlab = "",
            ylab = "Probability Density",
            main = "Mycteroperca olfax")
```

    ## Loading required package: sm

    ## Package 'sm', version 2.2-5.6: type help(sm) for summary information

    ## KernSmooth 2.23 loaded
    ## Copyright M. P. Wand 1997-2009

    ## [1] "Dependencies Installed:" "TRUE"                   
    ## [3] "TRUE"                    "TRUE"                   
    ## 
    ## Test of equal densities:  p-value =  0.216

``` r
Molfax_figM$graph <- Molfax_figM$graph+
  theme(plot.title = element_text(face="italic"))+
  geom_vline(xintercept = FishData$JuvLgth_m[FishData$ValidName == "Mycteroperca olfax"]*100, 
                      linetype = "dashed", color = "red")+
  annotate(geom = "text", x = c(57.5,57.5), y = c(0.065,0.075), 
           label = c(paste0("n = ", nrow(Molfax)), 
                     case_when(Molfax_figM$est$p == 0 ~ "p < 0.001",
                               TRUE ~ paste0("p = ", format(round(Molfax_figM$est$p, 2), 
                                                   nsmall = 2)))), size = 4)
Molfax_figM$graph
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
Molfax_figBio <- kde.compare(length = Molfax$Length_cm,
            group = Molfax$Bioregion,
            align = "no",
            nboot = 500,
            xlab = "",
            ylab = "Probability Density",
            main = "")
```

    ## [1] "Dependencies Installed:" "TRUE"                   
    ## [3] "TRUE"                    "TRUE"                   
    ## 
    ## Test of equal densities:  p-value =  0.052

``` r
Molfax_figBio$graph <- Molfax_figBio$graph + 
  geom_vline(xintercept = FishData$JuvLgth_m[FishData$ValidName == "Mycteroperca olfax"]*100, 
                      linetype = "dashed", color = "red")+
  annotate(geom = "text", x = c(57.5,57.5), y = c(0.065,0.075), 
           label = c(paste0("n = ", nrow(Molfax)), 
                     case_when(Molfax_figBio$est$p == 0 ~ "p < 0.001",
                               TRUE ~ paste0("p = ", format(round(Molfax_figBio$est$p, 2), 
                                                   nsmall = 2)))), size = 4)
Molfax_figBio$graph
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
#Mugil spp
Mugil_figM <- kde.compare(length = Mugil$Length_cm,
            group = Mugil$Method,
            align = "no",
            nboot = 500,
            xlab = "",
            ylab = "",
            main = "")
```

    ## [1] "Dependencies Installed:" "TRUE"                   
    ## [3] "TRUE"                    "TRUE"                   
    ## 
    ## Test of equal densities:  p-value =  0

``` r
Mugil_figM$graph <- Mugil_figM$graph + 
  labs(title = expression(paste(italic("Mugil"), " spp")))+
  geom_vline(xintercept = FishData %>% filter(Genus == "Mugil") %>% 
               summarise_at(vars(JuvLgth_m), funs(max(., na.rm = T)*100)) %>% as.numeric(), 
                      linetype = "dashed", color = "red")+
  annotate(geom = "text", x = c(72.5,72.5), y = c(0.165,0.19), 
           label = c(paste0("n = ", nrow(Mugil)), 
                     case_when(Mugil_figM$est$p == 0 ~ "p < 0.001",
                               TRUE ~ paste0("p = ", format(round(Mugil_figM$est$p, 2), 
                                                   nsmall = 2)))), size = 4)
Mugil_figM$graph
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
Mugil_figBio <- kde.compare(length = Mugil$Length_cm,
            group = Mugil$Bioregion,
            align = "no",
            nboot = 500,
            xlab = "Fork Length (cm)",
            ylab = "",
            main = "")
```

    ## [1] "Dependencies Installed:" "TRUE"                   
    ## [3] "TRUE"                    "TRUE"                   
    ## 
    ## Test of equal densities:  p-value =  0

``` r
Mugil_figBio$graph <- Mugil_figBio$graph + 
  geom_vline(xintercept = FishData %>% filter(Genus == "Mugil") %>% 
               summarise_at(vars(JuvLgth_m), funs(max(., na.rm = T)*100)) %>% as.numeric(), 
                      linetype = "dashed", color = "red")+
  annotate(geom = "text", x = c(72.5,72.5), y = c(0.175,0.20), 
           label = c(paste0("n = ", nrow(Mugil)), 
                     case_when(Mugil_figBio$est$p == 0 ~ "p < 0.001",
                               TRUE ~ paste0("p = ", format(round(Mugil_figBio$est$p, 2), 
                                                   nsmall = 2)))), size = 4)
Mugil_figBio$graph
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

``` r
#Lutjanus spp
Lutjanus_figM <- kde.compare(length = Lutjanus$Length_cm,
            group = Lutjanus$Method,
            align = "no",
            nboot = 500,
            xlab = "",
            ylab = "",
            main = "")
```

    ## [1] "Dependencies Installed:" "TRUE"                   
    ## [3] "TRUE"                    "TRUE"                   
    ## 
    ## Test of equal densities:  p-value =  0

``` r
Lutjanus_figM$graph <- Lutjanus_figM$graph + 
  labs(title = expression(paste(italic("Lutjanus"), " spp")))+
  geom_vline(xintercept = FishData %>% filter(Genus == "Lutjanus") %>% 
               summarise_at(vars(JuvLgth_m), funs(max(., na.rm = T)*100)) %>% as.numeric(), 
                      linetype = "dashed", color = "red")+
  annotate(geom = "text", x = c(85,85), y = c(0.0425,0.05), 
           label = c(paste0("n = ", nrow(Lutjanus)), 
                     case_when(Lutjanus_figM$est$p == 0 ~ "p < 0.001",
                               TRUE ~ paste0("p = ", format(round(Lutjanus_figM$est$p, 2), 
                                                   nsmall = 2)))), size = 4)
Lutjanus_figM$graph
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

``` r
Lutjanus_figBio <- kde.compare(length = Lutjanus$Length_cm,
            group = Lutjanus$Bioregion,
            align = "no",
            nboot = 500,
            xlab = "",
            ylab = "",
            main = "")
```

    ## [1] "Dependencies Installed:" "TRUE"                   
    ## [3] "TRUE"                    "TRUE"                   
    ## 
    ## Test of equal densities:  p-value =  0

``` r
Lutjanus_figBio$graph <- Lutjanus_figBio$graph + 
  geom_vline(xintercept = FishData %>% filter(Genus == "Lutjanus") %>% 
               summarise_at(vars(JuvLgth_m), funs(max(., na.rm = T)*100)) %>% as.numeric(), 
                      linetype = "dashed", color = "red")+
  annotate(geom = "text", x = c(85,85), y = c(0.15,0.175), 
           label = c(paste0("n = ", nrow(Lutjanus)), 
                     case_when(Lutjanus_figBio$est$p == 0 ~ "p < 0.001",
                               TRUE ~ paste0("p = ", format(round(Lutjanus_figBio$est$p, 2), 
                                                   nsmall = 2)))), size = 4)
Lutjanus_figBio$graph
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->

``` r
#Putting graphs together
#Method
MetFig <- ggarrange(Molfax_figM$graph, Mugil_figM$graph, Lutjanus_figM$graph, ncol = 3, 
                    common.legend = T, hjust = -1.5, vjust = 1, legend = "bottom", labels = "A")
MetFig+theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->

``` r
# MetFig <- annotate_figure(MetFig, bottom = text_grob("Method comparison", hjust = 0.5, 
#                                                      vjust = -4))
MetFig
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->

``` r
#Bioregion
BioFig <- ggarrange(Molfax_figBio$graph, Mugil_figBio$graph, Lutjanus_figBio$graph, ncol = 3, 
                    legend = "bottom", common.legend = T, labels = "B", hjust = -1.5, vjust = 1)
# BioFig <- annotate_figure(BioFig, left = text_grob("Bioregion comparison", rot = 90, hjust = 0.3))
BioFig
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-9.png)<!-- -->

``` r
#Final
FinalFig <- ggarrange(MetFig, BioFig, nrow = 2)
FinalFig
```

![](GalapagosMangroveAnalysisAllPts_files/figure-gfm/unnamed-chunk-5-10.png)<!-- -->

``` r
ggsave("Figures/LengthsKDE3.tiff", FinalFig, device = "tiff", dpi = 400, width = 22.26, 
       height = 20.73, units = "cm")
rm(Molfax_figM, Mugil_figM, Lutjanus_figM, Molfax_figBio, Mugil_figBio, Lutjanus_figBio, MetFig, BioFig)
```
