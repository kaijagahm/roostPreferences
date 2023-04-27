# Look at roost order

# Load packages
library(tidyverse)
library(sf)
library(abind)

# Load data
## Roost polygons
rp <- sf::st_read("data/roosts50_kde95_cutOffRegion.kml")

## Roost coordinates
load("data/roosts_seasons.Rda")

## Seasons, for season names (this is more data than I need but oh well)
load("data/seasons.Rda")
s <- map_chr(seasons, ~.x$seasonUnique[1])

# Get roost polygon IDs
## Make sure that each roost point intersects only one (or 0) roost (multi)polygon(s)
table(unlist(map(roosts_seasons, ~as.numeric(lengths(st_intersects(.x, rp)))))) # should show only 1s and 0s

## Get ID's of which polygon each roost location falls into
roostIDs <- map(roosts_seasons, ~as.numeric(st_intersects(x = .x, y = rp)))

roosts_seasons <- map2(roosts_seasons, roostIDs, ~.x %>% 
                         mutate(roostID = as.character(.y)) %>%
                         mutate(roostID = case_when(is.na(roostID) ~ paste(Nili_id, roost_date, sep = "_"), 
                                                    TRUE ~ roostID)))

# Simplify the data
rs <- map(roosts_seasons, ~.x %>% 
            st_drop_geometry() %>%
            dplyr::select(Nili_id, roost_date, roostID))

seqs <- map(rs, ~.x %>%
              # pivot to wider so we get the NA's in there
              pivot_wider(names_from = "roost_date", 
                          values_from = "roostID") %>%
              # pivot back to long so we keep the NAs
              pivot_longer(cols = -Nili_id, values_to = "roostID", 
                           names_to = "roost_date", values_drop_na = FALSE) %>%
              split(f = as.factor(.$Nili_id)) %>%
              map(., ~.x %>% pull(roostID)))

compare <- function(indsList, threshold = 10){
  # `threshold` is the minimum number of nights, or minimum proportion of nights, for which we must have roost data for both individuals in order to proceed with the comparison. Default is 10 nights. 
  nNights <- length(indsList[[1]])
  if(!all(map_dbl(indsList, length) == nNights)){
    stop("All elements of `indsList` must be the same length (same number of nights). NA's are fine.")
  }
  # if threshold is provided as a proportion, convert it to a number of nights
  if(threshold < 0){
    threshold <- floor(nNights*threshold)
  }else{
    threshold <- threshold # treat it as a number of nights
  }
  # create a matrix to store the data
  m <- matrix(data = NA, nrow = length(indsList), ncol = length(indsList))
  for(row in 1:length(indsList)){
    for(col in 1:length(indsList)){
      vec1 <- indsList[[row]]
      vec2 <- indsList[[col]]
      compared <- vec1 == vec2
      if(sum(!is.na(compared)) >= threshold){
        nSame <- sum(compared[!is.na(compared)])
        prop <- nSame/sum(!is.na(vec1))
        m[row, col] <- round(prop, 3)
      }else{
        m[row, col] <- NA
      }
    }
  }
  return(m)
}

n <- 50

seqsRandomized <- vector(mode = "list", length = length(seqs))
actuals <- vector(mode = "list", length = length(seqs))
for(i in 1:length(seqs)){
  dat <- seqs[[i]]
  inds <- names(dat)
  randomReps <- vector(mode = "list", length = n)
  for(j in 1:n){
    randomReps[[j]] <- map(dat, ~sample(.x, size = length(.x), replace = F))
  }
  
  actual <- compare(dat, threshold = 10) %>% 
    as.data.frame() %>% setNames(inds) %>% mutate(ID1 = inds) %>%
    pivot_longer(cols = -ID1, names_to = "ID2", values_to = "actual")
  randomized <- imap(randomReps, ~.x %>%
                       compare(., threshold = 10) %>% 
                       as.data.frame() %>%
                       setNames(inds) %>%
                       mutate(ID1 = inds) %>%
                       pivot_longer(cols = -ID1, 
                                    names_to = "ID2", 
                                    values_to = "value") %>%
                       mutate(rep = .y)) %>%
    purrr::list_rbind()
  seqsRandomized[[i]] <- randomized
  actuals[[i]] <- actual
}
# XXX start here
actuals <- map2(actuals, s, ~.x %>% mutate(season = .y)) %>%
  list_rbind()

randomized <- map2(seqsRandomized, s, ~.x %>% mutate(season = .y)) %>%
  purrr::list_rbind() %>%
  filter(ID1 != ID2) %>% # remove self
  left_join(actuals)

pvals <- randomized %>%
  group_by(ID1, ID2, actual, season) %>%
  summarize(pVal = sum(value > actual)/n()) %>%
  ungroup()

mns <- randomized %>%
  mutate(diff = actual-value) %>%
  group_by(ID1, ID2, actual, season) %>%
  summarize(mn = mean(value, na.rm = T),
            mnDiff = mean(diff, na.rm = T)) %>%
  mutate(diffMnActual = actual-mn) %>% # okay, the difference of the means and the mean of the differences is the same. Good to know. 
  ungroup()

all <- mns %>%
  left_join(pvals, by = c("ID1", "ID2", "season", "actual")) %>%
  mutate(dyad = paste(ID1, ID2, sep = "_"))

# Just look at one individual, abel, and all his associations
abel <- randomized %>%
  filter(ID1 == "abel", season == "2022_b")

(abelPrefs <- all %>%
    filter(ID1 == "abel", season == "2022_b") %>%
    ggplot()+
    geom_segment(aes(x = fct_reorder(ID2, desc(diffMnActual)), xend = fct_reorder(ID2, desc(diffMnActual)), y = mn, yend = actual))+
    geom_point(aes(x = fct_reorder(ID2, desc(diffMnActual)), y = mn), col = "skyblue")+
    geom_point(aes(x = fct_reorder(ID2, desc(diffMnActual)), y = actual), col = "red")+
    theme_classic()+
    theme(axis.text.x = element_blank())+
    ylab("Prop. nights together")+
    xlab("Other individual")+
    ggtitle("Abel"))
ggsave(abelPrefs, file = "fig/abelPrefs.png")

pamela <- randomized %>%
  filter(ID1 == "pamela", season == "2022_b")

(pamelaPrefs <- all %>%
    filter(ID1 == "pamela", season == "2022_b", !is.na(actual)) %>%
    ggplot()+
    geom_segment(aes(x = fct_reorder(ID2, desc(diffMnActual)), xend = fct_reorder(ID2, desc(diffMnActual)), y = mn, yend = actual))+
    geom_point(aes(x = fct_reorder(ID2, desc(diffMnActual)), y = mn), col = "skyblue")+
    geom_point(aes(x = fct_reorder(ID2, desc(diffMnActual)), y = actual), col = "red")+
    theme_classic()+
    theme(axis.text.x = element_blank())+
    ylab("Prop. nights together")+
    xlab("Other individual")+
    ggtitle("Pamela"))
ggsave(pamelaPrefs, file = "fig/pamelaPrefs.png")

# Ran through this with Sean. Questions:
# 1. Do individuals differ in the overall size of their preferences?
(indivPrefs <- all %>%
    group_by(ID1, season) %>%
    filter(n() > 5) %>% # only individuals that interacted with at least 5 others
    mutate(mnPref = mean(diffMnActual, na.rm = T)) %>%
    ungroup() %>%
    ggplot(aes(x = fct_reorder(ID1, desc(mnPref)), y = diffMnActual))+
    geom_boxplot(col = "gray30", fill = "gray30", outlier.size = 0.5, size = 0.5)+
    theme_classic()+
    geom_hline(yintercept = 0, col = "red", lty = 2)+
    facet_wrap(~season, nrow = 3)+ # ... kinda? but there's just a lot of variability generally.
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    ylab("Proportion differential (actual vs. mean random)")+
    xlab("Individual"))
ggsave(indivPrefs, file = "fig/indivPrefs.png", width = 9, height = 5, dpi = 800)

# 2. Do roost locations vary in busyness night to night? E.g. do all individuals prefer the same site on *certain nights* all at the same time?
rss <- map2(rs, s, ~.x %>% mutate(season = .y)) %>%
  purrr::list_rbind()

busy <- rss %>%
  group_by(season, roost_date) %>%
  mutate(indsTonight = length(unique(Nili_id))) %>%
  ungroup() %>%
  group_by(season, roostID, roost_date) %>%
  mutate(indsHereTonight = length(unique(Nili_id))) %>%
  mutate(propIndsHereTonight = indsHereTonight/indsTonight)

# proportion of individuals at the site tonight
(propIndsHereTonight <- busy %>%
    group_by(roostID) %>%
    filter(length(unique(roost_date)) > 3) %>% # only include roosts used by someone on at least 3 nights
    ggplot(aes(x = roost_date, y = propIndsHereTonight, col = roostID))+
    geom_line()+
    theme_minimal()+
    theme(legend.position = "none") +
    facet_wrap(~season, nrow = 3, scales = "free_x"))
ggsave(propIndsHereTonight, file = "fig/propIndsHereTonight.png", width = 9, height = 5, dpi = 800)

# number of individuals at the site tonight
(indsHereTonight <- busy %>%
    group_by(roostID) %>%
    filter(length(unique(roost_date)) > 3) %>% # only include roosts used by someone on at least 3 nights
    ggplot(aes(x = roost_date, y = indsHereTonight, col = roostID))+
    geom_line()+
    theme_minimal()+
    theme(legend.position = "none") +
    facet_wrap(~season, nrow = 3, scales = "free_x"))
ggsave(indsHereTonight, file = "fig/indsHereTonight.png", width = 9, height = 5, dpi = 800)

# Okay, we see that some roosts are consistently more populated than others. There are changes over time in which roosts are being used. But I don't immediately discern a pattern where different roosts keep turning into temporary hotspots, and then declining, etc. I think that's the kind of pattern that would be necessary in order for preferences to be driven according to roost desirability...?

# Of course, number of individuals at a roost in a given night isn't a direct measure of desirability. And this is all tangled up...

# 3. Do individuals' preferences for certain other individuals change through the seasons? (this may be a question that we need more data to answer, if dyads aren't present in all seasons)
# Let's start by just looking at the dyads that are present in all three seasons
durableDyads <- all %>%
  filter(!is.na(actual), !is.na(diffMnActual)) %>%
  dplyr::select(season, dyad) %>%
  distinct() %>%
  group_by(dyad) %>%
  filter(n() == 3) %>%
  pull(dyad)

ddData <- all %>%
  filter(dyad %in% durableDyads)

(dyadIncDec <- ddData %>%
    mutate(season = factor(season, levels = c("2022_b", "2022_nb", "2023_b"))) %>%
    group_by(dyad) %>%
    arrange(season, .by_group = T) %>%
    mutate("Net change" = case_when(diffMnActual[3] > diffMnActual[1] ~ "+",
                                    diffMnActual[3] < diffMnActual[1] ~ "-", 
                                    TRUE ~ NA)) %>%
    filter(!is.na(`Net change`)) %>%
    ungroup() %>%
    ggplot()+
    geom_line(aes(x = season, y = diffMnActual, group = dyad, col = `Net change`),
              linewidth = 0.1, alpha = 0.3)+
    theme_classic()+
    theme(legend.position = "none")+
    scale_color_manual(values = c("firebrick", "skyblue3"))+
    geom_hline(yintercept = 0, lty = 2)+
    theme(text = element_text(size = 18))+
    xlab("Season")+
    ylab("Proportion differential (actual vs. mean random)")+
    facet_wrap(~`Net change`))
ggsave(dyadIncDec, file = "fig/dyadIncDec.png", width = 9, height = 7, dpi = 800)

# 4. What differences are there between seasons, regardless of year? (Should )



