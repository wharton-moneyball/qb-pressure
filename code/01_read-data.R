library(tidyverse)
library(arrow)

# List of tracking files for weeks 1-8
tracking_files <- list.files(
  path = "data/bdb", 
  pattern = "week[1-8]\\.csv", 
  full.names = TRUE
)

# Read all tracking data and combine
# Using map_df to combine into a single dataframe
tracking_all <- tracking_files %>%
  map_df(~read_csv(.x, show_col_types = FALSE))

# Save to processed/ folder as parquet
write_parquet(tracking_all, "data/processed/bdb-tracking-all.parquet")

# Print summary for confirmation
glimpse(tracking_all)
