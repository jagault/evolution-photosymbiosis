library(here)
library(data.table)

# Set working directory
setwd(here("analysis_rerun/results/stree/species_table"))

st <- fread("stree-species-table.csv")

# View(st)

st[`Observed State` == "Z", .N]
st[`Observed State` == "AZ", .N]


st[`Zoox Volatile` > 0.75, .N]
st[`Zoox Volatile` > 0.75, Species]
st[`Zoox Volatile` > 0.75, unique(Family)]

st[`Zoox Stable` > 0.75, .N]
st[`Zoox Stable` > 0.75, Species]
st[`Zoox Stable` > 0.75, unique(Family)]

st[`Azoox Volatile` > 0.75, .N]
st[`Azoox Volatile` > 0.75, Species]
st[`Azoox Volatile` > 0.75, unique(Family)]

st[`Azoox Labile` > 0.75, .N]
st[`Azoox Labile` > 0.75, unique(Family)]

st[`Azoox Stable` > 0.75, .N]
st[`Azoox Stable` > 0.75, unique(Family)]