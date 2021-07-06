library(here)
library(data.table)

# Set working directory
setwd(here("analysis_rerun/results/mtree/species_table"))

st <- fread("mtree-species-table.csv")

# View(st)

st[`Observed State` == "Zooxanthellate", .N]
st[`Observed State` == "Azooxanthellate", .N]



st[`Zoox Stable` > 0.75, .N]
st[`Zoox Stable` > 0.75, Species]
st[`Zoox Stable` > 0.75, unique(Family)]

st[`Zoox Labile` > 0.75, .N]
st[`Zoox Labile` > 0.75, Species]
st[`Zoox Labile` > 0.75, unique(Family)]


st[`Azoox Stable` > 0.75, .N]
st[`Azoox Stable` > 0.75, unique(Family)]

st[`Azoox Labile` > 0.75, .N]
st[`Azoox Labile` > 0.75, unique(Family)]


