library(pedsuite)

db = forrel::NorwegianFrequencies[1:10]

# True pedigree: Father & son
ped1 = nuclearPed(fa = "AF", mo = "MO", ch = "CH") |>
  profileSim(ids = c("AF", "CH"), markers = db, seed = 1729) |>
  setMutmod(model = "stepwise", rate = list(male=0.002, female = 0.001),
            range = 0.1, rate2 = 1e-6)

# Unrelated hypothesis
ped2 = singletons(typedMembers(ped1)) |>
  transferMarkers(from = ped1, to = _)

dat = list(H1 = ped1, H2 = ped2)

# For comparison with Familias below
kinshipLR(dat)

# Write .fam file and open in Familias
fam = file.path("inst/extdata", "paternity.fam")
writeFam(dat, famfile = fam)

# Open in Familias -> Calculate LR -> Save
openFamilias(fam)
