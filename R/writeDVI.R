.writeDVI = function(x, famfile = "dvi.fam", params = NULL, dbOnly = FALSE,
                     openFam = FALSE, FamiliasPath = NULL, verbose = TRUE) {
  if(!is.list(x) || !all(c("am", "pm", "missing") %in% names(x)))
    stop2("Input must be a list with entries 'am', 'pm' and 'missing'")

  pm = x$pm
  am = x$am

  if(!is.list(pm) || !all(vapply(pm, is.singleton, logical(1))))
    stop2("Non-singleton PM samples are not currently supported")
  if(!is.list(am) || !all(vapply(am, is.ped, logical(1))))
    stop2("Disconnected AM families are not currently supported")

  peds = c(pm, am)
  if(!length(peds))
    stop2("Empty DVI dataset")

  mnames = unlist(lapply(peds, name), use.names = FALSE)
  if(anyNA(mnames) || any(!nzchar(mnames)))
    stop2("All markers must have names")

  # DVI error settings are stored in params; remaining settings go to writeFam
  params = params %||% list()
  dropout = params$dropout %||% 0
  dropin = params$dropin %||% 0
  typingError = params$typingError %||% 0
  params$version = sub(" \\(DVI\\)$", "", params$version %||% "3.4.1")

  if(!isNumber(dropout, minimum = 0, maximum = 1-1e-10))
    stop2("`dropout` must be a number in the interval [0,1) : ", dropout)
  if(!isNumber(dropin, minimum = 0, maximum = 1-1e-10))
    stop2("`dropin` must be a number in the interval [0,1) : ", dropin)
  if(!isNumber(typingError, minimum = 0, maximum = 1-1e-10))
    stop2("`typingError` must be a number in the interval [0,1) : ", typingError)

  dbParams = params
  dbParams[c("dvi", "dropin", "dropout", "typingError")] = NULL

  # Use the same marker harmonisation as writeFam, without merging person names
  peds = harmoniseMarkers(peds, verbose = FALSE)
  pm = peds[seq_along(pm)]
  am = peds[length(pm) + seq_along(am)]

  if(dbOnly)
    return(writeFam(peds[[1]], famfile = famfile, params = dbParams, dbOnly = TRUE,
                    openFam = openFam, FamiliasPath = FamiliasPath, verbose = verbose))

  if(!length(pm))
    stop2("DVI data without PM samples are not supported")
  if(length(am) && (is.null(names(am)) || any(!nzchar(names(am)))))
    stop2("All AM families must be named")

  miss = lapply(am, function(p) intersect(p$ID, x$missing))
  if(any(lengths(miss) != 1L) || length(x$missing) != length(am))
    stop2("All AM families must have exactly one missing person")

  # Reuse the existing database writer, then mark the file as DVI
  file = writeFam(peds[[1]], famfile = famfile, params = dbParams,
                  dbOnly = TRUE, verbose = FALSE)
  raw = readLines(file, encoding = "UTF-8")
  raw[3] = paste(raw[3], "(DVI)")
  writeLines(raw, file, useBytes = TRUE)

  # Append the DVI section
  fam = base::file(file, "a", encoding = "UTF-8")
  on.exit(close(fam))
  addline = function(...) cat(..., file = fam, sep = "\n")

  # PM and AM persons differ only in the bracket depth of their headings
  addPerson = function(p, idx, num, depth) {
    if(!p$SEX[idx] %in% 1:2)
      stop2("Unknown sex is not currently supported: ", p$ID[idx])

    # DVI genotypes use allele labels, not the indices used in the main section
    genos = lapply(p$MARKERS, function(m) {
      a = m[idx, ]
      if(!any(a > 0)) return(NULL)
      if(any(a == 0))
        stop2("Half-typed genotypes are not currently supported: ", p$ID[idx])
      paste0(c("SystemName= ", "Allele1= ", "Allele2= "),
             c(attr(m, "name"), attr(m, "alleles")[a]))
    })

    # Omit fully missing genotypes and count the remaining markers
    title = sprintf("Person %d%s", num, if(depth == 3) " " else "")
    addline(paste0(strrep("[", depth), title, strrep("]", depth)),
            paste0("Name= ", p$ID[idx]),
            paste0("Gender= ", c("Male", "Female")[p$SEX[idx]]),
            paste0(strrep("[", depth + 1), "DNA data", strrep("]", depth + 1)),
            paste0("nSystems= ", sum(lengths(genos) > 0)),
            unlist(genos, use.names = FALSE))
  }

  # Unidentified persons and one empty note per reference family
  addline("", "[DVI]", "[[Unidentified persons]]", paste0("nPersons= ", length(pm)))
  for(i in seq_along(pm))
    addPerson(pm[[i]], idx = 1L, num = i, depth = 3L)

  addline("[[Notes]]", paste0("nFamilies= ", length(am)), rep("", length(am)),
          "[[Reference Families]]", paste0("nFamilies= ", length(am)))

  # Missing person occurs in relationships, but not in the Persons section
  for(i in seq_along(am)) {
    p = am[[i]]
    amName = names(am)[i]

    idx = which(p$ID != miss[[i]])

    addline(sprintf("[[[Family %d ]]]", i), paste0("Name= ", amName),
            "[[[[Persons]]]]", paste0("nPersons= ", length(idx)))

    for(j in seq_along(idx))
      addPerson(p, idx = idx[j], num = j, depth = 5L)

    # List parent-child pairs after replacing the missing person's label
    labs = p$ID
    labs[labs == miss[[i]]] = "Missing person"
    nonfou = nonfounders(p, internal = TRUE)
    rel = rbind(p$FIDX[nonfou], nonfou, p$MIDX[nonfou], nonfou, deparse.level = 0)
    rel = paste0(rep(c("Parent= ", "Child= "), length.out = length(rel)),
                 labs[as.vector(rel)])
    nrel = paste0("nRelations= ", length(rel) / 2)

    # Familias format contain two copies with different pedigree names
    addline("[[[[Pedigrees]]]]", "nPedigrees= 2",
            "[[[[[Pedigree 1]]]]]", "Name= Reference pedigree", nrel, rel,
            "[[[[[Pedigree 2]]]]]", paste0("Name= ", amName), nrel, rel)
  }

  # Further options (mostly hardcoded for now)
  addline("[Prior]", "UseList= Yes", paste0("Size= ", length(pm)),
          "[Advanced]", "[[Direct/Identity]]", paste0("Dropin= ", dropin),
          paste0("Dropout= ", dropout), paste0("Error= ", typingError),
          "[[Other]]", "NumberOfDecimals= 10",
          "[[DVI]]", "QuickSearch= Enabled", "AllowedMismatches= 5",
          "[[Simulation]]", "SaveGenotypes= Disabled", "SaveComplete= Enabled",
          "ForceMAF= Disabled")

  close(fam)
  on.exit(NULL)

  if(verbose)
    cat("Written to file:", file, "\n")

  if(openFam) {
    if(verbose)
      cat("Trying to open in Familias\n")
    openFamilias(file, FamiliasPath, verbose = verbose)
  }

  invisible(file)
}
