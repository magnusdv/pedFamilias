#' Export data to .fam
#'
#' This function produces a `.fam` file readable by the Familias software (Egeland et al.,
#' 2000), containing all input pedigrees, their marker data and mutation models. The
#' option `openFam = TRUE` calls `openFamilias()` to open a fresh Familias session with
#' the produced file loaded.
#'
#' **Parameters**
#'
#' The `params` argument controls various Familias settings. Most entries below  may be
#' omitted, in which case suitable defaults are used. Single values are recycled over
#' markers or individuals as appropriate. The `params` argument is broadly compatible with
#' the `params` component produced by [readFam()] with `includeParams = TRUE`.
#'
#' *General parameters*
#'
#' * `version`: Familias version written to the file. Default: `"3.4.1"`.
#' * `dvi`: Logical indicating whether to write a file for the Familias DVI module.
#' Default: `FALSE`.
#'
#' *Database parameters*
#'
#' * `dbName`: Name of the marker database. Default: `"unknown"`.
#' * `dbSize`: Database size for each marker. Default: `1000`.
#' * `maf`: Minor allele frequency for each marker. Default: `0`.
#' * `theta`: Theta/kinship/Fst correction. Default: `0`.
#' * `dropoutValue`: Dropout probability for each marker. Default: `0`.
#' * `dropoutConsider`: Logical vector indicating which pedigree members should be
#' subject to dropout.
#'
#' * `dropout`: Convenience alternative to `dropoutConsider` and `dropoutValue`; either a
#' single dropout probability or a named vector identifying individuals with dropout. For
#' DVI files `dropout` has a different meaning (see below).
#'
#' *DVI parameters*
#'
#' The following probabilities pertain only to Direct/Identity matching in the DVI module:
#'
#' * `dropout`: Dropout probability. Default: `0`.
#' * `dropin`: Dropin probability. Default: `0`.
#' * `typingError`: Typing error probability. Default: `0`.
#'
#' **DVI data**
#'
#' For DVI files, set `params$dvi = TRUE` and supply a single list with the following
#' entries:
#'
#' * `pm`: A list of singletons; the victim samples.
#' * `am`: A named list of reference families. Each must be a connected `ped` object, and have
#' exactly one missing person.
#' * `missing`: A vector naming the missing persons.
#'
#' To ensure that the input data has the correct format, we recommend writing DVI files
#' through the wrapper `dvir::writeFamDVI()`.
#'
#' @param ... For ordinary files, one or several pedigrees. Each argument should be either
#'   a single `ped` object or a list of such. If the pedigrees are unnamed, they are
#'   assigned names "Ped 1", "Ped 2", etc. For DVI files, a single list with entries `pm`,
#'   `am` and `missing`.
#' @param famfile The name or path to the output file to be written. The extension ".fam"
#'   is added if missing.
#' @param params A list of Familias parameters controlling database settings, dropout and
#'   other options. See Details.
#' @param dbOnly A logical. If TRUE, no pedigree information is included; only the
#'   frequency database.
#' @param openFam A logical. If TRUE, an attempt is made to open the produced `.fam` file
#'   in an external Familias session. Only available on Windows systems with a working
#'   Familias installation.
#' @param FamiliasPath The path to the Familias executable. If empty, the following are
#'   tried in order: "Familias3.exe", "C:/Program Files (x86)/Familias3/Familias3.exe".
#' @param verbose A logical, by default TRUE.
#'
#' @return The file name is returned invisibly.
#'
#' @seealso [readFam()].
#'
#' @references Egeland et al. (2000). _Beyond traditional paternity and identification
#'   cases. Selecting the most probable pedigree._ Forensic Sci Int 110(1): 47-59.
#'
#' @examples
#'
#' # Create pedigree with 2 markers
#' x = nuclearPed() |>
#'   addMarker(geno = c("2/2", "1/3", "2/3"), alleles = 1:3,
#'             afreq = c(.3,.3,.4), name = "M1")
#'
#' # Write to .fam
#' tmp = writeFam(x, famfile = tempfile())
#'
#' # Read back in
#' y = readFam(tmp)
#'
#' stopifnot(identical(x, y))
#'
#'
#' ### With stepwise mutation model
#' x2 = setMutmod(x, model = "stepwise",
#'                rate = list(male = 0.001, female = 0.002),
#'                range = 0.1, rate2 = 0.0001)
#'
#' # Write and read
#' y2 = x2 |>
#'   writeFam(famfile = tempfile()) |>
#'   readFam()
#'
#' stopifnot(identical(x2, y2))
#'
#'
#' ### Read/write including detailed parameters
#' params = list(theta = 0.1, dbName = "myDB", dropout = c("3" = 0.01))
#' fam = writeFam(x2, famfile = tempfile(), params = params)
#'
#' dat = readFam(fam, includeParams = TRUE)
#'
#' # Pedigree is now in the `main` slot
#' stopifnot(identical(x2, dat$main))
#'
#' # The `dropout` parameter is converted to (and is equivalent to):
#' dat$params$dropoutConsider
#' dat$params$dropoutValue
#'
#'
#' ### Read/write frequency database
#'
#' # Write database as fam file
#' dbfam = writeFam(x2, famfile = tempfile(), dbOnly = TRUE)
#'
#' # Read back in: The result is a list of marker attributes
#' a = readFam(dbfam)
#'
#' # Attach to a pedigree and write to a new file
#' z = singleton(1) |> setMarkers(locusAttributes = a)
#' dbfam2 = writeFam(z, famfile = tempfile(), dbOnly = TRUE)
#'
#' stopifnot(identical(readLines(dbfam), readLines(dbfam2)))
#'
#'
#' ### DVI file
#'
#' pm = singleton("V1") |> addMarker(geno = "1/2", name = "M")
#'
#' am = nuclearPed(children = "MP") |>
#'   addMarker(geno = c("1/1", "1/2", NA), name = "M")
#'
#' x = list(pm = list(V1 = pm), am = list(F1 = am), missing = "MP")
#' writeFam(x, famfile = tempfile(fileext = ".fam"), params = list(dvi = TRUE))
#'
#' @export
writeFam = function(..., famfile = "ped.fam", params = NULL, dbOnly = FALSE,
                    openFam = FALSE, FamiliasPath = NULL, verbose = TRUE) {

  peds = list(...)

  # Safeguard against accidental filename included in `peds`
  if(length(peds) == 2 && is.character(peds[[2]]) && endsWith(peds[[2]], ".fam")) {
    famfile = peds[[2]]
    peds = peds[1]
  }

  # DVI data?
  if(isTRUE(params$dvi)) {
    if(length(peds) != 1)
      stop2("Multiple DVI inputs are not currently supported")
    return(.writeDVI(peds[[1]], famfile = famfile, params = params, dbOnly = dbOnly,
                     openFam = openFam, FamiliasPath = FamiliasPath, verbose = verbose))
  }

  # Regular file from here
  if(length(peds) == 1)
    peds = peds[[1]]
  if(is.ped(peds))
    peds = list(peds)

  # Ensure each entry is a pedlist
  peds = lapply(peds, function(p) if(is.ped(p)) list(p) else p)
  if(!all(isped <- sapply(peds, is.pedList))) {
    idx = which(!isped)[1]
    stop2(sprintf("Argument %d is not a pedigree, but: '%s'", idx,
                  class(peds[[idx]])))
  }

  # Flatten components and check marker names
  pnms = names(peds)
  ncomp = lengths(peds)
  flat = unlist(peds, recursive = FALSE, use.names = FALSE)

  MARKERS = unique.default(unlist(lapply(flat, name)))
  if(anyNA(MARKERS) || any(!nzchar(MARKERS)))
    stop2("All markers must have names")

  # Harmonise marker order and attributes
  flat = harmoniseMarkers(flat, verbose = FALSE)
  peds = split(flat, rep(seq_along(peds), ncomp))
  names(peds) = pnms

  # All unique individual names
  LABS = unique.default(unlist(lapply(peds, labels)))
  nind = length(LABS)
  nmar = length(MARKERS)

  # Check repeated individuals
  labs = unlist(lapply(flat, labels), use.names = FALSE)
  for(id in unique(labs[duplicated(labs)])) {
    cmp = flat[vapply(flat, function(x) id %in% labels(x), logical(1))]
    sx = vapply(cmp, getSex, numeric(1), id = id)
    gx = lapply(cmp, getAlleles, ids = id)

    if(!all(sx == sx[1]))
      stop2("Conflicting sex for repeated ID: ", id)
    if(!all(vapply(gx[-1], identical, logical(1), gx[[1]])))
      stop2("Conflicting genotypes for repeated ID: ", id)
  }

  # Extra param: Dropout
  dropoutConsider = params$dropoutConsider %||% FALSE
  dropoutValue = params$dropoutValue %||% 0
  if(length(dropoutConsider) == 1)
    dropoutConsider = setnames(rep_len(dropoutConsider, nind), LABS)
  if(length(dropoutValue) == 1)
    dropoutValue = setnames(rep_len(dropoutValue, nmar), MARKERS)

  # If shortcut "dropout" is used, overrule the others
  if(!is.null(dropoutInd <- params[["dropout"]])) {

    # Familias doesn't support individual dropout values
    if(length(udr <- unique.default(dropoutInd[dropoutInd > 0])) > 1)
      stop2("All nonzero dropout values must be equal: ", sort(udr))

    # Override dropoutValue
    dropoutValue[] = max(udr, 0)

    # Override params$dropoutConsider
    if(!is.null(dnms <- names(dropoutInd))) {
      if(anyNA(match(dnms, LABS)))
        stop2("Unknown ID in `dropout`: ", .mysetdiff(dnms, LABS))
      dropoutConsider[] = FALSE
      dropoutConsider[dnms] = dropoutInd > 0
    }
    else if(length(dropoutInd) == 1 && is.numeric(dropoutInd))
      dropoutConsider[] = dropoutInd > 0
  }

  # Make sure only typed individuals have dropout (otherwise Familias crashes!)
  untyped = unlist(lapply(peds, untypedMembers))
  dropoutConsider[untyped] = FALSE

  # Extra parameter: Database size for each marker
  dbSize = params$dbSize %||% 1000
  if(length(dbSize) == 1)
    dbSize = setnames(rep_len(dbSize, nmar), MARKERS)

  # Extra parameter: Minor allele frequency for each marker
  maf = params$maf %||% 0
  if(length(maf) == 1)
    maf = setnames(rep_len(maf, nmar), MARKERS)

  # Open output connection
  if(!endsWith(famfile, ".fam"))
    famfile = paste0(famfile, ".fam")
  fam = file(famfile, "w", encoding = "UTF-8")

  isClosed = FALSE
  on.exit(if(!isClosed) close(fam), add = TRUE)

  # Quick utilities
  addline = function(...) cat(..., file = fam, sep = "\n")
  quo = function(s) sprintf('"%s"', s %||% "")

  # Preamble
  version = params$version %||% "3.4.1"

  addline(quo(paste("Output from Familias, version", version)),
          quo(sprintf("(Actually produced by R/pedsuite, %s)", format(Sys.Date(), "%d %b %Y"))),
          version,
          quo(""),
          if(dbOnly) 0 else nind)

  # Individuals and genotypes ---------------------------------------------

  # Loop over indivs
  taken = rep(FALSE, length(LABS))
  names(taken) = LABS

  # Hack to skip everyone if `dbOnly = TRUE`
  if(dbOnly)
    taken[] = TRUE

  for(ped in peds) for(x in ped) for(id in labels(x)) {
    if(taken[id])
      next

    addline(quo(id),
            "#FALSE#",
            paste0(-1, if(dropoutConsider[id] > 0) " (Consider dropouts)"),
            "#FALSE#",
            ifelse(getSex(x, id) == 1, "#TRUE#", "#FALSE#"))

    idx = internalID(x, id)

    # 1-based list of triples (a1, a2, midx)
    alslist = lapply(x$MARKERS, function(m) {
      a = m[idx, ]
      if(any(a > 0)) {
        mIdx = match(attr(m, "name"), MARKERS)
        c(a, mIdx) # else NULL
      }
    })
    alsvec = unlist(alslist)

    # number of markers for this indiv
    k = length(alsvec) / 3
    addline(k)
    if(k > 0)
      addline(alsvec - 1)  # convert to 0-based!

    taken[id] = TRUE
  }

  # Pedigrees ---------------------------------------------------------------

  # Skip all pedigrees if `dbOnly = TRUE`
  if(dbOnly) {
    npeds = 0
  } else {
    npeds = length(peds)
    pednames = names(peds) %||% paste("Ped", 1:npeds)
  }

  addline('"Known relations"', 0,0,0, npeds)

  # Loop over pedigrees
  for(i in seq_len(npeds)) {
    ped = peds[[i]]
    nrels = sum(lengths(lapply(ped, nonfounders))) * 2

    addline(i - 1, quo(pednames[i]), 0, 0, nrels)

    for(x in ped) {
      nonfou = nonfounders(x, internal = TRUE)
      if(!length(nonfou))
        next

      # Parent-child relationships: 1-based index within component
      relx = rbind(x$FIDX[nonfou], nonfou,
                   x$MIDX[nonfou], nonfou, deparse.level = 0)

      # Convert to 1-based index within complete LABS
      relxx = match(x$ID[as.vector(relx)], LABS)

      # Convert to 0-based
      addline(relxx - 1)
    }
  }

  # Database -------------------------------------------------------

  addline(sprintf("#FALSE# (#Databases: 1 ;Theta/Kinship/Fst: %g )", params$theta %||% 0),
          nmar,
          "#TRUE#",
          quo(params$dbName %||% "unknown"))

  takenMark = rep(FALSE, length(MARKERS))
  names(takenMark) = MARKERS

  for(ped in peds) for(x in ped) for(m in x$MARKERS) {
    attrs = attributes(m)
    mname = attrs$name
    if(takenMark[mname])
      next

    # Extract mutation parameters
    hasMut = !is.null(attrs$mutmod)

    if(hasMut)
      mut = pedmut::getParams(attrs$mutmod, format = 1,
                              params = c("model", "rate", "range", "rate2"))
    else
      mut = data.frame(model = c("equal", "equal"), rate = 0, range = 0, rate2 = 0)

    mod = match(mut$model, c("equal", "proportional", "stepwise"))
    if(anyNA(mod))
      stop2("Unsupported mutation model for marker '", mname, "': ", unique(mut$model[is.na(mod)]))

    mod[mod == 3 & !is.na(mut$rate2)] = 5

    nals = length(attrs$alleles)

    # Write marker name and mutation parameters
    addline(quo(mname),
            mut$rate %na% 0,
            mod - 1,
            nals,
            mut$range %na% 0,
            mut$rate2 %na% 0)

    addline("#FALSE#", 0,
            sprintf("%d\t(DatabaseSize = %d , Dropout probability = %g , Minor allele frequency = %g )",
                    nals, dbSize[[mname]], dropoutValue[[mname]], maf[[mname]]))

    frvec = character(2*nals)
    frvec[2*(1:nals) - 1] = quo(attrs$alleles)
    frvec[2*(1:nals)] = attrs$afreq

    addline(frvec)
    takenMark[mname] = TRUE
  }

  close(fam)
  isClosed = TRUE
  pth = normalizePath(famfile)

  if(verbose)
    cat("Written to file:", pth, "\n")

  if(openFam) {
    if(verbose)
      cat("Trying to open in Familias\n")
    openFamilias(pth, FamiliasPath, verbose = verbose)
  }

  invisible(pth)
}

#' @rdname writeFam
#' @export
openFamilias = function(famfile = NULL, FamiliasPath = NULL, verbose = TRUE) {
  FamiliasPath = .checkFamilias(FamiliasPath)
  if(verbose)
    cat("Familias executable:", FamiliasPath, "\n")
  if(is.null(famfile))
    cmd = sprintf('"%s"', shQuote(FamiliasPath))
  else
    cmd = sprintf('"%s %s"', shQuote(FamiliasPath), shQuote(famfile))

  shell(cmd, wait = FALSE)
}

.checkFamilias = function(FamiliasPath = NULL) {
  if(Sys.info()["sysname"] != "Windows")
    stop2("Familias is only available on Windows systems.")
  if(is.null(FamiliasPath))
    FamiliasPath = c("Familias3.exe", "C:/Program Files (x86)/Familias3/Familias3.exe")
  for(p in FamiliasPath) {
    if(file.exists(p) && file.info(p)$exe != "no")
      return(p)
  }
  stop2("Could not find Familias executable: ", FamiliasPath)
}

`%na%` = function(x, y) {
  if(length(x) == 1) {
    if(is.na(x)) y else x
  }
  else {
    x[is.na(x)] = y
    x
  }
}

.internalAlleleMatrix = function(x) {
  mat = matrix(unlist(x$MARKERS),
               nrow = length(x$ID),
               dimnames = list(x$ID, NULL))
  mat
}

