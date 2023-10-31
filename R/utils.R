stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a single number, with optional range constraints
isNumber = function(x, minimum = NA, maximum = NA) {
  isTRUE(length(x) == 1 &&
           is.numeric(x) &&
           (is.na(minimum) || x >= minimum) &&
           (is.na(maximum) || x <= maximum))
}

# Faster alternative to suppressWarnings(as.numeric(v))
# NB: doesn't catch scientific notation 1e-4.
asNum = function(v) {
  num = grep("[^-0-9.]", v, invert = TRUE)
  u = rep(NA_real_, length(v))
  u[num] = as.numeric(v[num])
  u
}

# Faster alternative to suppressWarnings(as.integer(v))
# NB: doesn't catch scientific notation 1e+12.
asInt = function(v) {
  num = grep("[^-0-9]", v, invert = TRUE)
  u = rep(NA_real_, length(v))
  u[num] = as.integer(v[num])
  u
}

# round + toString
rst = function(v, digits = 3)
  toString(round(v, digits))

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

pluralise = function(noun, n) {
  if(n == 1) noun else sprintf("%ss", noun)
}

setnames = function(x, nms) {
  names(x) = nms
  x
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}


ftime = function(st, digits = 3) {
  format(Sys.time() - st, digits = digits)
}
