# pedFamilias 0.2.6

* `readFam()` gains argument `convert`, allowing the parsed Familias data to be returned without conversion to pedsuite objects.

* `readFam()` now explicitly rejects files made with the 'Familial searching' module.

* `writeFam()` now harmonises marker data across pedigree components before writing the file.

* `writeFam()` detects and stops with an error if the data contains unnamed markers.

* Fixed a bug causing URLs to be rejected when `verbose = FALSE`.

* Fixed parsing of parameters written in scientific notation.


# pedFamilias 0.2.5

This is a maintenance release fixing a few minor bugs and typos.

# pedFamilias 0.2.4

* In `readFam()` allow `path` to be an URL.
* Minor doc improvements.


# pedFamilias 0.2.2

* Fix parsing of 'extra' individuals.
* Use fallback mutation model if stabilization fails.


# pedFamilias 0.2.0

* In `readFam(..., includeParams = T)`, store `version` and `dvi` in list of parameters.
* Use `cat()` instead of `message()` for verbose info in the console.
* `readFam()` gains new parameter `deduplicate`, for removing redundant pedigrees in files created with the `dvi` module.


# pedFamilias 0.1.1

* Initial CRAN submission.
