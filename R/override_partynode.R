partynode = function (id, split = NULL, kids = NULL, surrogates = NULL, info = NULL)
{
  # if (!is.integer(id) || length(id) != 1) {
  #   id <- as.integer(id0 <- id)
  #   if (any(is.na(id)) || !isTRUE(all.equal(id0, id)) ||
  #       length(id) != 1)
  #     stop(sQuote("id"), " ", "must be a single integer")
  # }
  if (is.null(split) != is.null(kids)) {
    stop(sQuote("split"), " ", "and", " ", sQuote("kids"),
         " ", "must either both be specified or unspecified")
  }
  if (!is.null(kids)) {
    if (!(is.list(kids) && all(sapply(kids, inherits, "partynode"))) ||
        length(kids) < 2)
      stop(sQuote("kids"), " ", "must be an integer vector or a list of",
           " ", sQuote("partynode"), " ", "objects")
  }
  if (!is.null(surrogates)) {
    if (!is.list(surrogates) || any(!sapply(surrogates, inherits,
                                            "partysplit")))
      stop(sQuote("split"), " ", "is not a list of", " ",
           sQuote("partysplit"), " ", "objects")
  }
  node <- list(id = id, split = split, kids = kids, surrogates = surrogates,
               info = info)
  class(node) <- "partynode"
  return(node)
}
