setMethodS3("getChromosomes", "PSCBS", function(this, ...) {
  segs <- this$output;
  chromosomes <- sort(unique(segs$chromosome), na.last=TRUE);

  # Drop NA dividers
  if (length(chromosomes) > 1) {
    chromosomes <- chromosomes[!is.na(chromosomes)];
  }

  chromosomes;
})

setMethodS3("nbrOfChromosomes", "PSCBS", function(this, ...) {
  length(getChromosomes(this, ...));
})

setMethodS3("nbrOfSegments", "PSCBS", function(this, ...) {
  segs <- this$output;
  nrow(segs);
})

setMethodS3("extractByChromosome", "PSCBS", function(x, chromosome, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chromosome':
  chromosome <- Arguments$getInteger(chromosome, disallow=c("NaN", "Inf"));

  
  chromosomes <- getChromosomes(this);

  res <- this;

  fields <- c("data", "output");
  for (ff in seq(along=fields)) {
    field <- fields[[ff]];
    thisFF <- this[[field]];
    keep <- is.element(thisFF$chromosome, chromosome);
    thisFF <- thisFF[keep,,drop=FALSE];
    res[[field]] <- thisFF;
  } # for (ff ...)

  res;
})


setMethodS3("append", "PSCBS", function(x, other, addSplit=TRUE, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, "PSCBS");

  # Argument 'addSplit':
  addSplit <- Arguments$getLogical(addSplit);


  # Allocate results
  res <- this;

  fields <- c("data", "output");
  for (ff in seq(along=fields)) {
    field <- fields[[ff]];
    thisFF <- this[[field]];
    otherFF <- other[[field]];

    # Sanity check
    if (ncol(thisFF) != ncol(otherFF)) {
      throw(sprintf("Cannot merge %s objects. Arguments 'this' and 'other' has different number of columns in field '%s': %s != %s", class(this)[1], field, ncol(thisFF), ncol(otherFF)));
    }

    resFF <- thisFF;
    if (addSplit) {
      resFF <- rbind(resFF, NA);
    }
    resFF <- rbind(resFF, otherFF);

    res[[field]] <- resFF;
  } # for (ff ...)

  res;
})



############################################################################
# HISTORY:
# 2010-09-26
# o getChromosomes() no longer returns NA divers.
# 2010-09-24
# o Added append() and more for PSCBS objects.
############################################################################
