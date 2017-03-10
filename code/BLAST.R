# modified by Joshua Bloom 02232017

#######################################################################
# rBLAST - Interface to BLAST
# Copyright (C) 2015 Michael Hahsler and Anurag Nagar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

.findExecutable <- function(exe, interactive=TRUE) {
  path <- Sys.which(exe)
  if(all(path=="")) {
    if(interactive) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    return(character(0))
  }

  path[which(path!="")[1]]
}

makeblast  <- function(db = NULL) {
  if(is.null(db)) stop("No BLAST database specified!")
  db <- file.path(normalizePath(dirname(db)), basename(db))
  if(length(Sys.glob(paste(db, "*", sep="")))<1) stop("BLAST database does not exit!")

  structure(list(db = db))
}

printBLAST <- function(x, info=TRUE, ...) {
  cat("BLAST Database\nLocation:", x$db, "\n")

  if(info) {
    out <- system(paste(.findExecutable("blastdbcmd"), "-db", x$db,
      "-info"), intern=TRUE)
    cat(paste(out, collapse="\n"))
    cat("\n")
  }
}

blast_help <- function() {
  system(paste(.findExecutable(c("blastn")),
    "-help"))
}

#
predictBLAST <- function(object, newdata, BLAST_args="", custom_format ="",
  ...) {

  db <- object$db
  x <- newdata

  ## get temp files and change working directory
  wd <- tempdir()
  dir <- getwd()
  temp_file <- basename(tempfile(tmpdir = wd))
  on.exit({
    #cat(temp_file, "\n")
    file.remove(Sys.glob(paste(temp_file, "*", sep="")))
    setwd(dir)
  })
  setwd(wd)

  infile <- paste(temp_file, ".fasta", sep="")
  outfile <- paste(temp_file, "_BLAST_out.txt", sep="")

  writeXStringSet(x, infile, append=FALSE, format="fasta")

  system(paste(.findExecutable("blastn"), "-db", db,
    "-query", infile, "-out", outfile, '-outfmt "10 std', custom_format,
    '"', BLAST_args))

  c_names <- c("QueryID",  "SubjectID", "Perc.Ident",
      "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
      "S.start", "S.end", "E", "Bits")

  ## rdp output column names
  if(custom_format == "") {
  }else{
    c_names <- c(c_names, unlist(strsplit(custom_format, split = " +")))
  }

  ## read and parse rdp output
  if(is(try(cl_tab <- read.table(outfile, sep=","), silent=TRUE), "try-error")) {
    warning("BLAST did not return a match!")
    cl_tab <- data.frame(matrix(ncol=length(c_names), nrow=0))
  }

  if(ncol(cl_tab) != length(c_names)) stop("Problem with format (e.g., custom_format)!")
  colnames(cl_tab) <- c_names

  cl_tab
}



