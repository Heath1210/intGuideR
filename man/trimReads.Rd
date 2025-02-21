% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trimReads.R
\name{trimReads}
\alias{trimReads}
\title{Reads filtration and trimming}
\usage{
trimReads(
  path2fq1,
  path2fq2,
  path2mg,
  outdir = NULL,
  dsODN1 = "ATATGTTAATAACGGTAT",
  dsODN2 = "ATGACAACTCAATTAAAC",
  linker = "CTCCGCTTAAGGGACT"
)
}
\arguments{
\item{path2fq1}{the file path to fq1 file.}

\item{path2fq2}{the file path to fq2 file.}

\item{path2mg}{the file path to merged fastq file.}

\item{outdir}{the file fold where output files put in}

\item{dsODN1}{the sequence of dsODN closest to integrated host genome with
18bp length recommended (forward)}

\item{dsODN2}{the sequence of dsODN closest to integrated host genome with
18bp length recommended (reverse)}

\item{linker}{the sequence of adaptor linker at the end, near the genome part.
Default linker is modified from INSPIIRED pipeline.
If you use other linker, change to corresponding sequence.}
}
\description{
Filter reads matching primer and linker, dump others.
   Adjust reads to same direction that linker is on the left side.
   Trim off linker and primer sequence followed by adding UMI, barcode
   to reads ID which has been simplied but informative enough to
   distinguish each reads. Save fasta files for further analysis.
}
