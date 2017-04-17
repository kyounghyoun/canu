
# Canu / minialign

This is a fork of [Canu](https://github.com/marbl/canu) assembler. See the original repository for information of the Canu pipeline. The original README.md is also reserved in README.canu.md.

## Overview

[Minialign](https://github.com/ocxtal/minialign) long-read alignment tool is experimentally incorporated in the Canu assembler pipeline as an all-versus-all overlapper. The precise alignment calculation stage in minialign is expected to help the whole pipeline generate better consensus sequences as reported for Canu / BLASR combination.

## Installation

```bash
# install minialign
$ git clone https://github.com/ocxtal/minialign
$ cd minialign && make
$ make install PREFIX=$PREFIX
# install canu
$ git clone https://github.com/ocxtal/canu
$ cd canu/src && make
$ cp ../Linux-amd64/bin/* $PREFIX/bin		# canu binaries and minialign must be installed in the same bin directory
```

## Run

One of `corOverlapper=minialign`, `obtOverlapper=minialign`, or `utgOverlapper=minialign` is specified to invoke the minialign pipeline. By the default, `(w, k) = (5, 14)` is used in the correction and `(w, k) = (10, 16)` in the trimming and unitigging stages. To change these parameters, pass `{tag}WindowSize` and `{tag}MerSize` as options.

```bash
# example
$ canu -correct -p asm -d asm genomeSize=4.7m corOverlapper=minialign -pacbio-raw reads.fa
$ canu -trim-assemble -p asm -d asm genomeSize=4.7m obtOverlapper=minialign utgOverlapper=minialign -pacbio-corrected asm/correctedReads.fasta
```

## Results

Stats and a dotplot of a D.melanogaster sample are shown below.

