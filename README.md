
# Canu / minialign

Note: This is a fork of [Canu](https://github.com/marbl/canu) assembler. See the original repository for information of the Canu pipeline. The original README.md is also preserved in README.canu.md.

## Overview

[Minialign](https://github.com/ocxtal/minialign) long-read alignment tool is experimentally incorporated in the Canu assembler pipeline as an all-versus-all overlapper. The precise alignment calculation of the minialign is expected to reduce false positives in the overlapping stage and help the whole pipeline generate better consensus sequences.

## Installation

```bash
# install minialign
$ git clone https://github.com/ocxtal/minialign	# minialign-0.5.1 or later is required
$ cd minialign && make
$ make install PREFIX=$PREFIX
# install canu
$ git clone https://github.com/ocxtal/canu
$ cd canu/src && make
$ cp ../Linux-amd64/bin/* $PREFIX/bin		# canu and minialign binaries must be installed in the same bin directory
```

## Run

One of `corOverlapper=minialign`, `obtOverlapper=minialign`, or `utgOverlapper=minialign` must be specified to invoke minialign in the pipeline. By the default, `(w, k) = (5, 14)` is used in the correction and `(10, 16)` in the trimming and unitigging stages. To change these parameters, pass `{tag}WindowSize` and `{tag}MerSize` as options.

```bash
# example
$ canu -correct -p asm -d asm genomeSize=4.7m corOverlapper=minialign -pacbio-raw reads.fa
$ canu -trim-assemble -p asm -d asm genomeSize=4.7m obtOverlapper=minialign utgOverlapper=minialign -pacbio-corrected asm/correctedReads.fasta
```

## Results

Not yet.

