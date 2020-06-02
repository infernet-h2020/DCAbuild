# DCAbuild

This package produces the set of parameters needed by DCAlign to align a full set of sequences.


Installation
============

The package relies on 2 external dependencies:

1. the HMMER suite can be installed at http://hmmer.org/ (make sure that the executables
are in your binary PATH). A `make install` should do it). Everything has
been tested with HMMER 3.1b1, and 3.1b2, but later versions should be
(hopefully) backward compatibles.

2. The MAFFT package can be installed at https://mafft.cbrc.jp/alignment/software/mafft-7.453-with-extensions-src.tgz
(make sure that the executables are in your binary PATH)

Once HMMER has been installed, the `Epistasis` package can be installed as a
local package cloning the repo in a local folder:

```
include("DCAbuild.jl")

```
## Usage
```
DCAbuild.build_model("../test/PF00684/PF00684seed.ins",
                     "../test/PF00684/PF00684_full_length_sequences.fasta", :amino, 67)
```

This is a long run (it re-aligns all the seed sequences (~ 1500) for 81 times (i.e. all possible values of the gap penalties).

## Faster run for debug
Pick only `Mtest` sequences for determining the gap penalties
```
DCAbuild.build_model("../test/PF00684/PF00684seed.ins",
                     "../test/PF00684/PF00684_full_length_sequences.fasta", :amino, 67, Mtest = Mtest)
```
