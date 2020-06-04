# DCAbuild

This package produces the set of parameters needed by DCAlign to align a full set of sequences.


Installation
============

The package relies on 2 external dependencies:

1. the HMMER suite can be installed at http://hmmer.org/ (make sure that the executables are in your binary PATH). A `make install` should do it). Everything
has been tested with HMMER 3.1b1, and 3.1b2, but later versions should be
(hopefully) backward compatibles.

2. The MAFFT package can be installed at https://mafft.cbrc.jp/alignment/software/mafft-7.453-with-extensions-src.tgz
(make sure that the executables are in your binary PATH)

Once HMMER and MAFFT have been installed, the `DCAbuild` package can be installed either as local package cloning the repo in a local folder, or entering in the PackageManager (typing the `]` key) the follwing 4 packages in the exact header order shown below:

```
(@v1.?) pkg> add https://github.com/carlobaldassi/GaussDCA
(@v1.?) pkg> add https://github.com/pagnani/PottsGauge
(@v1.?) pkg> add https://github.com/pagnani/PlmDCA
(@v1.?) pkg> add https://github.com/anna-pa-m/DCAbuild (remember to add the correct name)
```

## Usage

1. Single-core:
```
julia> using DCAbuild
julia> DCAbuild.build_model("../test/PF00684/PF00684seed.ins",
                     "../test/PF00684/PF00684_full_length_sequences.fasta", :amino, 67)
```
2. Multi-core
Start julia from the shell with `julia -p ncore` where `ncore` is the number of
cores (tyicall `ncore ≤` number of physical core of the computer)
```
julia> @everwhere using DCAbuild
julia> DCAbuild.build_model("../test/PF00684/PF00684seed.ins",
                     "../test/PF00684/PF00684_full_length_sequences.fasta", :amino, 67)
```


This is a long run (it re-aligns all the seed sequences (~ 1500) for 81 times (i.e. all possible values of the gap penalties).

##

## Faster run for debug
Pick only `Mtest` sequences for determining the gap penalties
```
DCAbuild.build_model("../test/PF00684/PF00684seed.ins",
                     "../test/PF00684/PF00684_full_length_sequences.fasta", :amino, 67, Mtest = Mtest)
```
