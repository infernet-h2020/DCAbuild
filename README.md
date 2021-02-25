DCAbuild
============

This package produces the set of parameters needed by [DCAlign](https://github.com/infernet-h2020/DCAlign) to align a full set of sequences.


Installation
============

The package relies on two external dependencies:

1. the HMMER suite can be installed at http://hmmer.org/ (make sure that the executables are in your binary PATH). A `make install` should do it). Everything has been tested with HMMER 3.1b1, and 3.1b2, but later versions should be (hopefully) backward compatibles.

2. The MAFFT package can be installed at https://mafft.cbrc.jp/alignment/software/mafft-7.453-with-extensions-src.tgz (make sure that the executables are in your binary PATH)

Once HMMER and MAFFT have been installed, the `DCAbuild` package can be installed either as a local package cloning the repo in a local folder, or entering in the PackageManager (typing the `]` key) the following five packages in the exact header order shown below:

```
(@v1.?) pkg> add https://github.com/carlobaldassi/GaussDCA
(@v1.?) pkg> add https://github.com/pagnani/PottsGauge
(@v1.?) pkg> add https://github.com/pagnani/PlmDCA
(@v1.?) pkg> dev https://github.com/infernet-h2020/DCAlign
(@v1.?) pkg> add https://github.com/infernet-h2020/DCAbuild
(@v1.?) julia> using DCAbuild
```

Usage
============

## Pre-processing: Get aligned seed in a proper format 

`DCAlign` receives as input parameters the Direct Couplings Analysis (DCA) model associated with the seed sequences together with a set of gap and insertions penalties learned from the seed alignment. `DCAbuild` provides the necessary input files from a well-aligned seed and the information on the insertions statistics. Depending on the seed sequences, we propose here two routines to get an alignment of the seed in FASTA and *.ins* format (the latter is the usual FASTA format where together with the aligned symbols, we show insertions as lowercase symbols).

### Pfam seed

Suppose that we need to learn the parameters describing a protein family whose seed alignment can be download from [Pfam](https://pfam.xfam.org/). We suggest to download the aligned seed from the section *Alignments* in Stockholm or FASTA format, with mixed notation "." and "-", (an example of such files is provided in `test/PF00684/PF00684_seed.txt`) and to run
```
DCAbuild.align_seed_pfam(input_file, output_file)
```
where `input_file` is the seed in Stockholm format and `output_file` is the label used to name the output file which will be `output_file.fasta` and `output_file.ins`. In the example, running

```
DCAbuild.align_seed_pfam("../test/PF00684/PF00684_seed.txt", "../test/PF00684/PF00684_seed")
```
we get `PF00684_seed.fasta` and `PF00684_seed.ins`. This routine makes use of the `hmmbuild` suite with the option `-O` which provides the annotated and possibly modified multiple sequence alignment of the input seed.

### Unaligned seed 

Suppose that we are giving a set of unaligned sequences belonging to a seed. To align them, we propose to use the `MAFFT` (with the accurate `linsi` option) package; the columns of the output multiple sequences alignment will be filtered by removing the columns that show more than a fraction `1.0 - thr_ins` of gaps (the default value is `thr_ins = 0.5`). The function which performs these tasks is
```
DCAbuild.align_seed_mafft(input_file, output_file)
```
where `input_file` contains the unaligned sequences in FASTA format and `output_file` is the label used to name the output file `output_file.fasta` and `output_file.ins`. As an example, we provide a set of synthetically generated sequences in `../test/Synthetic_prof/Train_prof_short.full`. Running

```
DCAbuild.align_seed_mafft("../test/Synthetic_prof/Train_prof_short.full", "../test/Synthetic_prof/Train_prof", thr_ins = 0.2)

```
we get `Train_prof.fasta` and `Train_prof.ins` where each column of the final MSA in `Train_prof.fasta` contains at most 20 % of gaps.

## Get the DCA model and the gap/insertion penalties

To get the set of parameters characterizing the seed sequences, we need to run `DCAbuild.build_model` which takes as input

+ `fileseed` : the multiple sequence alignment of the seed in *.ins* format
+ `filefull` : the set of unaligned sequences of the seed
+ `ctype` : a `Symbol` variable specifying the type of symbols used, `:amino` for amino-acids and `:nbase` for RNA
+ `L` : the length of the aligned sequences. This parameter is provided by `DCAbuild.align_seed_mafft` and `DCAbuild.align_seed_pfam`
This function produces:
+ `filename_ins` : a file containing the penalties of opening and extending an insertion in all the positions of the multiple sequence alignment (see [DCAlign](https://github.com/infernet-h2020/DCAlign) for further details). Default: `LambdaOpen_LambdaExt.dat`
+ `filename_par` : a file storing the DCA parameters learned from [PlmDCA](https://github.com/pagnani/PlmDCA). Default: `Parameters_PlmDCA.dat`
+ `filename_gap` : a file containing the values of the external and internal gap penalties (see [DCAlign](https://github.com/infernet-h2020/DCAlign) for further details). Default: `Gap_Ext_Int.dat`

### Example

Using the files obtained through the pre-processing and the unaligned set of seed sequences, we can learn the parameters characterizing the `PF00684` family using

```
DCAbuild.build_model("../test/PF00684/PF00684seed.ins",
                     "../test/PF00684/PF00684_full_length_sequences.fasta", :amino, 66)
```

Pay attention that this is a very long run: to infer the gap penalties, it re-aligns all the seed sequences (~ 1500) for 81 times (i.e. all possible values of the gap penalties). We can choose of accelerating the process by reducing the number of sequences used within the training of the gap penalties setting the value of `Mtest` used for this task. In this case the run would be:

```
DCAbuild.build_model("../test/PF00684/PF00684seed.ins",
                     "../test/PF00684/PF00684_full_length_sequences.fasta", :amino, 66, Mtest = Mtest)
```

## RNA sequences

The code is not yet able to cope with RNA sequences, we will provide soon an RNA counterpart of this procedure. Feel free to contact us if you are interested in using [DCAlign](https://github.com/infernet-h2020/DCAlign) for this kind of data.
