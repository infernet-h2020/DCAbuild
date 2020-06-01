# DCAbuild

This package produces the set of parameters needed by DCAlign to align a full set of sequences. 

## Install (to be fixed)
Enter in `src/` and within `julia` write:
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

