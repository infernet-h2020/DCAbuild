module DCAbuild
using AlignPotts
using PlmDCA, FastaIO
using Statistics, Printf, DelimitedFiles
#using ExtractMacro, OffsetArrays
#export Seq, palign
#import Base.show
#using Distributed


include("insertions.jl")      # all insertions related functions
include("utils.jl")	      # i/o functions, 
include("build_model.jl")     # main script for learning a Potts model, insertion penalties and gap penalties

end #end module
