module DCAbuild
#using AlignPotts
using PlmDCA, FastaIO
using Statistics, Printf, DelimitedFiles, Logging
#using ExtractMacro, OffsetArrays
#export Seq, palign
#import Base.show
#using Distributed


include("insertions.jl")      # all insertions related functions
include("utils.jl")	      # i/o functions,
include("build_model.jl")     # main script for learning a Potts model, insertion penalties and gap penalties


function __init__()
	if Sys.which("mafft") === nothing
		 error("Please install mafft binaries at https://mafft.cbrc.jp/alignment/software/mafft-7.453-with-extensions-src.tgz and make the executables available to the system")
	end
	if Sys.which("hmmbuild") === nothing
		error("Please install hmmer suite at http://eddylab.org/software/hmmer/hmmer.tar.gz and make the executables available to the system")
	end
end

end #end module
