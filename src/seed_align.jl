## Pipeline file to build an alignment from a full length sequence
struct SimpleAlign
    Header::Vector{String}
    Seqs::Vector{String}
    function SimpleAlign(H::Vector{String},S::Vector{String})
        l = length(S[1])
        length(H) == length(S) || throw(DimensionMismatch("header and seqs have different lenghts"))
        for i in eachindex(S)
            length(S[i]) == l || throw(DimensionMismatch("not an alignment"))
        end
        new(H,S)
    end
end

SimpleAlign(d::SummaryAlign) = SimpleAlign(d.Header, d.Seqs)

function idx_insert(d::SummaryAlign, thr_ins)
    gapfreq =  vec(1.0 .- sum(reshape(d.Pi,(d.q-1,d.N)),dims=1))
    return findall(gapfreq .> thr_ins)
end

function insertify(d::SummaryAlign, idx)
    headers, seqs = d.Header,d.Seqs
    newseqs = similar(seqs)
    for i in eachindex(seqs)
        charseq = collect(seqs[i])
        for j in idx
            _c = charseq[j]
            charseq[j] = isletter(_c) ? lowercase(_c) : '.'
        end
        newseqs[i] = String(charseq)
    end
    SimpleAlign(headers,newseqs)
end



function writefasta(fou::String, d::Dict{String,String})
    FastaWriter(fou) do fw
        for (k,v) in d
            writeentry(fw,k,v)
        end
    end
end

function writefasta(fou::String,d::SimpleAlign)
    FastaWriter(fou) do fw
        h,s = d.Header,d.Seqs
        for i in eachindex(h)
            writeentry(fw,h[i],s[i])
        end
    end
end

"""
    function create_aligned_input(
            fin::String,
            fou::String;
            nthreads::Integer = 8,
            thr_ins::Real = 1.0
        )
Read `fin` full-length seed in fasta format and produces `fou.fasta` alignment
and `fou.hmm` hmmer hmm. The pipeline runs `mafft-linsi`, then constructs a hmm
through `hmmbuild` and realign `fin` with `hmmalign`.
"""
function create_aligned_input(
    fin::String,
    fou::String;
    nthreads::Integer = 8,
    thr_ins::Real = 1.0,
)
    mktempdir() do tmpdir
        @info "working in $tmpdir"
        linsiout = joinpath(tmpdir, "all.linsi.fasta")
        @info "Aligning $fin with mafft-linsi. This may take a while ... "
        linsi(fin, linsiout, nthreads = nthreads)
        @info "alignment done!"
        rescor = summary_align(linsiout)
        tmp_align = if thr_ins < 1.0 # make insertion by hand
            idxins = idx_insert(rescor, thr_ins)
            insertify(rescor, idxins)
        else
            SimpleAlign(rescor)
        end
        file_tmp_align = joinpath(tmpdir, "tmp.linsi.fasta")
        writefasta(file_tmp_align, tmp_align)
        hmmbuild(file_tmp_align, fou * ".hmm", nthreads = 1)
        file_tmp_stk = joinpath(tmpdir, "tmp.linsi.stk")
        hmmalign(fou * ".hmm", file_tmp_align, file_tmp_stk)
        d = read_stockholm(file_tmp_stk)
        writefasta(fou * ".fasta", d)
    end
end
