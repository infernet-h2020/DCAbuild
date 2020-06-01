
function filter_insertions(seed::Dict)

    tmpname = tempname() #PlmDCA takes a file as input
    #tmpname = "tmp.fasta"
    tmpfile = open(tmpname, "w")
    for name in keys(seed)
        @printf(tmpfile, ">%s\n", name)
	seq = seed[name]
        s = ""
        for i in 1:length(seq)
	    if isuppercase(seq[i]) || seq[i] == '-'
                s = s * seq[i]
            end
        end
        @printf(tmpfile, "%s\n", s)
    end
    flush(tmpfile)
    return tmpname

end


function print_results(filename_ins::String, l_o::Array{Float64,1}, l_e::Array{Float64,1}, filename_par::String, PlmData::Any, ctype::Symbol, L::Int64)


    println("### Printing insertions penalties in ", filename_ins, " ###")
    file_ins = open(filename_ins, "w")
    for i in 1:L
        @printf(file_ins, "%.4f\t%.4f\n", l_o[i], l_e[i])
    end

    filep = open(filename_par, "w")
    if ctype == :amino
        q = 21
    elseif ctype == :nbase
        q = 5
    end


    println("### Printing DCA model in ", filename_par, " ###")
    for i in 1:L, j in i+1:L, a in 1:q, b in 1:q
    	@printf(filep, "J %d %d %d %d %f\n", a,b,i,j,PlmData.Jtensor[a,b,i,j])
    end
    flush(filep)
    for i in 1:L, a in 1:q
    	@printf(filep, "h %d %d %f\n", a,i,PlmData.htensor[a,i])
    end
    flush(filep)

end


function readfull(filefull; ctype::Symbol="amino")
    ffull = FastaIO.FastaReader(filefull)
    dheader = Dict{String,String}()
    for (_name,seq) in ffull
        name=String(split(_name)[1])
        if occursin('|', name)
            (aux,name) = split(name, "|")
        end
        #if occursin('/', name)
        #    (name,aux) = split(name,"/")
        #end
        if !haskey(dheader,name)
            dheader[name] = seq
        end
    end
    close(ffull)
    return dheader
end

function extract_full_seq(full::Dict, seed::Dict, Mtest::Int64)

    dheader = Dict{String, String}()

    for (_name, seq) in seed
        name = String(split(_name)[1])
        if occursin('/', name)
            (name,aux) = split(name, "/")
        end
        try
            seqfull = full[name]
            if !haskey(dheader, name)
                dheader[name] = seqfull
            end
        catch
            println("Error: ", name, " not found in full length sequences")
        end
    end

    tmpname = tempname()
    fulltmp = open(tmpname, "w")

    if Mtest == length(seed)
        for (name, seq) in dheader
            @printf(fulltmp, ">%s\n", name)
            @printf(fulltmp, "%s\n", seq)
        end
    else
        idx = rand(1:length(seed), Mtest)
        names = collect(keys(dheader))
        for i in 1:Mtest
            a = idx[i]
            name = names[a]
            seq = dheader[name]
            @printf(fulltmp, ">%s\n", name)
            @printf(fulltmp, "%s\n", seq)
        end
    end
    flush(fulltmp)

    return tmpname
end


function compute_average_dist(trueseed::Dict, tmpseed::Dict)

    i = 1
    d = zeros(length(trueseed))
    hd, tmp1, tmp2, tmp3 = 0, 0, 0, 0
    for (_name,seq) in trueseed
        L = length(seq)
        if occursin('/',_name)
            (name, hit) = split(_name, "/")
            (b,en) = split(hit, "-")
            b = parse(Int64,b)
            en = parse(Int64, en)
        else
            b = 0
            en = L
            name = _name
        end
        th = 0.10*en
        for (tmpname, tmpseq) in tmpseed
            if occursin('/', tmpname)
                (tmpname, hit) = split(tmpname, "/")
                (tmpb, tmpen) = split(hit, "-")
                tmpb = parse(Int64, tmpb)
                tmpen = parse(Int64, tmpen)
            else
                tmpb = 0
                tmpen = L
            end
            if occursin(tmpname, name)
                if abs(b-tmpb) < th && abs(en - tmpen) < th
                    hd, tmp1, tmp2, tmp3 = AlignPotts.hammingdist(seq, tmpseq)
                    continue
                end
            end
        end
        d[i] = hd / (L*1.0)
        i = i + 1
    end

    return mean(d)

end
