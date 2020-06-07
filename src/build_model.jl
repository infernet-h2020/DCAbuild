# take unaligned seed -> make a msa


function build_model(fileseed::String, filefull::String, ctype::Symbol, L::Int64;
		    filename_ins::String="LambdaOpen_LambdaExt.dat",
		    filename_par::String="Parameters_PlmDCA.dat",
		    Mtest::Int64=0)

	if ctype != :amino && ctype != :nbase
		error("Wrong second argument: choose between :amino and :nbase")
		return
	end
	println("### Reading seed ###")
	if ctype == :nbase
		print_pos = true
	else
		print_pos = false
	end
	seed = AlignPotts.readfull(fileseed, ctype=ctype, pos = true)
	println("### Inferring insertions penalties ###")
	l_o, l_e = infer_ins_pen(seed, L)

	println("### Inferring a Potts model using PlmDCA ###")
	aligntmp = filter_insertions(seed)
	println("Temporary FASTA alignment in ", aligntmp)
	PlmData = PlmDCA.plmdca(aligntmp, theta=0.2)
	print_results(filename_ins, l_o, l_e, filename_par, PlmData, ctype, L)

	println("### Finding gap penalties ###")
	println("WARNING: Reasonable values are obtained when using many (> 500) sequences")

	full = AlignPotts.readfull(filefull, ctype=ctype, pos = print_pos)
	if Mtest == 0
		Mtest = length(seed)
		println("Using all seed sequences to get the gap penalties")
	else
		println("Using ", Mtest, " out of ", length(seed), " to get the gap penalties")
	end
	fulltmp = extract_full_seq(full, seed, Mtest, ctype=ctype)
	println("Temporary full length sequences in ", fulltmp)

	if ctype == :amino
		q = 21
	elseif ctype == :nbase
		q = 5
	end

	# mu = 0.00:0.50:4.00
	# muint = 0.00:0.50:4.00
	mu = 0.00:1.00:3.00
	muint = 0.00:1.00:3.00
	# mu = parameter for external gap
	# muint = internal gaps
	d = zeros(length(mu),length(muint)) |> SharedArray
	aseed = AlignPotts.readfull(aligntmp, ctype=ctype, pos = true)
	@sync @distributed for a in 1:length(mu)
		for b in 1:length(muint)
			filename_out = tempname()
			filename_flag = tempname()
			AlignPotts.align_all(q, L, filename_par, fulltmp, filename_ins, mu[a], muint[b]; typel=:plm, filename_flag=filename_flag, filename_align=aligntmp, filename_ins=fileseed, filename_out=filename_out)
			tmpseed = AlignPotts.readfull(filename_out, ctype=ctype, pos = true)
			d[a,b] = compute_average_dist(aseed, tmpseed)
			rm(filename_out)
			rm(filename_flag)
		end
	end
	display(d)
	rm(fulltmp)
	rm(aligntmp)
end
