# take unaligned seed -> make a msa


function build_model(fileseed::String, filefull::String, ctype::Symbol, L::Int64;
		    filename_ins::String="LambdaOpen_LambdaExt.dat",
		    filename_par::String="Parameters_PlmDCA.dat",
		    filename_gap::String="Gap_Ext_Int.dat",
		    Mtest::Int64=0,
		    verbose::Bool=true)

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
	seed = DCAlign.readfull(fileseed, ctype=ctype, pos = true)
	println("### Inferring insertions penalties ###")
	l_o, l_e = infer_ins_pen(seed, L)

	println("### Inferring a Potts model using PlmDCA ###")
	aligntmp = filter_insertions(seed)
	println("Temporary FASTA alignment in ", aligntmp)
	PlmData = PlmDCA.plmdca(aligntmp, theta=0.20)


	println("### Finding gap penalties ###")
	println("WARNING: Reasonable values are obtained when using many (> 500) sequences")

	full = DCAlign.readfull(filefull, ctype=ctype, pos = print_pos)
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

	print_results(filename_ins, l_o, l_e, filename_par, PlmData, ctype, L, filename_gap, 0.0, 0.0) # tmp
	mu = 0.00:0.50:4.00
	muint = 0.00:0.50:4.00
	d = zeros(length(mu),length(muint))
	aseed = DCAlign.readfull(aligntmp, ctype=ctype, pos = true)
	for a in 1:length(mu)
		for b in 1:length(muint)
			println("#### Aligning ", Mtest, " sequences using (μext, μint) = (", mu[a], ", ", muint[b],")")
			filename_out = tempname()
			filename_flag = tempname()
			DCAlign.align_all(q, L, filename_par, fulltmp, filename_ins, mu[a], muint[b]; typel=:plm, filename_flag=filename_flag, filename_align=aligntmp, filename_ins=fileseed, filename_out=filename_out, verbose = verbose, maxiter = 300)
			tmpseed = DCAlign.readfull(filename_out, ctype=ctype, pos = true)
			d[a,b] = compute_average_dist(aseed, tmpseed)
			rm(filename_out)
			rm(filename_flag)
		end
	end
	#println(d)
	aux = argmin(d)
	muext_best = mu[aux[1]]
	muint_best = muint[aux[2]]
	print_results(filename_ins, l_o, l_e, filename_par, PlmData, ctype, L, filename_gap, muext_best, muint_best)
	rm(fulltmp)
	rm(aligntmp)
	println("Done!")

end
