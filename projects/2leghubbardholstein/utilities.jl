using ITensors.HDF5

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames 
        n = Name(fn)
        d = getfield(struc, fn)

        ## DEBUGGING ## 
        if n=="optimized_basis"#typeof(d) == Vector{ITensor}
            try 
                h5open(path, "r+") do file
                    for i in 1:length(d)
                        write(file, n*"_"*string(i), d[i])
                    end
                end
            catch e
                if !isfile(path)
                    h5open(path, "w") do file
                        for i in 1:length(d)
                            write(file, n*"_"*string(i), d[i])
                        end
                    end
                else
                    h5open(path, "r+") do file
                        for i in 1:length(d)
                            if haskey(file, n*"_"*string(i))
                                delete_object(file, n*"_"*string(i))
                                file[n*"_"*string(i)] = d[i]
                            end
                        end
                    end
                end
            end
        else
            try 
                h5open(path, "r+") do file
                    write(file, n, d)
                end
            catch
                if !isfile(path)
                    h5open(path, "w") do file
                        write(file, n, d) 
                    end
                else
                    h5open(path, "r+") do file
                        if haskey(file, n)
                            delete_object(file, n)
                            file[n] = d
                        end
                    end
                end
            end
        end
    end
end

parse_numbers(s, delimiter) = parse(Float64, split(s, delimiter, keepempty=false)[end])

function read_ITensor_vector(d, f, prefix, delimiter)
    # Read in the gates 
    Ts = []
    idx = []
    for key in keys(d)
        if startswith(key, prefix)
            T = read(f, key, ITensor)
            push!(Ts, T)
            push!(idx, parse_numbers(key, delimiter))
        end
    end
    return Ts[sortperm(idx)]
end

function load_structs(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)

    # Read in the states 
    d["ground_state"] = read(f, "ground_state", MPS)
    d["phi_t"] = read(f, "phi_t", MPS)
    d["psi_t"] = read(f, "psi_t", MPS)
    d["mpo"] = read(f, "mpo", MPO)

    # Read in the ITensor vectors 
    if d["DMRG_LBO"]
        optimized_basis = read_ITensor_vector(d, f, "optimized_basis_", "_")
    else
        optimized_basis = nothing
    end

    close(f)

    # Make the structs
    sites = siteinds(d["ground_state"])
    params = Parameters(d["N"], d["t"], d["U"], d["ω"], d["g0"], d["g1"],d["λ"],
                        d["doping"],d["max_phonons"], d["init_phonons"],
                        d["DMRG_numsweeps"], d["DMRG_noise"],d["DMRG_maxdim"],
                        d["DMRG_cutoff"], d["DMRG_LBO"], d["max_LBO_dim"],d["min_LBO_dim"],
                        d["mid"], d["T"], d["τ"], 
                        d["TEBD_cutoff"], d["TEBD_maxdim"], d["TEBD_LBO"])
    hubbholst = HubbardHolsteinModel(params, sites)
    dmrg_results = DMRGResults(d["ground_state"], d["ground_state_energy"], 
                            d["ground_state_entropy"], optimized_basis)
    equilibrium_corr = EquilibriumCorrelations(d["spin"], d["charge"], 
                                                d["sSC"], d["pSC"], d["dSC"])
    tebd_results = TEBDResults(d["entropy"], d["self_overlap"],
                                d["corrs"], d["phi_t"], d["psi_t"])

    return params, hubbholst, dmrg_results, equilibrium_corr, tebd_results
end