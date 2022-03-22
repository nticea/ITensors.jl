using ITensors.HDF5

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames 
        n = Name(fn)
        d = getfield(struc, fn)
        try 
            h5open(path, "r+") do file
                write(file, n, d)
            end
        catch
            h5open(path, "w") do file
                write(file, n, d) 
            end
        end
    end
end

function load_structs(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    d["ground_state"] = read(f, "ground_state", MPS)
    d["phi_t"] = read(f, "phi_t", MPS)
    d["psi_t"] = read(f, "psi_t", MPS)
    close(f)

    # Make the structs
    params = Parameters(d["N"], d["t"], d["U"], d["ω"], d["g0"], d["g1"], d["doping"],
                        d["DMRG_numsweeps"], d["DMRG_noise"],d["DMRG_maxdim"],
                        d["DMRG_cutoff"], d["mid"], d["T"], d["τ"], 
                        d["TEBD_cutoff"], d["TEBD_maxdim"])
    hubbholst = HubbardHolsteinModel(params)
    dmrg_results = DMRGResults(d["ground_state"], d["ground_state_energy"], 
                            d["ground_state_entropy"])
    equilibrium_corr = EquilibriumCorrelations(d["spin"], d["charge"], 
                                                d["sSC"], d["pSC"], d["dSC"])
    tebd_results = TEBDResults(d["entropy"], d["self_overlap"], d["phonon_flux"],
                                d["corrs"], d["phi_t"], d["psi_t"])

    return params, hubbholst, dmrg_results, equilibrium_corr, tebd_results
end