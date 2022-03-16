using JLD

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames 
        n = Name(fn)
        d = getfield(struc, fn)
        jldopen(path, "w") do file
            write(file, n, d) 
        end
    end
end