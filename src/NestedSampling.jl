module NestedSampling

using Libdl

import Base: run

export MultiNest, run, multinest_dumper, multinest_loglike, multinest_default_dumper

function __init__()
    lib = Libdl.find_library("libmultinest",String[])
    println(lib)
    dl = Libdl.dlopen_e(lib)
    dl == C_NULL && error("Cannot find MultiNest Library")
    global const multinest_ptr = Libdl.dlsym(dl,:__nested_MOD_nestrun)
    multinest_ptr == C_NULL && error("Cannot link to Multinest Library")
end

include("multinest.jl")

# package code goes here

end # module
