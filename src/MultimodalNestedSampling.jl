module MultimodalNestedSampling

import Base: run

export MultiNest, run, multinest_dumper, multinest_loglike, multinest_default_dumper

lib = Libdl.find_library("libmultinest",String[])
println(lib)
dl = Libdl.dlopen_e(lib)
dl == C_NULL && error("Cannot find MultiNest Library")
const multinest_ptr = Libdl.dlsym(dl,:__nested_MOD_nestrun)
multinest_ptr == C_NULL && error("Cannot link to Multinest Library")

include("multinest.jl")

# package code goes here

end # module
