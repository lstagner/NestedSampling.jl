# Multinest Dumper Routine
function multinest_default_dumper(physlive::Matrix{T},posterior::Matrix{T},
                          paramconstr::Matrix{T}, maxloglike::T,
                          logz::T, inslogz::T, logzerr::T) where T <: AbstractFloat
end

function multinest_dumper(nsamples_::Ptr{Cint},
                          nlive_::Ptr{Cint},
                          npar_::Ptr{Cint},
                          physlive_::Ptr{Ptr{Cdouble}},
                          posterior_::Ptr{Ptr{Cdouble}},
                          paramconstr_::Ptr{Ptr{Cdouble}},
                          maxloglike_::Ptr{Cdouble},
                          logZ_::Ptr{Cdouble},
                          insLogZ_::Ptr{Cdouble},
                          logZerr_::Ptr{Cdouble},
                          nested_::Ptr{Nothing})

    nsamples = unsafe_load(nsamples_)
    nlive = unsafe_load(nlive_)
    npar = unsafe_load(npar_)
    physlive = unsafe_wrap(Array,unsafe_load(physlive_),(nlive, Int32(npar+1)))
    posterior = unsafe_wrap(Array,unsafe_load(posterior_),(nsamples, Int32(npar+2)))
    paramconstr = unsafe_wrap(Array,unsafe_load(paramconstr_),(npar, Int32(4)))
    maxloglike = unsafe_load(maxloglike_)
    logZ = unsafe_load(logZ_)
    insLogZ = unsafe_load(insLogZ_)
    logZerr = unsafe_load(logZerr_)
    nested = unsafe_pointer_to_objref(nested_)::MultiNest
    nested.dumper(physlive, posterior, paramconstr, maxloglike, logZ, insLogZ, logZerr)
    return
end

# Multinest Specific Settings
struct MultiNest
    loglike::Function
    root::String
    ins::Cint
    mmodal::Cint
    ceff::Cint
    efr::Cdouble
    mode_ztol::Cdouble
    ndims::Cint
    npar::Cint
    nclspar::Cint
    updint::Cint
    maxmodes::Cint
    pwrap::Vector{Cint}
    seed::Cint
    outfile::Cint
    initmpi::Cint
    logzero::Cdouble
    maxiter::Cint
    dumper::Function
end

function MultiNest(loglike::Function, nparams::Int, name::String;
                   prior::Function = identity,
                   importance_sampling::Bool = true,
                   multimodal::Bool = true,
                   max_modes::Int = 100,
                   mode_ztol::Float64 = -1e90,
                   nderived::Int = 0,
                   ncluster_dims::Int = 0,
                   mode_npar::Int = (nparams - nderived),
                   constant_efficiency::Bool = false,
                   efficiency_rate::AbstractFloat = 0.8,
                   seed::Int = -1,
                   update_every::Int = 1000,
                   write_output_files::Bool = true,
                   work_dir::String = ".",
                   init_mpi::Bool = true,
                   periodic::Vector{Bool} = fill(false, nparams - nderived),
                   maxiter::Int = 0,
                   dumper::Function = multinest_default_dumper)

    ndims = nparams - nderived

    # Verify Inputs
    nparams > 0 || error("nparams must be greater than zero")
    nparams > nderived || error("nderived must be less than nparams")
    ncluster_dims < ndims || error("ncluster_dims must be less than (nparams - nderived)")
    0 < efficiency_rate <= 1 || error("efficiency_rate must be in the range (0,1]")
    max_modes > 0 || error("max_modes must be greater than 0")
    length(periodic) == ndims || error("size(periodic) != nparams - nderived")
    maxiter >= 0 || error("max_iter must be non-negative")
    update_every > 0 || error("update_every must be non-negative")
    seed >= -1 || error("seed must be positive integer")

    MultiNest(x->loglike(prior(x)), rpad(joinpath(work_dir,name*"_"), 100, ' '), Cint(importance_sampling),
              Cint(multimodal), Cint(constant_efficiency), Cdouble(efficiency_rate),
              Cdouble(mode_ztol), Cint(ndims), Cint(nparams), Cint(ncluster_dims),
              Cint(update_every), Cint(max_modes), Cint.(periodic), Cint(seed),
              Cint(write_output_files), Cint(init_mpi), Cdouble(nextfloat(-Inf)),
              Cint(maxiter), dumper)

end

function multinest_loglike(cube_::Ptr{Cdouble}, ndim_::Ptr{Cint},
                           npar_::Ptr{Cint}, lnew_::Ptr{Cdouble},
                           nested_::Ptr{Nothing})
    ndim = unsafe_load(ndim_)
    npar = unsafe_load(npar_)
    cube = unsafe_wrap(Array, cube_, npar)
    nested = unsafe_pointer_to_objref(nested_)::MultiNest
    lnew = nested.loglike(cube)
    unsafe_store!(lnew_, lnew)
    return
end

function run(S::MultiNest; nlive::Int=25*S.ndims, ztol=0.1, resume=false, verbosity=0)

    loglike = @cfunction(multinest_loglike, Nothing, (Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Nothing}))

    dumper = @cfunction(multinest_dumper, Nothing, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                                         Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                                         Ptr{Ptr{Cdouble}}, Ptr{Cdouble},
                                         Ptr{Cdouble}, Ptr{Cdouble},Ptr{Cdouble},
                                         Ptr{Nothing}))

    root_dir = dirname(S.root)
    println(root_dir)
    isdir(root_dir) || mkdir(root_dir)

    ccall(multinest_ptr, Nothing, ( Ref{Cint}, Ref{Cint}, Ref{Cint},
          Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}, Ref{Cint},
          Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ptr{UInt8},
          Ref{Cint}, Ptr{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},
          Ref{Cdouble}, Ref{Cint}, Ptr{Nothing}, Ptr{Nothing}, Any ),
          S.ins, S.mmodal, S.ceff, nlive, ztol, S.efr,
          S.ndims, S.npar, S.nclspar, S.maxmodes, S.updint, S.mode_ztol,
          S.root, S.seed, S.pwrap, verbosity, resume, S.outfile, S.initmpi,
          S.logzero, S.maxiter, loglike, dumper, S)
end
