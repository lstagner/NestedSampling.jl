# PolyChord Specific Settings
struct PolyChord
    loglike::Function
    prior::Function
    ndims::Cint
    nderived::Cint
    num_repeats::Cint
    do_clustering::Bool
    precision_crit::Cdouble
    base_dir::String
    file_root::String
    write_resume::Bool
    read_resume::Bool
    write_live::Bool
    write_dead::Bool
    write_stats::Bool
    equals::Bool
    posteriors::Bool
    cluster_posteriors::Bool
    update_files::Cint
    boost_posterior::Cdouble
end

function PolyChord(loglike::Function, nparams::Int, name::String;
                   prior::Function = identity,
                   multimodal::Bool = true,
                   nderived::Int = 0,
                   update_every::Int = 1000,
                   write_output_files::Bool = true,
                   work_dir::String = ".")


end

function polychord_loglike(cube_::Ptr{Cdouble}, ndims::Cint,
                           pars_::Ptr{Cdouble}, nderived::Cint)
    cube = unsafe_wrap(Array, cube_, ndims)
    return lnew
end

function run(S::PolyChord; nlive::Int=25*S.ndims, ztol=0.1, resume=true, verbosity=0)

    loglike = cfunction(multinest_loglike, Cdouble, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint))

    root_dir = dirname(S.root)
    println(root_dir)
    isdir(root_dir) || mkdir(root_dir)

    ccall(multinest_ptr, Void, ( Ref{Cint}, Ref{Cint}, Ref{Cint},
          Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}, Ref{Cint},
          Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ptr{UInt8},
          Ref{Cint}, Ptr{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},
          Ref{Cdouble}, Ref{Cint}, Ptr{Void}, Ptr{Void}, Any ),
          S.ins, S.mmodal, S.ceff, nlive, ztol, S.efr,
          S.ndims, S.npar, S.nclspar, S.maxmodes, S.updint, S.mode_ztol,
          S.root, S.seed, S.pwrap, verbosity, resume, S.outfile, S.initmpi,
          S.logzero, S.maxiter, loglike, dumper, S)
end

