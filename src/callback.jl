export TEvoCallback,
    NoTEvoCallback,
    LocalMeasurementCallback,
    SpecCallback,
    measurement_ts

"""
A TEvoCallback can implement the following methods:

- apply!(cb::TEvoCallback, psi ; t, kwargs...): apply the callback with the
current state `psi` (e.g. perform some measurement)

- checkdone!(cb::TEvoCallback, psi; t, kwargs...): check whether some criterion
to stop the time evolution is satisfied (e.g. convergence of cbervable, error
too large) and return `true` if so.

- callback_dt(cb::TEvoCallback): time-steps at which the callback needs access
for the wave-function (e.g. for measurements). This is used for TEBD evolution
where several time-steps are bunched together to reduce the cost.
"""
abstract type TEvoCallback end

struct NoTEvoCallback <: TEvoCallback
end

apply!(cb::NoTEvoCallback,args...; kwargs...) = nothing
checkdone!(cb::NoTEvoCallback,args...; kwargs...) = false
callback_dt(cb::NoTEvoCallback) = 0


const Measurement = Vector{Vector{Float64}}

struct LocalMeasurementCallback <: TEvoCallback
    ops::Vector{String}
    sites::Vector{<: Index}
    measurements::Dict{String, Measurement}
    ts::Vector{Float64}
    dt_measure::Float64
end

function LocalMeasurementCallback(ops,sites,dt_measure)
    return LocalMeasurementCallback(ops,
                                    sites,
                                    Dict(o => Measurement[] for o in ops),
                                    Vector{Float64}(),
                                    dt_measure)
end

measurement_ts(cb::LocalMeasurementCallback) = cb.ts
ITensors.measurements(cb::LocalMeasurementCallback) = cb.measurements
callback_dt(cb::LocalMeasurementCallback) = cb.dt_measure
ops(cb::LocalMeasurementCallback) = cb.ops
sites(cb::LocalMeasurementCallback) = cb.sites

function Base.show(io::IO, cb::LocalMeasurementCallback)
    println(io, "LocalMeasurementCallback")
    println(io, "Operators: ", ops(cb))
    if length(measurement_ts(cb))>0
        println(io, "Measured times: ", callback_dt(cb):callback_dt(cb):measurement_ts(cb)[end])
    else
        println(io, "No measurements performed")
    end
end

function measure_localops!(cb::LocalMeasurementCallback,
                          wf::ITensor,
                          i::Int)
    for o in ops(cb)
        m = dot(wf, noprime(op(sites(cb),o,i)*wf))
        imag(m)>1e-8 && (@warn "encountered finite imaginary part when measuring $o")
        measurements(cb)[o][end][i]=real(m)
    end
end

function apply!(cb::LocalMeasurementCallback, psi; t, sweepend, sweepdir, bond, alg, kwargs...)
    prev_t = length(measurement_ts(cb))>0 ? measurement_ts(cb)[end] : 0

    # perform measurements only at the end of a sweep (TEBD: for finishing
    # evolution over the measurement time interval; TDVP: when sweeping left)
    # and at measurement steps.
    # For TEBD algorithms we want to perform measurements only in the final
    # sweep over odd bonds. For TDVP we can perform measurements to the right of
    # each bond when sweeping back left.
    if (t-prev_t≈callback_dt(cb) || t==prev_t) && sweepend && (bond % 2 ==1 || !(alg isa TEBDalg))
        if t != prev_t
            push!(measurement_ts(cb), t)
            foreach(x->push!(x,zeros(length(psi))), values(measurements(cb)) )
        end
        wf = psi[bond]*psi[bond+1]
        measure_localops!(cb,wf,bond+1)
        if alg isa TEBDalg
            measure_localops!(cb,wf,bond)
        elseif bond==1
            measure_localops!(cb,wf,bond)
        end
    end
end


checkdone!(cb::LocalMeasurementCallback,args...) = false

struct SpecCallback <: TEvoCallback
    truncerrs::Vector{Float64}
    current_truncerr::Base.RefValue{Float64}
    entropies::Measurement
    bonddims::Vector{Vector{Int64}}
    bonds::Vector{Int64}
    ts::Vector{Float64}
    dt_measure::Float64
end

function SpecCallback(dt,psi::MPS,bonds::Vector{Int64}=collect(1:length(psi)-1))
    bonds = sort(unique(bonds))
    if maximum(bonds) > length(psi)-1 || minimum(bonds)<1
        throw("bonds must be between 1 and $(length(psi)-1)")
    end
    return SpecCallback(Vector{Float64}(), Ref(0.0), Measurement(),
                        Vector{Vector{Int64}}(),bonds, Vector{Float64}(),dt)
end


measurement_ts(cb::SpecCallback) = cb.ts
ITensors.measurements(cb::SpecCallback) = Dict("entropy"=> cb.entropies,
                                      "bonddim"=>cb.bonddims,
                                      "truncerrs"=>cb.truncerrs)

callback_dt(cb::SpecCallback) = cb.dt_measure
bonds(cb::SpecCallback) = cb.bonds

function Base.show(io::IO, cb::SpecCallback)
    println(io, "SpecCallback")
    if length(measurement_ts(cb))>0
        println(io, "Measured times: ", callback_dt(cb):callback_dt(cb):measurement_ts(cb)[end])
    else
        println(io, "No measurements performed")
    end
end

function apply!(cb::SpecCallback, psi; t, sweepend,bond,spec,sweepdir, kwargs...)
    cb.current_truncerr[] += truncerror(spec)
    prev_t = length(measurement_ts(cb))>0 ? measurement_ts(cb)[end] : 0
    if (t-prev_t≈callback_dt(cb) || t==prev_t) && sweepend
        if t != prev_t
            push!(measurement_ts(cb), t)
            push!(cb.bonddims,zeros(Int64,length(cb.bonds)))
            push!(cb.entropies,zeros(length(cb.bonds)))
        end

        if bond in bonds(cb)
            i = findfirst(x->x==bond,bonds(cb))
            cb.bonddims[end][i] = length(eigs(spec))
            cb.entropies[end][i] = entropy(spec)
        end
        if sweepdir=="right" && bond==length(psi)-1
            push!(cb.truncerrs, cb.current_truncerr[])
        elseif sweepdir=="left" && bond==1
            push!(cb.truncerrs, cb.current_truncerr[])
        end
    end
end

checkdone!(cb::SpecCallback,args...) = false
