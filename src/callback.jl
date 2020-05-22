export TEvoCallback,
       NoTEvoCallback

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

apply!(cb::TEvoCallback,args...; kwargs...) = nothing
checkdone!(cb::TEvoCallback,args...; kwargs...) = false
callback_dt(cb::TEvoCallback) = 0

struct NoTEvoCallback <: TEvoCallback
end

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
measurements(cb::LocalMeasurementCallback) = cb.measurements
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

function apply!(cb::LocalMeasurementCallback, psi; t, sweepend, sweepdir, bond, kwargs...)
    # perform measurements only at the end of a sweep and at measurement steps
    prev_t = length(measurement_ts(cb))>0 ? measurement_ts(cb)[end] : 0
    if (t-prev_t==callback_dt(cb) || t==prev_t) && sweepend && sweepdir == "left"
        if t != prev_t
            push!(measurement_ts(cb), t)
            foreach(x->push!(x,zeros(length(psi))), values(measurements(cb)) )
        end
        measure_localops!(cb,psi[bond]*psi[bond+1],bond+1)
        bond==1 && measure_localops!(cb,psi[bond]*psi[bond+1],bond)
    end
end

struct ErrorCallback <: TEvoCallback
end

