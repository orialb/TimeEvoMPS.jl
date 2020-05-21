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

struct MeasurementCallback <: TEvoCallback
    ops::Vector{String}
    sites::Vector{<: Index}
    measurements::Dict{String, Measurement}
    ts::Vector{Float64}
    dt_measure::Float64
end

measurement_ts(cb::MeasurementCallback) = o.ts
measurements(cb::MeasurementCallback) = o.measurements
callback_dt(cb::MeasurementCallback) = o.dt_measure

function measure_localops!(cb::MeasurementCallback,
                          wf::ITensor,
                          i::Int)
    for o in ops(cb)
        m = dot(wf, noprime(op(sites(cb),o,i)*wf))
        imag(m)>1e-8 && (@warn "encountered finite imaginary part when measuring $o")
        measurements(cb)[o][end][i]=real(m)
    end
end

function apply!(cb::MeasurementCallback, psi; t, kwargs...)
    if t-measurement_ts(o)[end]==callback_dt(o)
    end
end

struct ErrorCallback <: TEvoCallback
end

