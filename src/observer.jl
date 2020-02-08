export TEvoObserver,
       NoTEvoObserver

abstract type TEvoObserver end

observe!(o::TEvoObserver,args...; kwargs...) = nothing
checkdone!(o::TEvoObserver,args...; kwargs...) = false
measurement_dt(o::TEvoObserver) = 0

struct NoTEvoObserver <: TEvoObserver
end

const Measurement = Vector{Vector{Float64}}

struct MeasurementObserver <: TEvoObserver
  ops::Vector{String}
  sites::Vector{<: Index}
  measurements::Dict{String, Measurement}
  ts::Vector{Float64}
  dt_measure::Float64
end

measurement_ts(o::MeasurementObserver) = o.ts
measurements(o::MeasurementObserver) = o.measurements
measurement_dt(o::MeasurementObserver) = o.dt_measure

function measure_localops!(obs::MeasurementObserver,
                          wf::ITensor,
                          i::Int)
  for o in ops(obs)
    m = dot(wf, noprime(op(sites(obs),o,i)*wf))
    imag(m)>1e-8 && (@warn "encountered finite imaginary part when measuring $o")
    measurements(obs)[o][end][i]=real(m)
  end
end

function observe!(o::MeasurementObserver, psi; t, kwargs...)
    if t-measurement_ts(o)[end]==measurement_dt(o)
    end
end

struct ErrorObserver <: TEvoObserver
end

