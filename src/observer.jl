export TEvoObserver,
       NoTEvoObserver

abstract type TEvoObserver end

observe!(o::TEvoObserver,args...; kwargs...) = nothing
checkdone!(o::TEvoObserver,args...; kwargs...) = false
measurement_step(o::TEvoObserver) = 0

struct NoTEvoObserver <: TEvoObserver
end

