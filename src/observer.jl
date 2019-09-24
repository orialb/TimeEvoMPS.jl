abstract type TEvoObserver end

observe!(o::TEvoObserver,args...) = nothing
checkdone!(o::TEvoObserver,args...) = false
measurement_step(o::TEvoObserver) = 0

struct NoTEvoObserver <: TEvoObserver
end

