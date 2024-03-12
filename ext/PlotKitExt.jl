
module PlotKitExt

using PlotKit



# define stuff here
PlotKit.setup3d() = println("working")


include("basic3d.jl")
using Basic3D

include("axes3d.jl")
using Axes3D

include("plot3d.jl")
using Plot3D

include("ctools.jl")
using Ctools

include("corert.jl")
using .CoreRT

include("rtloop.jl")
using .RTloop

include("raytracer.jl")
using .Raytracer

include("surface.jl")
using .Surface

include("rtcairo.jl")
using .RTcairo


##############################################################################
function reexport(m)
    for a in names(m)
        eval(Expr(:export, a))
    end
end


reexport(Basic3D)
reexport(Axes3)
reexport(Plot3)
reexport(Surface)
reexport(RTcairo)




end



