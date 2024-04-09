
module PlotKit3D


include("basic3d.jl")
using .Basic3D

include("axes3d.jl")
using .Axes3D

include("axisdrawables3d.jl")
using .AxisDrawables3D

include("plot3d.jl")
using .Plot3D

include("corert.jl")
using .CoreRT

include("shapes.jl")
using .Shapes

include("shapesurface.jl")
using .ShapeSurface

include("shapetransform.jl")
using .ShapeTransform

include("shapeintersection.jl")
using .ShapeIntersection

include("rtloop.jl")
using .RTloop

include("raytracer.jl")
using .RayTracer

include("rtcairo.jl")
using .RTcairo

include("surface3d.jl")
using .Surface3D



##############################################################################
function reexport(m)
    for a in names(m)
        eval(Expr(:export, a))
    end
end


reexport(Basic3D)
reexport(Axes3D)
reexport(AxisDrawables3D)
reexport(Plot3D)
reexport(CoreRT)
reexport(Shapes)
reexport(ShapeSurface)
reexport(ShapeTransform)
reexport(ShapeIntersection)
reexport(RTloop)
reexport(RayTracer)
reexport(Surface3D)
reexport(RTcairo)



end



