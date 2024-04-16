
module Surface3D

using PlotKit

using ..Axes3D: Box3, Box3f, AxisOptions3, Axis3, drawaxis3, get_hexagon
using ..CoreRT: Grid, transform, addlimits, hitshape
using ..ShapeSurface: Surface
using ..RTloop: Lighting
using ..RayTracer: Camera, raytrace_main
using ..RTcairo: cairoimagefrommatrix
using ..AxisDrawables3D: AxisDrawable3

export SurfaceOptions, parse_raytrace_options, Raytrace, surf, drawraytrace

#
# A Surface is a subtype of Shape, something that can be raytraced.
#

const pk = PlotKit


Base.@kwdef mutable struct SurfaceOptions
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    zmin = missing
    zmax = missing
    sampleheight = true
    samplegradient = true
    texture1 = Grid(:red)
    texture2 = Grid(:green)
    #rto = RaytraceOptions()
    # below was in RaytraceOptions
    lighting = Lighting()
    refine = true
    shadows = false
    axisoptions3 = AxisOptions3()
    renderwidth = missing
    renderheight = missing
end


# why is Raytrace separate from RaytraceOptions?
# We should just unify these.
Base.@kwdef mutable struct Raytrace
    box
    shapes
    axis3 = nothing
    # below was in RaytraceOptions
    lighting = Lighting()
    refine = true
    shadows = false
    axisoptions3 = AxisOptions3()
    renderwidth = missing
    renderheight = missing
end


function PlotKit.draw(ad::AxisDrawable3, rt::Raytrace)
    #drawraytrace(ctx, axis3, shapes, box, rto::RaytraceOptions)
    renderwidth = ifnotmissing(rt.renderwidth, rt.axis3.width)
    renderheight = ifnotmissing(rt.renderheight, rt.axis3.height)
    camera = Camera(rt.axis3.axismap3, rt.axis3.width,
                    rt.axis3.height, renderwidth, renderheight)

    shapes_in_cube_coords = [transform(s, rt.axis3) for s in rt.shapes]
    for s in shapes_in_cube_coords
        addlimits(s, rt.axis3.axismap3.cube)
    end
     

    X = raytrace_main(camera, rt.lighting,
                      shapes_in_cube_coords,
                      rt.refine, rt.shadows)
    
    image = cairoimagefrommatrix(X)
    drawimage(ad.ctx, camera, image)
end


#function raytrace(shapes, box, rto::RaytraceOptions)
#    axis3 = Axis3(box, rto.axisoptions3)
#    d = Drawable(axis3.width, axis3.height)
#    rt = Raytrace(box, shapes, axis3, rto)
#    return rt
#end
# 
# 
# function extra_raytrace_options!(rto; kw...)
#     setoptions!(rto.axisoptions3, "axisoptions3_", kw...)
#     setoptions!(rto.axisoptions3.tickbox, "tickbox_", kw...)
#     setoptions!(rto.axisoptions3.axisbox, "axisbox_", kw...)
#     setoptions!(rto.axisoptions3.ticks, "ticks_", kw...)
#     setoptions!(rto.axisoptions3.axisstyle3, "axisstyle3_", kw...)
# end
#     
#     
# function parse_raytrace_options(; kw...)
#     rto = RaytraceOptions()
#     setoptions!(rto, "", kw...)
#     extra_raytrace_options!(rto; kw...)
#     return rto
# end
# 

function Raytrace(shapes, box; kw...)
    rt = Raytrace(; shapes, box, allowed_kws(Raytrace, kw)...)
    setoptions!(rt.axisoptions3, "axisoptions3_", kw...)
    setoptions!(rt.axisoptions3.tickbox, "tickbox_", kw...)
    setoptions!(rt.axisoptions3.axisbox, "axisbox_", kw...)
    setoptions!(rt.axisoptions3.ticks, "ticks_", kw...)
    setoptions!(rt.axisoptions3.axisstyle3, "axisstyle3_", kw...)
    rt.axis3 = Axis3(box, rt.axisoptions3)
#    rt = raytrace(shapes, box, rto)
    #d = Drawable(axis3.width, axis3.height)
    #rt = Raytrace(box, shapes, axis3, rto)
    
    return rt
end


#function PlotKit.draw(rt::Raytrace)
#    d = Drawable(rt.axis3.width, rt.axis3.height)
#    drawaxis3(d.ctx, rt.axis3)
#    drawraytrace(d.ctx, rt.axis3, rt.shapes, rt.box, rt.rto)
#    return d
#end

function PlotKit.draw(rt::Raytrace)
    ad = AxisDrawable3(rt.axis3)
    drawaxis3(ad)
    #drawraytrace(ad.ctx, rt.axis3, rt.shapes, rt.box, rt.rto)
    draw(ad, rt)
    return ad
end

    
##############################################################################
# surface


function getlimitsfromfn(box2, zfun::Function)
    x = range(box2.xmin, box2.xmax, length = 100)
    y = range(box2.ymin, box2.ymax, length = 100)
    zmin = Inf
    zmax = -Inf
    for xi in 1:length(x)
        for yi in 1:length(y)
            z = zfun(x[xi],y[yi])
            if z < zmin
                zmin = z
            end
            if z > zmax
                zmax = z
            end
        end
    end
    return zmin, zmax
end


#
# SurfaceOptions isn't public. But we have it so that
# we don't need to define a function with a long list of keyword options.
# Instead, we keep the defaults in the SurfaceOptions struct.
#
function surf(zfun, dzfun; kw...)
    so = SurfaceOptions(; allowed_kws(SurfaceOptions, kw)...)
    setoptions!(so.axisoptions3, "axisoptions3_", kw...)
    setoptions!(so.axisoptions3.tickbox, "tickbox_", kw...)
    setoptions!(so.axisoptions3.axisbox, "axisbox_", kw...)
    setoptions!(so.axisoptions3.ticks, "ticks_", kw...)
    setoptions!(so.axisoptions3.axisstyle3, "axisstyle3_", kw...)

    box2 = getbox(so)
    if ismissing(so.zmin) || ismissing(so.zmax)
        computed_zmin, computed_zmax = getlimitsfromfn(box2, zfun)
        so.zmin = ifnotmissing(so.zmin, computed_zmin)
        so.zmax = ifnotmissing(so.zmax, computed_zmax)
    end
    box3 = Box3(so.xmin, so.xmax, so.ymin, so.ymax, so.zmin, so.zmax)
    surface = Surface(zfun, dzfun; so.sampleheight, so.samplegradient,
                      so.texture1, so.texture2)
    Raytrace([surface], box3;
             lighting = so.lighting,
             refine = so.refine,
             shadows = so.shadows,
             axisoptions3 = so.axisoptions3,
             renderwidth = so.renderwidth,
             renderheight= so.renderheight)
end

end
