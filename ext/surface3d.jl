
module Surface3D

using PlotKit

using ..Axes3D: Box3, Box3f, AxisOptions3, Axis3, drawaxis3, get_hexagon
using ..CoreRT: Grid, transform, addlimits, hitshape
using ..ShapeSurface: Surface
using ..RTloop: Lighting
using ..RayTracer: Camera, raytrace_main
using ..RTcairo: cairoimagefrommatrix
using ..AxisDrawables3D: AxisDrawable3

export Surface, SurfaceOptions, parse_raytrace_options, raytrace, surf, drawraytrace

const pk = PlotKit

Base.@kwdef mutable struct RaytraceOptions
    lighting = Lighting()
    refine = true
    shadows = false
    axisoptions3 = AxisOptions3()
    renderwidth = missing
    renderheight = missing
end

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
    rto = RaytraceOptions()
end


Base.@kwdef mutable struct Raytrace
    box
    shapes
    axis3
    rto::RaytraceOptions
end


function drawraytrace(ctx, axis3, shapes, box, rto::RaytraceOptions)
    renderwidth = ifnotmissing(rto.renderwidth, axis3.width)
    renderheight = ifnotmissing(rto.renderheight, axis3.height)
    camera = Camera(axis3.axismap3, axis3.width, axis3.height, renderwidth, renderheight)

    shapes_in_cube_coords = [transform(s, axis3) for s in shapes]
    for s in shapes_in_cube_coords
        addlimits(s, axis3.axismap3.cube)
    end
     

    X = raytrace_main(camera, rto.lighting,
                      shapes_in_cube_coords,
                      rto.refine, rto.shadows)
    
    image = cairoimagefrommatrix(X)
    drawimage(ctx, camera, image)
end


function raytrace(shapes, box, rto::RaytraceOptions)
    axis3 = Axis3(box, rto.axisoptions3)
    d = Drawable(axis3.width, axis3.height)
    rt = Raytrace(box, shapes, axis3, rto)
    return rt
end

function extra_raytrace_options!(rto; kw...)
    setoptions!(rto.axisoptions3, "axisoptions3_", kw...)
    setoptions!(rto.axisoptions3.tickbox, "tickbox_", kw...)
    setoptions!(rto.axisoptions3.axisbox, "axisbox_", kw...)
    setoptions!(rto.axisoptions3.ticks, "ticks_", kw...)
    setoptions!(rto.axisoptions3.axisstyle3, "axisstyle3_", kw...)
end
    
    
function parse_raytrace_options(; kw...)
    rto = RaytraceOptions()
    setoptions!(rto, "", kw...)
    extra_raytrace_options!(rto; kw...)
    return rto
end


function raytrace(shapes, box; kw...)
    rto = parse_raytrace_options(; kw...)
    axis3 = Axis3(box, rto.axisoptions3)
    rt = raytrace(shapes, box, rto)
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
    drawraytrace(ad.ctx, rt.axis3, rt.shapes, rt.box, rt.rto)
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


function surf(zfun, dzfun; kw...)
    so = SurfaceOptions()
    setoptions!(so, "", kw...)
    setoptions!(so.rto, "raytrace_", kw...)
    extra_raytrace_options!(so.rto; kw...)

    box2 = getbox(so)
    if ismissing(so.zmin) || ismissing(so.zmax)
        computed_zmin, computed_zmax = getlimitsfromfn(box2, zfun)
        so.zmin = ifnotmissing(so.zmin, computed_zmin)
        so.zmax = ifnotmissing(so.zmax, computed_zmax)
    end
    box3 = Box3(so.xmin, so.xmax, so.ymin, so.ymax, so.zmin, so.zmax)
    surface = Surface(zfun, dzfun; so.sampleheight, so.samplegradient,
                      so.texture1, so.texture2)
    raytrace([surface], box3, so.rto)
end

end
