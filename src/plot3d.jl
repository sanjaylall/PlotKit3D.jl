
module Plot3D

using PlotKit

#export mesh, 

export Mesh, mesh, mesh_height_color_fn, mesh_height_fn, sample_mesh
using Cairo

using ..Basic3D: Vec3, Array23, Array32, Array33, dot
using ..Axes3D: Axis3, Box3, Box3f, drawaxis3, ctxfromaxis
using ..AxisDrawables3D: AxisDrawable3


##############################################################################


Base.@kwdef mutable struct Mesh
    axis3
    x
    y
    Z
    patchcolor = Color(0.8, 0.8, 1.0)
    linecolor = Color(:black)
    linewidth = 0.3
    xlineindices = nothing
    ylineindices = nothing
end

##############################################################################

ifnotnothing(x::Nothing, y) = y
ifnotnothing(x, y) = x

atij(i, j, f::Function) = f(i,j)
atij(i, j, f) = f

# cfun is a function, mapping (x,y) to color of patch
function drawmesh(ad::AxisDrawable3, me::Mesh)
    x = me.x
    y = me.y
    Z = me.Z

    xlineindices = ifnotnothing(me.xlineindices, 1:size(Z,1))
    ylineindices = ifnotnothing(me.ylineindices, 1:size(Z,2))

    depth = ad.axis3.axismap3.depth

    pt(i,j) =  Vec3(x[i],y[j],Z[i,j])

    # top left index i,j
    patch(i,j) =  [pt(i,j), pt(i,j+1), pt(i+1,j+1), pt(i+1,j), i, j]

    patches = [patch(i,j) for i=1:length(x)-1 for j=1:length(y)-1]
    
    sfun(p) = -dot(depth,  (p[1] + p[2] + p[3] + p[4])/4 )

    sort!(patches, by = sfun)

    for p in patches
        i = p[5]
        j = p[6]

        patchcolor = atij(x[i], y[j], me.patchcolor)

        Cairo.set_line_cap(ad.ctx, Cairo.CAIRO_LINE_CAP_SQUARE)
        Cairo.set_line_join(ad.ctx, Cairo.CAIRO_LINE_JOIN_BEVEL)

        line(ad, [p[1], p[2], p[3], p[4]]; fillcolor = patchcolor, closed = true)
       
        line(ad, p[1], p[2], linestyle = LineStyle(i in xlineindices ? me.linecolor : patchcolor, me.linewidth))
        line(ad, p[3], p[4], linestyle = LineStyle(i+1 in xlineindices ? me.linecolor : patchcolor , me.linewidth))
        line(ad, p[4], p[1], linestyle = LineStyle(j in ylineindices ? me.linecolor : patchcolor, me.linewidth))
        line(ad, p[2], p[3], linestyle = LineStyle(j+1 in ylineindices ? me.linecolor : patchcolor, me.linewidth))
    end



end

boundingbox3(x,y,z) = Box3(minimum(x), maximum(x),
                           minimum(y), maximum(y),
                           minimum(z), maximum(z))

function PlotKit.draw(me::Mesh; kw...)
    ad = AxisDrawable3(me.axis3)
    drawaxis3(ad)
    drawmesh(ad, me)
    return ad
end

# f may be a function or a matrix
function Mesh(x, y, f; axis3 = nothing, kw...)
    Z = sample_mesh(x, y, f)
    b = boundingbox3(x, y, Z)
    axis3 = ifnotnothing(axis3, Axis3(b; kw...))
    me = Mesh(; axis3, x, y, Z,  allowed_kws(Mesh, kw)...)
    return me
end

sample_mesh(x, y, Z) = Z

function sample_mesh(x, y, f::Function)
    Z = zeros(length(x), length(y))
    for xi in 1:length(x)
        for yi in 1:length(y)
            Z[xi,yi] = f(x[xi],y[yi])
        end
    end
    return Z
end


# accepts height and color functions
#function mesh_height_color_fn(x, y, f::Function, cfun::Function; kwargs...)
#    Z = sample_mesh(x, y, f)
#    Mesh(x, y, Z; cfun, kwargs...)
#end

# accepts a function instead of a Z matrix
#function mesh_height_fn(x, y, f::Function; patchcolor = Color(0.8,0.8,1), kwargs...)
#    cfun(i,j) = patchcolor
#    mesh_height_color_fn(x, y, f, cfun; kwargs...)
#end

##############################################################################    





end


