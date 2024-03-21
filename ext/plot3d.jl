
module Plot3D

using PlotKit

#export mesh, mesh_height_color_fn, mesh_height_fn

export Mesh
using Cairo

using ..Basic3D: Vec3, Array23, Array32, Array33, dot
using ..Axes3D: Axis3, Box3, Box3f, drawaxis3, ctxfromaxis



##############################################################################


Base.@kwdef mutable struct Mesh
    axis3
    x
    y
    Z
    cfun
end

##############################################################################



# cfun is a function, mapping (x,y) to color of patch
function PlotKit.mesh(ctx, ax3, x, y, Z, cfun;
              linecolor = Color(:black),
              linewidth = 0.3,
              xlineindices = 1:size(Z,1),
              ylineindices = 1:size(Z,2), kwargs...
              )


    am3 = ax3.axismap3
    f(q)  = ctxfromaxis(am3, q)

    pt(i,j) =  Vec3(x[i],y[j],Z[i,j])

    # top left index i,j
    patch(i,j) =  [pt(i,j), pt(i,j+1), pt(i+1,j+1), pt(i+1,j), i, j]

    patches = []
    for i=1:length(x)-1
        for j=1:length(y)-1
            push!(patches, patch(i,j))
        end
    end

    function sfun(p)
        c = (p[1] + p[2] + p[3] + p[4])/4
        return -dot(am3.depth, c)
    end
    sort!(patches, by = sfun)

    function pline(p1, p2, col)
        move_to(ctx, f(p1))
        line_to(ctx, f(p2))
        Cairo.set_line_width(ctx, linewidth)
        source(ctx, col)        
        Cairo.stroke(ctx)
    end
    
    for p in patches
        i = p[5]
        j = p[6]

        patchcolor = cfun(x[i], y[j])

        Cairo.set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_SQUARE)
        Cairo.set_line_join(ctx, Cairo.CAIRO_LINE_JOIN_BEVEL)
        move_to(ctx, f(p[1]))
        line_to(ctx, f(p[2]))
        line_to(ctx, f(p[3]))
        line_to(ctx, f(p[4]))
        Cairo.close_path(ctx)
        source(ctx, patchcolor)
        Cairo.fill(ctx)
       
        pline(p[1], p[2],  i in xlineindices ? linecolor : patchcolor )
        pline(p[3], p[4], i+1 in xlineindices ? linecolor : patchcolor )
        pline(p[4], p[1], j in ylineindices ? linecolor : patchcolor )
        pline(p[2], p[3], j+1 in ylineindices ? linecolor : patchcolor )
    end



end




boundingbox3(x,y,z) = Box3(minimum(x), maximum(x),
                           minimum(y), maximum(y),
                           minimum(z), maximum(z))


function PlotKit.draw(me::Mesh; kw...)
    d = Drawable(me.axis3.width, me.axis3.height)
#    drawaxis(d.ctx, me.axis3.ax2)
    drawaxis3(d.ctx, me.axis3)
    mesh(d.ctx, me.axis3, me.x, me.y, me.Z, me.cfun; kw...)
    return d
end

function PlotKit.mesh(x, y, Z, cfun; kwargs...)
    b = boundingbox3(x, y, Z)
    axis3 = Axis3(b; kwargs...)
    me = Mesh(axis3, x, y, Z, cfun)
    return me
end


# no color function
function extra_mesh(x, y, Z; patchcolor = [0.8,0.8,1], kwargs...)
    cfun(i,j) = patchcolor
    mesh(x, y, Z, cfun; kwargs...)
end

function PlotKit.sample_mesh(x, y, f)
    Z = zeros(length(x), length(y))
    for xi in 1:length(x)
        for yi in 1:length(y)
            Z[xi,yi] = f(x[xi],y[yi])
        end
    end
    return Z
end


# accepts height and color functions
function PlotKit.mesh_height_color_fn(x, y, f::Function, cfun::Function; kwargs...)
    Z = sample_mesh(x, y, f)
    mesh(x, y, Z, cfun; kwargs...)
end

# accepts a function instead of a Z matrix
function PlotKit.mesh_height_fn(x, y, f::Function; patchcolor = Color(0.8,0.8,1), kwargs...)
    cfun(i,j) = patchcolor
    mesh_height_color_fn(x, y, f, cfun; kwargs...)
end

##############################################################################    





end


