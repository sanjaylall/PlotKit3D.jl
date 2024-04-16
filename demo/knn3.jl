
module Knn3

using PlotKit
using PlotKit3D
using LinearAlgebra
import Zygote
import Random

plotpath(x) = joinpath(ENV["HOME"], "plots/", plotname(x))

# main 5,11,15, 29, 36

filename() = "softknn"

eye(n) = Matrix(I(n))

# needs a function filename()
function plotname(x)
    a=findlast(isequal('.'), x)
    basename = x[1:a-1]
    extension = x[a+1:end]
    basename2 = replace(basename, '.' => "__")
    fname = filename()
    return  "$(fname)_$(basename2).$(extension)"
end

function get_surf(rho, pts, h)
    normsq(x) = sum(x.*x)
    norm(x) = sqrt(normsq(x))
    n = size(pts,1)

    function weight(x,xi)
        num = exp(-normsq(x-xi)/(rho*rho))
        den = 0.0
        for i=1:n
            den = den + exp(-normsq(x-pts[i,:])/(rho*rho)) 
        end
        return num/den
    end
    function zfn(x,y)
        yh = 0.0
        for i=1:n
            yh = yh + weight([x,y],pts[i,:])*h[i]
        end
        return yh
    end
    dzfn(x,y) = Zygote.gradient(zfn, x,y)

    return zfn, dzfn
end

function get_shapes(rho, pts, h)

    n = size(pts,1)

    
    red     = Material(0.8*Color(:red), 0.4, 0.45, 0.35, 8, 0.0, 0.0)
    blue    = Material(0.5*Color(:blue), 0.4, 0.45, 0.35, 8, 0.0, 0.0)
    black   = Material(Color(:white), 0.0, 0.0, 0.0, 0, 0.0, 0.0)
    green   = Material(Color(:green), 0.2, 0.45, 0.35, 8, 0.0, 0.0)
    
    toptexture    = Grid(red, black)
    bottomtexture = Grid(green, black)
    blueuniform   = Uniform(blue)

    sampleheight = true
    samplegradient = true

    zfn, dzfn = get_surf(rho, pts, h)

    surface = Surface(zfn, dzfn, toptexture, bottomtexture, sampleheight, samplegradient)

    shapes = Array{Shape,1}(undef,0)
    push!(shapes, surface)
    for i=1:n
        r = 0.125
        ell = Ellipsoid(eye(3)/(r*r),[pts[i,1],pts[i,2], h[i]], blueuniform)
        push!(shapes, ell)
    end

    return shapes
    
end


function main(n=36,r=false)
    xlimits = [-3,3]
    ylimits = [-3,3]

    pts = [ -2 2.5
          -2 0
          0 -2
          2.5 0
          1 2
          2 2]
    
    n = size(pts,1)
    h = [-2.0
         -1.2
         -2
         -0.6
         0.2
         1.0]

    rho = 1

    box = Box3(-3,3,-3,3,-3,3)

    light1 = Light(Vec3(-30,0,20), 0.8*Color(:white))
    light2 = Light(Vec3(-30,-20,-30), 0.5*Color(:white))
    lighting = Lighting([light1,light2], 0.8*Color(:white))


    for rho in [2, 1, 0.5]
        shapes = get_shapes(rho, pts, h)
        rt = Raytrace(shapes, box; axisoptions3_fontsize=30, axisoptions3_width=400,
                  axisoptions3_height=400)
        d = draw(rt)
        save(d, plotpath("knn3d_soft_$(rho).pdf"))
    end
    


    
end


end

