

module DecTree

using PlotKit
using PlotKit3D
using LinearAlgebra

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)
eye(n) = Matrix(I(n))

function bar(xmin,xmax,ymin,ymax,zmin,zmax)
    red     = Material(0.8*Color(:red), 0.4, 0.45, 0.35, 8, 0.0, 0.0)
    blue    = Material(0.5*Color(:blue), 0.4, 0.45, 0.35, 8, 0.0, 0.0)
    black   = Material(Color(:white), 0.0, 0.0, 0.0, 0, 0.0, 0.0)
    green   = Material(Color(:green), 0.2, 0.45, 0.35, 8, 0.0, 0.0)

    reduniform     = Uniform(red)
    greenuniform   = Uniform(green)
    redgrid        = Grid(red, black,0.01,0.1,0.1)

    A = [1.0  0.0  0.0;
        -1.0  0.0  0.0;
         0.0  1.0  0.0;
         0.0 -1.0  0.0;
         0.0  0.0  1.0;
         0.0  0.0 -1.0;]
    b = [xmax,-xmin,ymax,-ymin,zmax,-zmin]

    textures = [greenuniform, greenuniform, greenuniform, greenuniform, redgrid, greenuniform]
    
    poly = Polytope(A, b, textures)
    return poly
end

        

# main
function main()

    light1 = Light(Vec3(-30,0,50), Color(:white))
    light2 = Light(Vec3(-30,-20,-30), 0.5*Color(:white))
    lighting = Lighting([light1,light2], 0.8*Color(:white))

    shapes = Shape[]
    push!(shapes, bar(0,0.2, 0,0.8,  -3, 0.25))
    push!(shapes, bar(0.2,1, 0,0.8,  -3, 0.5))
    push!(shapes, bar(0,1, 0.8,1,  -3, 0.6))
          
    box = Box3(0,1,0,1,0,1)
    rt = Raytrace(shapes, box; lighting, axisoptions3_fontsize=30, axisoptions3_width=800,
                  axisoptions3_height=800)
    d = draw(rt)
    save(d, plotpath("tree10.pdf"))

end


end

