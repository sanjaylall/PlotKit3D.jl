

module Knn2

using PlotKit
using PlotKit3D
using LinearAlgebra

const pk3 = PlotKit3D

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)

eye(n) = Matrix(I(n))

normsq(a) = sum(a.*a)
norm(a) = sqrt(normsq(a))

function voronoik(c1,c2,p,zbot,ztop, tol)
    xmin = -3
    ymin = -3

    n = size(p,1)
    cx1 = c1[1]
    cy1 = c1[2]

    cx2 = c2[1]
    cy2 = c2[2]
    
    function f_in(a)
        x = a.x
        y = a.y
        z = a.z

        if z > ztop || z < zbot
            return false
        end
        for i=1:n
            qx = p[i,1]
            qy = p[i,2]
            if normsqdiff(x, y, qx, qy) < normsqdiff(x, y, cx1, cy1)
                return false
            end
            if normsqdiff(x, y, qx, qy) < normsqdiff(x, y, cx2, cy2)
                return false
            end
        end
        return true
    end

    function f_nrm(a)
        x = a.x
        y = a.y
        z = a.z

        if abs(z-ztop) < tol
            return Vec3(0,0,1),2
        end
        if abs(z-zbot) < tol
            return Vec3(0,0,-1),2
        end
        if abs(x-xmin) < tol
            return Vec3(-1,0,0),1
        end
        if abs(y-ymin) < tol
            return Vec3(0,-1,0),1
        end
        for i=1:n
            qx = p[i,1]
            qy = p[i,2]
            if abs(sqrt(normsqdiff(x, y, qx, qy)) - sqrt(normsqdiff(x, y, cx1, cy1))) < 2*tol
                nx = qx-cx1
                ny = qy-cy1
                return Vec3(nx, ny, 0), 1
            end
            if abs(sqrt(normsqdiff(x, y, qx, qy)) - sqrt(normsqdiff(x, y, cx2, cy2))) < 2*tol
                nx = qx - cx2
                ny = qy - cy2
                return Vec3(nx, ny, 0), 1
            end
        end
        println("failed");
        return Vec3(1,0,0),1
    end
    return f_in, f_nrm
end

normsqdiff(x1,y1, x2,y2) = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)


function main()

    light1 = Light(Vec3(-30,0,50), Color(:white))
    light2 = Light(Vec3(-30,-20,-30), 0.5*Color(:white))
    lighting = Lighting([light1,light2], 0.8*Color(:white))

    red     = Material(0.8*Color(:red), 0.4, 0.45, 0.35, 8, 0.0, 0.0)
    blue    = Material(0.5*Color(:blue), 0.4, 0.45, 0.35, 8, 0.0, 0.0)
    black   = Material(Color(:white), 0.0, 0.0, 0.0, 0, 0.0, 0.0)
    green   = Material(Color(:green), 0.2, 0.45, 0.35, 8, 0.0, 0.0)
    
    reduniform     = Uniform(red)
    blueuniform     = Uniform(blue)
    greenuniform     = Uniform(green)
    
    shapes = Array{Shape,1}(undef,0)

    zbot = -3

    p = [ -2 2.5
          -2 0
          0 -2
          2.5 0
          1 2
          2 2]
    
    n = size(p,1)
    h = [-2.0
         -1.2
         -2
         -0.6
         0.2
         1.0]

    tol1 = 0.001
    tol2 = 3*tol1
    
    for i=1:n
        for j = i+1:n
            others1 = filter(x -> x!= i, 1:n)
            others  = filter(x -> x!= j, others1)

            f_in, f_nrm = voronoik(p[i,:], p[j,:], p[others,:], zbot, (h[i]+h[j])/2, tol2)
            arb = ArbitrarySolid(f_in, f_nrm, tol1, greenuniform, Grid(red,black))
            push!(shapes, arb)
        end
    end

    for i=1:n
        r = 0.125
        ell = Ellipsoid(eye(3)/(r*r),[p[i,1],p[i,2], h[i]], blueuniform)
        push!(shapes, ell)
    end
    
    box = Box3(-3,3,-3,3,-3,3)
    rt = Raytrace(shapes, box; axisoptions3_fontsize=30, axisoptions3_width=400,
                  axisoptions3_height=400)
    d = draw(rt)
    save(d, plotpath("knn3d2.pdf"))
    
end


end
