
module Knn1

using PlotKit
using PlotKit3D
using LinearAlgebra

const pk3 = PlotKit3D

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)

function eye(n)
    A = zeros(n,n)
    for i=1:n
        A[i,i]=1
    end
    return A
end

function hemisphere()
    function f_in(x,y,z)
        if x*x + y*y + z*z < 1 && y>0
            return true
        end
        return false
    end

    function f_nrm(x,y,z)
        if y < 0.02 && y > -0.02
            return 0,1,0,1
        end
        return x,y,z,1
    end
    return f_in, f_nrm
end

normsq(a) = sum(a.*a)
norm(a) = sqrt(normsq(a))

function voronoi()
    c = [0,0, 0]
    p = [2 1 0
         4 5 1
         -3 -3 -1
         1 0 0
         -1 4  0]
    
    p = p /3
    n = size(p,1)
    
    function f_in(x,y,z)
        a = [x,y,z]
        for i=1:n
            q = p[i,:]
            if norm(a-q) < norm(a-c)
                return false
            end
        end
        return true
    end

    function f_nrm(x,y,z)
        a = [x,y,z]
        z = zeros(n)
        for i=1:n
            q = p[i,:]
            z[i] = abs(norm(a-q) - norm(a-c)) 
        end
        y,j = findmin(z)
        q = p[j,:]
        nrm = q-c
        return nrm[1], nrm[2], nrm[3],1
    end
    return f_in, f_nrm
end

normsqdiff(x1,y1, x2,y2) = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)

function voronoi2d(c,p,zbot,ztop, tol)

    xmin = -3
    ymin = -3
    n = size(p,1)
    cx = c[1]
    cy = c[2]
    
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
            if normsqdiff(x, y, qx, qy) < normsqdiff(x, y, cx, cy)
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
            return Vec3(0,0,1), 2
        end
        if abs(z-zbot) < tol
            return Vec3(0,0,-1), 2
        end
        if abs(x-xmin) < tol
            return Vec3(-1,0,0), 1
        end
        if abs(y-ymin) < tol
            return Vec3(0,-1,0), 1
        end
        
        a = [x,y]
        z = zeros(n)
        for i=1:n
            q = p[i,:]
            z[i] = abs(normsq(a-q) - normsq(a-c)) 
        end
        y,j = findmin(z)
        q = p[j,:]
        nrm = q-c
        return Vec3(nrm[1], nrm[2], 0.0), 1
    end

        
    return f_in, f_nrm
end


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

    #f_in, f_nrm = hemisphere()
    

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
        others = filter(x -> x!= i, 1:n)
        f_in, f_nrm = voronoi2d(p[i,:], p[others,:], zbot, h[i], tol2)
        arb = ArbitrarySolid(f_in, f_nrm, tol1, greenuniform, Grid(red,black))
        push!(shapes, arb)
    end

    for i=1:n
        r = 0.125
        ell = Ellipsoid(eye(3)/(r*r),[p[i,1],p[i,2], h[i]], blueuniform)
        push!(shapes, ell)
    end
    
    
    box = Box3(-3,3,-3,3,-3,3)

    # winsize=(1000,1000)
    # p = Cplot3.axes3(box,  fontsize=30, winsize=winsize, dest="knn3d.pdf")
    # camera = Camera(p.ax3, winsize)
    # for s in shapes
    #     init!(s, p.ax3)
    # end
    # X = raytrace_main(camera, lighting, shapes, winsize, true, false)
    # image = cairoimagefrommatrix(X)
    # drawimage(p.ctx, p.ax3, image)
    # Cplot.close(p)


    rt = raytrace(shapes, box; axisoptions3_fontsize=30, axisoptions3_width=500,
                  axisoptions3_height=500)
    d = draw(rt)
    save(d, plotpath("knn3d.pdf"))
    
end


end

