
module ShapeSurface

using LinearAlgebra

using ..Basic3D: Vec3, Array23, Array32, Array33
using ..Axes3D: Box3, Box3f
using ..CoreRT: Texture, Shape, CoreRT, nohit, Hitdata, gettexture

export Surface


##############################################################################

function unnormalize(xn::Float64, xmin::Float64, xmax::Float64)
    tx = (xmax - xmin) / 2
    cx = xmin + tx
    xu = cx + tx * xn
end

function normalize(xu::Float64, xmin::Float64, xmax::Float64)
    tx = (xmax - xmin) / 2
    cx = xmin + tx
    xn = (xu - cx) / tx
end

function unnormalize(xn, yn, cube::Box3f) 
    xu = unnormalize(xn, cube.xmin, cube.xmax)
    yu = unnormalize(yn, cube.ymin, cube.ymax)
    return xu, yu
end

function normalize(xu, yu, cube::Box3f)
    xn = normalize(xu, cube.xmin, cube.xmax)
    yn = normalize(yu, cube.ymin, cube.ymax)
    return xn, yn
end


struct Height
    zfun::Function
    sampleheight::Bool
    m::Int64
    H::Array{Float64,2}    # height samples
    cube::Box3f
end

function Height(zfun, sampleheight::Bool, cube::Box3)
    m = 100
    H = createheightsamples(zfun, m, Box3f(cube))
    return Height(zfun, sampleheight, m, H, cube)
end

struct GradientHeight
    dzfun::Function        # gradient of height, in axis coords
    samplegradient::Bool
    m::Int64
    Gx::Array{Float64,2}   # samples of dz/dx
    Gy::Array{Float64,2}   # samples of dz/dy
    cube::Box3f
end

function GradientHeight(dzfun, samplegradient::Bool, cube::Box3)
    m = 100
    Gx, Gy = creategradientsamples(dzfun, m, Box3f(cube))
    return GradientHeight(dzfun, samplegradient, m, Gx, Gy, Box3f(cube))
end

##############################################################################
# Surface
#
# A surface is defined by a function z = h(x,y)
#
mutable struct Surface{S1<:Texture,S2<:Texture} <: Shape
    zfun::Function         # z = height(x,y) defines surface, in axis coords
    dzfun::Function        # gradient of height, in axis coords
    texture1::S1           # top
    texture2::S2           # bottom
    sampleheight::Bool
    samplegradient::Bool

    height::Height         # object for caching zfun
    gradientheight::GradientHeight 
    lipschitz::Float64     # Lipschitz constant for f
    m::Int64               # 2m+1 by 2m+1 sample points
    H::Array{Float64,2}    # height samples
    Gx::Array{Float64,2}   # samples of dz/dx
    Gy::Array{Float64,2}   # samples of dz/dy

    
    function Surface{S1,S2}(zf, dz, texture1::S1, texture2::S2, sampleheight, samplegradient) where {S1<:Texture,S2<:Texture}
        new(zf, dz, texture1, texture2, sampleheight, samplegradient)
    end
end

function Surface(zfun, dzfun, texture1::S1, texture2::S2,
                 sampleheight, samplegradient) where {S1<:Texture,S2<:Texture}
    Surface{S1,S2}(zfun, dzfun, texture1, texture2, sampleheight, samplegradient)
end

# defaults
function Surface(zfun, dzfun; sampleheight = true, samplegradient = true,
                 texture1 = nothing, texture2 = nothing)
    if isnothing(texture1)
        texture1 = Grid(:red)
    end
    if isnothing(texture2)
        texture2 = Grid(:green)
    end
    s = Surface(zfun, dzfun, texture1, texture2,
                      sampleheight, samplegradient)
end

##############################################################################
# height interp

function (sh::Height)(xu::Float64, yu::Float64)
    if sh.sampleheight
        xn, yn = normalize(xu, yu, sh.cube)
        return binterp(sh.H, 2 + sh.m + xn * sh.m, 2 + sh.m + yn * sh.m)
    end
    # This type q::Float64 is here, because the closure used for transforming
    # zfun hides the type information
    # and so dzfun becomes type unstable.
    # Putting a type on q prevents gradient() becoming
    # type unstable also.
    q::Float64 = sh.zfun(xu, yu)
    return q
end


function (sh::GradientHeight)(xu::Float64, yu::Float64)
    if sh.samplegradient
        xn, yn = normalize(xu, yu, sh.cube)
        dzdx = binterp(sh.Gx, 2+sh.m+xn*sh.m, 2+sh.m+yn*sh.m)
        dzdy = binterp(sh.Gy, 2+sh.m+xn*sh.m, 2+sh.m+yn*sh.m)
        return dzdx, dzdy
    end
    # This is here, because the closure used for transforming
    # dzfun hides the type information
    # and so dzfun becomes type unstable.
    # Putting a type on q prevents gradient() becoming
    # type unstable also.
    q::Tuple{Float64, Float64} = sh.dzfun(xu, yu)
    return q
end


##############################################################################

# bilinear interpolation
function binterp(H, i, j)
    ia = floor(i)
    ja = floor(j)

    it = i - ia
    jt = j - ja

    return  (1-it)*H[Int64(ia), Int64(ja)]*(1-jt) +  it*H[Int64(ia)+1, Int64(ja)]*(1-jt) + (1-it)*H[Int64(ia), Int64(ja)+1]*jt +  it*H[Int64(ia)+1, Int64(ja)+1]*jt
    
end

# 
# function gradient(s::Surface,  x::Float64, y::Float64)
#     if s.samplegradient
#         dzdx = binterp(s.Gx, 2+s.m+x*s.m, 2+s.m+y*s.m)
#         dzdy = binterp(s.Gy, 2+s.m+x*s.m, 2+s.m+y*s.m)
#         return dzdx, dzdy
#     end
#     # This is here, because the closure used for transforming dzfun hides the type information
#     # and so dzfun becomes type unstable.
#     # Putting a type on q prevents gradient() becoming
#     # type unstable also.
#     q::Tuple{Float64, Float64} = s.dzfun(x, y)
#     return q
# end
# 
#function height(s::Surface, x::Float64, y::Float64)
#    if s.sampleheight
#        return binterp(s.H, 2+s.m+x*s.m, 2+s.m+y*s.m)
#    end
#    q::Float64 = s.zfun(x, y)
#    return  q
#end

function surfacenormal(s::Surface, x::Float64, y::Float64)
    dzdx, dzdy = s.gradientheight(x, y)
    n = sqrt(1+dzdx*dzdx + dzdy*dzdy)
    return  Vec3(-dzdx/n, -dzdy/n, 1/n)
end


function isdiagonal(T::Array{Float64,2})
    for i=1:size(T,1)
        for j=1:size(T,2)
            if i != j && T[i,j] != 0
                return false
            end
        end
    end
    return true
end

# return time, side, point
#
# if intersection occurs at time > tmax returns -1.0
#
function intersectsurface(s::Surface, o::Vec3, d::Vec3,  tmin::Float64, tmax::Float64)
    t = tmin
    olddiff = 0.0
    L = s.lipschitz
    while t <= tmax
        q = o + d*t
        h = s.height(q.x, q.y)
        diff = h - q.z
        if abs(diff) < 0.001
            if olddiff > 0
                return t, 2, q
            end
            return t, 1, q
        end
        olddiff = diff
        t = t + abs(diff)/L
    end
    return -1.0, 1, Vec3(0,0,0)
end

##############################################################################
# initialization




function lipschitz(dz::Function)
    # compute lipschitz bound
    L::Float64 = 0.0
    for x = -1:0.01:1
        for y = -1:0.01:1
            g1, g2 = dz(x, y)
            g = sqrt(g1*g1 + g2*g2)
            if g > L 
                L = g
            end
        end
    end
    return L
end

function createheightsamples(zf::Function, m::Int64, cube::Box3f)
    # want samples from -m to m
    # so number of samples is 
    N = 2*m+1  
    # add one extra row/col each side for interpolation
    H = zeros(N+2,N+2)
    for i=0:N+1
        for j=0:N+1
            xn = (i-m-1)/m  # i = 1 maps to xn = -1,  and i = N maps to xn = 1
            yn = (j-m-1)/m
            x, y = unnormalize(xn, yn, cube) 
            H[i+1,j+1] = zf(x, y)
        end
    end
    return H
end


function creategradientsamples(dz::Function, m::Int64, cube::Box3f)
    N = 2*m+1  
    Gx = zeros(N+2,N+2)
    Gy = zeros(N+2,N+2)
    for i=0:N+1
        for j=0:N+1
            xn = (i-m-1)/m  # i = 1 maps to xn = -1,  and i = N maps to xn = 1
            yn = (j-m-1)/m
            x, y = unnormalize(xn, yn, cube) 
            gx, gy =  dz(xn, yn)
            Gx[i+1,j+1] = gx
            Gy[i+1,j+1] = gy
        end
    end
    return Gx, Gy
end

function initialize(s::Surface)
    # add 1 to lipschitz since we are solving for root g(t) = f(r(t))-r(t)_3
    # and so need lipschitz of g
    s.lipschitz = lipschitz(s.dzfun) + 1
    s.m = 100  # -m to m
#    s.H = createheightsamples(s.zfun, s.m)
#    s.Gx, s.Gy = creategradientsamples(s.dzfun, s.m)
end

##############################################################################
# transforms

# if xa in axis coords, then to get cube coords xc
#   xc = T * xa + c  where T = diagm(scale) c = origin
#
# Here T is diagm(tx, ty, tz) and c = [qx, qy, qz]
#
# These closures are slow, but because heights
# and gradients are cached it doesn't affect performance
# in a noticeable way.
function transform_fn(f, tx, ty, tz, qx, qy, qz)
    function f_new(cx, cy)
        ax = (cx - qx) / tx
        ay = (cy - qy) / ty
        cz = f(ax, ay)
        az::Float64 = tz * cz + qz
        return az
    end
    return f_new
end

function transform_grad_fn(dz, tx, ty, tz, qx, qy, qz)
    function dz_new(cx, cy)
        ax = (cx - qx) / tx
        ay = (cy - qy) / ty
        g1, g2 = dz(ax, ay)
        return tz * g1 / tx,  tz * g2 / ty
    end
    return dz_new
end



##############################################################################
# Surface public functions

function CoreRT.hitshape(s::Surface, o::Vec3, d::Vec3, tmin::Float64, tmax::Float64)
    if tmin < 0
        return nohit(s)
    end
    if tmax < tmin
        return nohit(s)
    end
    tenter, face, point = intersectsurface(s, o, d, tmin, tmax)
    if tenter >= tmin
        tleave = tenter
        return tenter, tleave, Hitdata(s, point, face)
    end
    return nohit(s)

end


function CoreRT.transform(s::Surface, T::Array{Float64,2}, c::Array{Float64,1})
    @assert isdiagonal(T)
    zfun_new = transform_fn(s.zfun, T[1,1], T[2,2], T[3,3], c[1], c[2], c[3])
    dzfun_new = transform_grad_fn(s.dzfun, T[1,1], T[2,2], T[3,3], c[1], c[2], c[3])
    snew = Surface(zfun_new, dzfun_new, s.texture1, s.texture2, s.sampleheight, s.samplegradient)
    return snew
end

# need to account for limits here
function CoreRT.addlimits(s::Surface, cube::Box3)
    initialize(s)
    s.height = Height(s.zfun, s.sampleheight, cube)
    s.gradientheight = GradientHeight(s.dzfun, s.samplegradient, cube)
end



function CoreRT.hitnormal(s::Surface, hitdata)
    point = hitdata.point
    if hitdata.face == 1
        return surfacenormal(s, point.x, point.y)
    end
    return -surfacenormal(s, point.x, point.y)
end

function CoreRT.hittexture(s::Surface, hitdata)
    if hitdata.face == 1
        return gettexture(s, s.texture1, hitdata.point)
    end
    return gettexture(s, s.texture2, hitdata.point)
end



#
# q is point on surface
#
function CoreRT.getdistancefromstripe(s::Surface, q::Vec3, spacingx, spacingy)
    # decide width of stripe normalized by gradient
    # otherwise stripes will be thicker on slopes
    # this is an approximation
    g1,g2 = s.gradientheight(q.x, q.y)
    sx = sqrt(g1*g1+1)
    sy = sqrt(g2*g2+1)
    distfromstripe = min(abs(sx*(q.x / spacingx - round(q.x / spacingx))),
                         abs(sy*(q.y / spacingy - round(q.y / spacingy))))
    return distfromstripe
end




end

