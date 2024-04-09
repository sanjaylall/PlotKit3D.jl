
module Shapes

using LinearAlgebra
using PlotKit

using ..Basic3D: Vec3, Array23, Array32, Array33
using ..Axes3D: Box3, Box3f
using ..CoreRT: Texture, Shape, CoreRT, nohit, Hitdata, gettexture

export Ellipsoid

const pk = PlotKit

# interface to Shape objects
# getdistancefromstripe is only required to support grid texture
#export hitnormal,  hittexture,  hitshape,  init!,   getdistancefromstripe
#export transform



# Shape Objects
export ArbitrarySolid, Polytope, Unitcube

# special case
export intersectunitcube



# Unitcube has no init! and no getdistancefromstripe
# Ellipsoid has no getdistancefromstripe
##############################################################################
# Ellipsoid

# set is (x-c)^T A (x-c) < 1
mutable struct Ellipsoid{S<:Texture} <: Shape
    A::Array33
    q::Vec3
    texture::S
end

Ellipsoid(A, q, texture::S) where {S<:Texture} = Ellipsoid{S}(A, q, texture)
Ellipsoid(A::Array{T}, q, texture::S) where {T<:Number, S<:Texture} = Ellipsoid{S}(Array33(A), q, texture)

#ellipsoid(a...) = Ellipsoid(a...)

function getparams(p::Ellipsoid)
    A = [vec(p.A.a)  vec(p.A.b) vec(p.A.c)]
    return A, vec(p.q)
end
    
function intersectellipsoid(A::Array33, q::Vec3, o, d)
    a = dot(d, A*d)
    b = 2*dot(d, A*(o-q))
    c = dot(o-q, A*(o-q))-1
    desc = b*b - 4 *a*c
    if desc < 0
        # no roots
        return -1.0, -1.0, 0
    end
    s = sqrt(desc)
    tenter = (-b-s)/(2*a)
    tleave = (-b+s)/(2*a)
    if tenter < 0
        return 0.0, tleave, 0
    end
    return tenter, tleave, 1
end

##############################################################################
# public Ellipsoid

function CoreRT.hitshape(s::Ellipsoid, o::Vec3, d::Vec3, tmin::Float64, tmax::Float64)
    if tmin < 0
        return nohit(s)
    end
    if tmax < tmin
        return nohit(s)
    end
    tenter, tleave, face= intersectellipsoid(s.A, s.q, o, d)
    if tleave >= tmin && tenter <= tmax
        thit = max(tenter, tmin)
        point = o + thit*d
        return thit, tleave, Hitdata(s, point, face)
    end
    # did not hit
    return nohit(s)
end

function CoreRT.hitnormal(s::Ellipsoid, hitdata)
    return normalize(s.A*(hitdata.point - s.q))
end

function CoreRT.hittexture(s::Ellipsoid, hitdata)
    return gettexture(s, s.texture, hitdata.point)
end


# if xa in axis coords, then to get cube coords xc
#   xc = T * xa + c  where T = diagm(scale) c = origin
#
# we have (xa - q)^T A (xa - q) < 1
#
# which gives new constraints on xc
#
function CoreRT.transform(p::Ellipsoid, T::Array{Float64,2}, c::Array{Float64,1})
    A, q = getparams(p)
    newA = inv(T')*A*inv(T)
    newq = T*q + c
    return Ellipsoid(newA, newq, p.texture)
end



##############################################################################
# ArbitrarySolid


mutable struct ArbitrarySolid{S1<:Texture,S2<:Texture} <: Shape
    f_isinside  # function that checks if inside
    f_normal    # function that returns outward normal at point.
                # gradient should not be normalized (for efficiency)
    tol         # tolerance 
    texture1::S1
    texture2::S2
end

#arbitrarysolid(a...; kw...) = ArbitrarySolid(a...; kw...)

function intersectarbitrary(s::ArbitrarySolid, o::Vec3, d::Vec3,  tmin::Float64, tmax::Float64)
    t = tmin
    olddiff = 0.0
    while t <= tmax
        q = o + d*t
        if s.f_isinside(q)
            return t
        end
        t = t + s.tol
    end
    return -1.0
end

##############################################################################
# public ArbitrarySolid


function CoreRT.transform(p::ArbitrarySolid, T::Array{Float64,2}, c::Array{Float64,1})
    invTfast = Array33(inv(T))
    invTTfast = Array33(Matrix(inv(T)'))
    cfast = Vec3(c)
    f = xc -> invTfast*(xc - cfast)
    fin_new = xc -> p.f_isinside(f(xc))
    
    function fnrm_new(xc)
        g, face = p.f_normal(f(xc))
        return invTTfast * g, face
    end
    
    return ArbitrarySolid(fin_new, fnrm_new, p.tol, p.texture1, p.texture2)
end



function CoreRT.hitshape(s::ArbitrarySolid, o::Vec3, d::Vec3, tmin::Float64, tmax::Float64)
    if tmin < 0
        return nohit(s)
    end
    if tmax < tmin
        return nohit(s)
    end
    tenter = intersectarbitrary(s, o, d, tmin, tmax)
    if tenter >= 0
        point = o + tenter*d

        # set tleave to be tenter since we don't know what it does
        tleave = tenter

        # get color
        n,face = s.f_normal(point)
        return tenter, tleave, Hitdata(s, point, face)
    end
    # did not hit
    return nohit(s)
end


function CoreRT.hitnormal(s::ArbitrarySolid, hitdata)
    point = hitdata.point
    n, face = s.f_normal(point)
    return normalize(n)
end

function CoreRT.hittexture(s::ArbitrarySolid, hitdata)
    if hitdata.face == 1
        return gettexture(s, s.texture1, hitdata.point)
    else
        return gettexture(s, s.texture2, hitdata.point)
    end
end



# q is point on surface
#
function CoreRT.getdistancefromstripe(s::ArbitrarySolid, q::Vec3, spacingx, spacingy)
    distfromstripe = min(abs(xa/spacingx - round(xa/spacingx)),
                          abs(ya/spacingy - round(ya/spacingy)))
    return distfromstripe
end



##############################################################################
# polytopes

mutable struct Polytope{S<:Texture} <: Shape
    At::Vector{Vec3}  # stored as a list of rows
    b::Vector{Float64}
    textures::Array{S,1}
end

# do we need this?
Polytope(A, b, textures::Array{S,1}) where {S<:Texture} = Polytope{S}(A,b,textures)

# Convert from sing Matrix form to faster Vec3 representation
Polytope(A::Matrix{T}, b, textures::Array{S,1}) where {T <: Number, S<:Texture} = Polytope(Vec3[z[:] for z in eachrow(A)], b, textures)

#polytope(a...; kw...) = Polytope(a...; kw...)

function getparams(p::Polytope)
    # convert to usual matrix form
    A = reduce(vcat, vec(q)' for q in p.At)
    return A, p.b
end


# returns
#
#
#   tenter = -1.0   if ray misses polytope
#   tenter =  0.0   if ray origin is inside polytope
#   tenter >  0     if ray hits polytope then tenter is the time at which it does so
#
#   tleave = -1.0        if ray misses polytope
#   tleave >= tenter     otherwise
#
#   face is the number of the face
#   face = 0 if the ray starts inside the polytope
#
#
# returns -1 if no intersection between ray o+td and polytope
# else return time t when ray enters polytope
function intersectpolytope(pAt::Vector{Vec3},  pb::Array{Float64, 1}, o::Vec3, d::Vec3)
    # number of faces
    m = length(pAt)
        
    tenter = -Inf
    tleave = Inf
    face = 1
    for i=1:m
        vn = pb[i] - dot(pAt[i], o)
        vd = dot(pAt[i], d)
        if vd == 0
            # missed face i since ray is parallel to it
            if vn < 0
                # ray is outside halfspace i, hence misses polytope
                return -1.0,-1.0, i
            end
        elseif vd > 0
            # ray is traveling away from face i from outside
            # or towards it from the inside 
            #
            t = vn/vd
            if t < 0
                # vn/vd will be negative if the ray is traveling away from the face
                # from the outside
                return -1.0, -1.0, i
            end
            tleave = min(t, tleave)
        else
            # ray is traveling toward face i from outside
            # or away from it from inside
            # find time of intersection
            t = vn/vd
            if t>tenter
                tenter = t
                face = i
            end
        end
        if tleave < tenter
            return -1.0, -1.0, face
        end
    end
    if tenter < 0
        # we are inside the polytope
        #println("inside")
        return 0.0, tleave, 0
    end
    # we will hit polytope face at time tenter >= 0, and leave at tleave >= tenter
    return tenter, tleave, face
end


##############################################################################
# public Polytope

#  returns -1.0, -1.0 if there is no intersection for  tmin <= t <= tmax
#
#  returns  tenter, tleave    if there is an intersection
#           with   tmin <= tenter <= tmax
#           and    tenter <= tleave
#
#  face = 0 if hit is on interior
#  
#
function CoreRT.hitshape(p::Polytope, o::Vec3, d::Vec3, tmin::Float64, tmax::Float64)
    if tmin < 0
        return nohit(p)
    end
    if tmax < tmin
        return nohit(p)
    end
    tenter, tleave, face = intersectpolytope(p.At, p.b, o, d)
    if tleave >= tmin && tenter <= tmax
        thit = max(tenter, tmin)
        point = o + thit*d
        return thit, tleave, Hitdata(p, point, face)
    end
    # did not hit
    return nohit(p)
end

function CoreRT.hitnormal(p::Polytope, hitdata)
    face = hitdata.face
    return normalize(p.At[face])
end

function CoreRT.hittexture(p::Polytope, hitdata)
    return gettexture(p, p.textures[hitdata.face], hitdata.point)
end

# if xa in axis coords, then to get cube coords xc
#   xc = T * xa + c  where T = diagm(scale) c = origin
#
# 
# so if we have A*xa <=  b
# that is  A*(inv(T)* (xc-c) ) <= b
# that is  A*inv(T)*xc <= b + A*inv(D)*c
#
# so At -> inv(T')*At
# and b ->  b + A*inv(T)*c
# 
function CoreRT.transform(p::Polytope, T::Array{Float64,2}, c::Array{Float64,1})
    A, b = getparams(p)
    newA = A*inv(T)
    newb = p.b + newA*c
    return Polytope(newA, newb, p.textures)
end


##############################################################################
# Unitcube



# returns
#
#
#   tenter = -1.0   if ray misses polytope
#   tenter =  0.0   if ray origin is inside polytope
#   tenter >  0     if ray hits polytope then tenter is the time at which it does so
#
#   tleave = -1.0        if ray misses polytope
#   tleave >= tenter     otherwise
#
#   face is the number of the face the ray hits
#   face = 0 if the ray starts inside the polytope
#
#
# returns -1 if no intersection between ray o+td and polytope
# else return time t when ray enters polytope
#
function intersectunitcube(o::Vec3, d::Vec3)

    tenter = -Inf
    tleave = Inf
    face = 1

    # face x = 1
    # vn is normal signed distance of origin from face
    # vn > 0 means origin on inside-side of face
    # vn < 0 means origin is outside cube, separated from it by the facial plane
    vn = 1.0 - o.x
    vd = d.x
    if vd == 0
        # x-component of velocity is zero
        if  vn < 0
            # x component of position > 1
            # so we don't hit the cube
            return -1.0, -1.0, 1
        end
    elseif vd > 0
        # x-component of velocity > 0
        # so set t to be time at which
        # ray crosses the x=1 plane
        t = vn/vd
        if t < 0
            # origin is on outside side of x=1 face, heading away from it
            # so we don't hit the cube
            return -1.0, -1.0, 1
        end
        # origin is on inside side of x=1 face, heading toward it
        # t is the time of intersection
        # so we know tleave <= t
        tleave = min(t, tleave)
    else
        t = vn/vd
        # now we know tenter >= t
        tenter = max(t, tenter)
    end
    if tleave < tenter
        return -1.0, -1.0, face
    end


    
    vn = 1.0 + o.x
    vd = -d.x
    if vd == 0
        if vn < 0
            return -1.0, -1.0, 2
        end
    elseif vd > 0
        t = vn/vd
        if t < 0
            return -1.0, -1.0, 2
        end
        tleave = min(t, tleave)
    else
        t = vn/vd
        if t > tenter
            face = 2
            tenter = t
        end
    end
    if tleave < tenter
        return -1.0, -1.0, face
    end

    

    vn = 1.0 - o.y
    vd = d.y
    if vd == 0
        if vn < 0
            return -1.0, -1.0, 3
        end
    elseif vd > 0
        t = vn/vd
        if t < 0
            return -1.0, -1.0, 3
        end
        tleave = min(t, tleave)
    else
        t = vn/vd
        if t > tenter
            face = 3
            tenter = t
        end
    end
    if tleave < tenter
        return -1.0, -1.0, face
    end

    

    vn = 1.0 + o.y
    vd = -d.y
    if vd == 0
        if vn < 0
            return -1.0, -1.0, 4
        end
    elseif vd > 0
        t = vn/vd
        if t < 0
            return -1.0, -1.0, 4
        end
        tleave = min(t, tleave)
    else
        t = vn/vd
        if vn/vd > tenter
            face = 4
            tenter = t
        end
    end
    if tleave < tenter
        return -1.0, -1.0, face
    end


    
    vn = 1.0 - o.z
    vd = d.z
    if vd == 0
        if vn < 0
            return -1.0, -1.0, 5
        end
    elseif vd > 0        
        t = vn/vd
        if t < 0
            return -1.0, -1.0, 5
        end
        tleave = min(t, tleave)
    else
        t = vn/vd
        if t > tenter
            face = 5
            tenter = t
        end
    end
    if tleave < tenter
        return -1.0, -1.0, face
    end


    
    vn = 1.0 + o.z
    vd = -d.z
    if vd == 0
        if vn < 0
            return -1.0, -1.0, 6
        end
    elseif vd > 0
        t = vn/vd
        if t < 0
            return -1.0, -1.0, 6
        end
        tleave = min(t, tleave)
    else
        t = vn/vd
        if t > tenter
            face = 6
            tenter = t
        end
    end
    if tleave < tenter
        return -1.0, -1.0, face
    end

    if tenter < 0
        return 0.0, tleave, 0
    end
    return tenter, tleave, face
end



mutable struct Unitcube{S<:Texture} <: Shape
    texture::S
end




##############################################################################
# public Unitcube

# hitshape only returns hits for which
# tmin <= thit <= tmax
# but tleave may be greater than tmax
#

function CoreRT.hitshape(s::Unitcube, o::Vec3, d::Vec3, tmin::Float64, tmax::Float64)
    if tmin < 0
        return nohit(s)
    end
    if tmax < tmin
        return nohit(s)
    end
    tenter, tleave, face = intersectunitcube(o, d)
    if tleave >= tmin && tenter <= tmax
        thit = max(tenter, tmin)
        point = o + thit*d
        return thit, tleave, Hitdata(s, point, face)
    end
    # did not hit
    return nohit(s)
end
    

function CoreRT.hitnormal(s::Unitcube, hitdata)
    At = (Vec3(1,0,0),
          Vec3(-1,0,0),
          Vec3(0,1,0),
          Vec3(0,-1,0),
          Vec3(0,0,1),
          Vec3(0,0,1))
    return At[hitdata.face]
end

function CoreRT.hittexture(s::Unitcube, hitdata)
    return gettexture(s, s.texture, hitdata.point)
end

end
