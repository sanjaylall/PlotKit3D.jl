
module ShapeTransform

using ..Basic3D: Vec3, Array23, Array32, Array33
using ..Axes3D: Box3, Box3f
using ..CoreRT: Texture, Shape, CoreRT


#export  hitnormal,  hittexture,  hitshape




##############################################################################
# Transform



mutable struct Transform{S1<:Shape} <: Shape
    shape::S1
    T::Array33
    Ttrans::Array33
    Tinv::Array33
    Tinvtrans::Array33
    c1::Vec3
    c::Vec3
    # points move from x to Tx+c
    # so that if a shape is defined by p(x) < 0
    # then the moved shape is p(T^{-1}(x-c)) < 0
    function Transform{S1}(shape::S1, Tm::Array{Float64,2}, c::Array{Float64,1})  where {S1<:Shape}
        T = Array33(Tm)
        Ttrans = Array33(Tm')
        Tinv = Array33(inv(Tm))
        Tinvtrans = Array33(inv(Tm)')
        c1 = Vec3(inv(Tm)*c)
        new(shape, T, Ttrans, Tinv, Tinvtrans, c1, Vec3(c))
    end
    function Transform{S1}(shape::S1, Tm::Array{Float64,2})  where {S1<:Shape}
        T = Array33(Tm)
        Ttrans = Array33(Tm')
        Tinv = Array33(inv(Tm))
        Tinvtrans = Array33(inv(Tm)')
        new(shape, T, Ttrans, Tinv, Tinvtrans, Vec3(0,0,0), Vec3(0,0,0))
    end
end

Transform(shape::S1, T, c) where {S1<:Shape} = Transform{S1}(shape,T,c)
Transform(shape::S1, T) where {S1<:Shape} = Transform{S1}(shape,T)

##############################################################################

function CoreRT.hitshape(s::Transform, o::Vec3, d::Vec3, tmin::Float64, tmax::Float64)
    d2 = s.Tinv*d
    ds = norm(d2)

    tenter, tleave, hd = hitshape(s.shape, s.Tinv*o - s.c1, (1/ds)*d2, tmin*ds, tmax*ds)
    return tenter/ds, tleave/ds, Hitdata(s, s.T*hd.point + s.c, hd.face, hd)
end

#
# and we need to know the true shape
#
function CoreRT.hitnormal(s::Transform, hitdata)
    n = hitnormal(hitdata.child.shape, hitdata.child)
    return normalize(s.Tinvtrans*n)
end

function CoreRT.hittexture(s::Transform, hitdata)
    return hittexture(hitdata.child.shape, hitdata.child)
end




end
