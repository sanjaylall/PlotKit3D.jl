
module ShapeIntersection

using ..Basic3D: Vec3, Array23, Array32, Array33
using ..Axes3D: Box3, Box3f
using ..CoreRT: Texture, Shape, CoreRT


##############################################################################
# Constructive Solid Geometry

mutable struct Intersection{S1<:Shape,S2<:Shape} <: Shape
    s1::S1
    s2::S2
    function Intersection{S1,S2}(s1::S1, s2::S2) where {S1<:Shape,S2<:Shape}
        new(s1,s2)
    end
end

#Intersection{S1<:Shape,S2<:Shape}(s1::S1, s2::S2) = Intersection{S1,S2}(s1,s2)

# public
function CoreRT.hitshape(s::Intersection, o::Vec3, d::Vec3,
                        tmin::Float64, tmax::Float64)
    tenter1, tleave1, hitdata1 = hitshape(s.s1, o, d, tmin, tmax)
    if tenter1 < 0
        return nohit(s)
    end
    tenter2, tleave2, hitdata2 = hitshape(s.s2, o+tenter1*d, d, 0.0, min(tleave1, tmax)-tenter1)
    if tenter2 < 0
        return nohit(s)
    end
    if hitdata2.face == 0
        return tenter2+tenter1, min(tleave1, tleave2+tenter1), Hitdata(s, hitdata1.point,
                                                                       hitdata1.face, hitdata1)
    end
    return tenter2+tenter1, min(tleave1, tleave2+tenter1),  Hitdata(s, hitdata2.point,
                                                                    hitdata2.face, hitdata2)
end


function CoreRT.hitnormal(s::Intersection, hd)
    n::Vec3 =  hitnormal(hd.child.shape, hd.child)
    return n
end

function CoreRT.hittexture(s::Intersection, hd)
    shd = hd.child  # abstract type Hitdata
    ss = shd.shape  # hence abstract type Shape

    t::Material, r::Bool = hittexture(ss, shd)
    return t, r
end


end
