

module RTcairo

using PlotKit

using ..RayTracer: Camera
# This file contains the code that the raytracer uses
# to transfer the raytraced image onto a Cairo surface
#



##############################################################################
# Utility functions



# apply the bitmap to the axes
function PlotKit.drawimage(ctx, camera::Camera, p::Pik)
    hexagon = get_hexagon(camera.am3)
    #  drawimage(ctx, p, Point(0,0))   # see whole image
    drawimage_to_mask(ctx, p, hexagon,
                      camera.renderwidth/camera.width, camera.renderheight/camera.height)
end


fftoi(a::Float64) = round(UInt32, a * 255)

pack(r::Float64, g::Float64, b::Float64, a::Float64) = fftoi(a) << 24 + fftoi(r) << 16 + fftoi(g) << 8 + fftoi(b)

# pixelmap(r,g,b) = pack(r,g,b,1.0)
function pixelmap(r,g,b)
    if r + g + b < 3
        return pack(r, g, b, 1.0)
    end
    return pack(1.0, 1.0, 1.0, 0.05)
end
    
function cairoimagefrommatrix(X::Array{Float64,3})
    Z = UInt32[pixelmap(X[x,y,1], X[x,y,2], X[x,y,3]) for x = 1:size(X,1), y=1:size(X,2)]
    pik = Pik(Z, size(X,1), size(X,2))
end
    

end


