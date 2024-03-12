


# function color(name::AbstractString)
#     c=Dict{AbstractString,Any}()
    
#     c["white"] = [1 1 1]
#     c["red"]   = [1 0 0]
#     c["green"]  = [0 1 0]
#     c["darkgreen"]  = [0 100/255 0]

#     c["blue"]   = [0 0 1]
#     c["cyan"] = [0 1 1]
#     c["magenta"]   = [1 0 1]
#     c["yellow"] = [1 1 0]
#     c["black"] = [0 0 0]

#     return c[name]
# end

# color(c::Float64) = c*[1 1 1]
# color(c::Tuple) = c
# color(c::Array) = c

# # alpha
# colora(c::AbstractString) = [color(c) 1]
# colora(c::Array) = c
# colora(c::Tuple) = c
