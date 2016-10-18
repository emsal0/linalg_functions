function project(vec1::Array{Float64}, vec2::Array{Float64})::Array{Float64}
    return (dot(vec1,vec2)/dot(vec2,vec2)) * vec2
end

function unitvec(vec::Array{Float64})::Array{Float64}
    return (1/norm(vec)) * vec
end

function gram_schmidt(vecs...)
    u = [unitvec(vecs[1])]

    for i in vecs[2:end]
        new_u = i
        for j in u
            p = project(j, i)
            new_u -= p
        end
        push!(u, unitvec(new_u))
    end

    return u
end
