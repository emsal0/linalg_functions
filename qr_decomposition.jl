function qr_decompose(matx::Array{Float64,2})::Tuple{Array{Float64,2}, Array{Float64, 2}}
    vecs = [matx[:, i] for i in 1:size(matx)[2]]
    gs_vecs = gram_schmidt(vecs...)
    unitvecs = [unitvec(vec) for vec in gs_vecs]

    Q = zeros(size(matx))
    for i in 1:size(Q)[2]
        Q[:, i] = gs_vecs[i]
    end

    R = zeros(size(matx))
    for i in 1:size(R)[2]
        for j in 1:i
            term = dot(gs_vecs[j], vecs[i])
            R[j, i] = term
        end
    end

    return (Q,R)
end
