using LinearAlgebra
using Statistics
using Clustering
using Distances

include("entropy.jl")


#alighn and report rmsd or TMscore
# function kabsch_umeyama(m1::Matrix,m2::Matrix)

#     function translate_to_centroid(coords)
#         # Normalises the molecular coordinates by centering them.
#         center = [mean(coords[:,1]);mean(coords[:,2]);mean(coords[:,3])]
#         centroid = transpose(center)
#         translated_geom = broadcast(-,coords,centroid)
#         return translated_geom
#     end
    
#     function cross_covariance_matrix(Pmatrix,Qmatrix)
#         # Cross covariance matrix gives measure of variability between two matrices
#         CCmatrix = transpose(Pmatrix) * Qmatrix
#         return CCmatrix
#     end
    
#     function optimal_rotation_matrix(CCmatrix)
#         # Returns 3x3 matrix that can be applied to P to get Q
#         ORmatrix = sqrt(transpose(CCmatrix)*CCmatrix)*inv(CCmatrix)
#         return ORmatrix
#     end    

#     normalisedP = translate_to_centroid(m1)
#     normalisedQ = translate_to_centroid(m2)

#     xcov = cross_covariance_matrix(normalisedP,normalisedQ)
#     orot = optimal_rotation_matrix(xcov)

#     num_atoms = size(normalisedP)[1]
#     rotated = zeros(Float64,num_atoms,3)

#     rotated  = (orot*normalisedP')'

#     RMSD_value = norm(rotated-normalisedQ)

#     #TMscore
#     distances = colwise(euclidean, rotated', normalisedQ')
#     L0 = num_atoms
#     d0 = 1.24 * cbrt(L0 - 15) - 1.8
#     tmscore = (1/L0) * sum([1 / (1 + ((d/d0) ^2)) for d in distances])
    
#     return RMSD_value
# end

function kabsch_umeyama(P::AbstractMatrix, Q::AbstractMatrix)
    # Ensure matrices have the same dimensions
    if size(P) != size(Q)
        error("Input matrices must have the same dimensions.")
    end

    N, _ = size(P)

    # 1. Center the point sets
    centroid_P = sum(P, dims=1) ./ N
    centroid_Q = sum(Q, dims=1) ./ N
    P_centered = P .- centroid_P
    Q_centered = Q .- centroid_Q

    # 2. Compute the covariance matrix
    H = P_centered' * Q_centered

    # 3. Perform Singular Value Decomposition (SVD)
    svd_result = svd(H)
    U, _, V = svd_result.U, svd_result.S, svd_result.V

    # 4. Calculate the optimal rotation matrix, correcting for reflection
    d = sign(det(V * U'))
    correction_matrix = diagm([1, 1, d])
    R = V * correction_matrix * U'

    # 5. Apply the rotation and translation to align P with Q
    # P_aligned = (R * P_centered')' .+ centroid_Q
    # For older Julia versions without broadcasting `.+` for matrices
    P_aligned = (R * P_centered')' .+ repeat(centroid_Q, N, 1)


    # 6. Calculate the RMSD
    diff = P_aligned .- Q
    rmsd = sqrt(sum(diff.^2) / N)

    return rmsd
end



#calculate all-vs-all kabsch rmsd for fragments in the matrix
function fxity_kabsh(xyzcoords, cutoff = 1.0)    
    try        
        nfrags = length(xyzcoords)  # Change this to the desired size
        matrix = zeros(Float64, nfrags, nfrags)

        for i = 1:nfrags # Fill the upper triangle
            for j = i+1:nfrags  # Ensure j >= i for the upper triangle
                matrix[i, j] = kabsch_umeyama(xyzcoords[i], xyzcoords[j])
            end
        end

        matrix += matrix' #make a symmetric matrix
        aver_rmsd = sum(matrix) / (nfrags * nfrags)

        cl = hclust(matrix, linkage=:complete)
        results = cutree(cl, h=cutoff) 
        nclusts = length(unique(results))
        norm_nclusts = nclusts / nfrags

        foldexity = entropy_shannon(results, 1)

        return foldexity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix 
        
    catch err
        print(err)
        return 0
    end
end


