using LinearAlgebra
using Statistics
using Clustering
using Distances

include("entropy.jl")


#alighn and report rmsd or TMscore
function kabsch_umeyama(m1::Matrix,m2::Matrix)

    function translate_to_centroid(coords)
        # Normalises the molecular coordinates by centering them.
        center = [mean(coords[:,1]);mean(coords[:,2]);mean(coords[:,3])]
        centroid = transpose(center)
        translated_geom = broadcast(-,coords,centroid)
        return translated_geom
    end
    
    function cross_covariance_matrix(Pmatrix,Qmatrix)
        # Cross covariance matrix gives measure of variability between two matrices
        CCmatrix = transpose(Pmatrix) * Qmatrix
        return CCmatrix
    end
    
    function optimal_rotation_matrix(CCmatrix)
        # Returns 3x3 matrix that can be applied to P to get Q
        ORmatrix = sqrt(transpose(CCmatrix)*CCmatrix)*inv(CCmatrix)
        return ORmatrix
    end    

    normalisedP = translate_to_centroid(m1)
    normalisedQ = translate_to_centroid(m2)

    xcov = cross_covariance_matrix(normalisedP,normalisedQ)
    orot = optimal_rotation_matrix(xcov)

    num_atoms = size(normalisedP)[1]
    rotated = zeros(Float64,num_atoms,3)

    rotated  = (orot*normalisedP')'

    RMSD_value = norm(rotated-normalisedQ)

    #TMscore
    distances = colwise(euclidean, rotated', normalisedQ')
    L0 = num_atoms
    d0 = 1.24 * cbrt(L0 - 15) - 1.8
    tmscore = (1/L0) * sum([1 / (1 + ((d/d0) ^2)) for d in distances])
    
    return RMSD_value
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


function dihedaral(xyzmatrix::Matrix{T}) where {T}

    if size(xyzmatrix) != (4, 3)

        error("Input matrix must be of size (4, 3)!")

    end

    p1 = xyzmatrix[1,1:3]
    p2 = xyzmatrix[2,1:3]
    p3 = xyzmatrix[3,1:3]
    p4 = xyzmatrix[4,1:3]

    b1 = p2 .- p1
    b2 = p3 .- p2
    b3 = p4 .- p3

    # Normalize the vectors
    b1 /= norm(b1)
    b2 /= norm(b2)
    b3 /= norm(b3)

    # Compute normals
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)

    # Compute the angle
    x = dot(n1, n2)
    y = dot(cross(n1, n2), b2)

    angle = atan(y, x)
    return angle * (180 / π)  # Convert to degrees

end