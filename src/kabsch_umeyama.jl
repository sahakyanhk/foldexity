using LinearAlgebra
using Statistics
using Clustering
using Distances


#kabsch-umeyama
function translate_to_centroid(coord_matrix)
    # Normalises the molecular coordinates by centering them.
    center = [mean(coord_matrix[:,1]);mean(coord_matrix[:,2]);mean(coord_matrix[:,3])]
    centroid = transpose(center)
    translated_geom = broadcast(-,coord_matrix,centroid)
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


function align_rmsd(m1,m2)

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
function fxity_kabsh(coordmatrix, cutoff = 1.0)    
    try        
        n = length(coordmatrix)  # Change this to the desired size
        matrix = zeros(Float64, n, n)

        for i = 1:n # Fill the upper triangle
            for j = i+1:n  # Ensure j >= i for the upper triangle
                matrix[i, j] = align_rmsd(coordmatrix[i], coordmatrix[j])
            end
        end

        matrix += matrix' #make a symmetric matrix
        aver_rmsd = sum(matrix) / (n * n)

        cl = hclust(matrix, linkage=:complete)
        results = cutree(cl, h=cutoff) 
        nclusts = length(unique(results))
        norm_nclusts = nclusts / n

        fxity = shannon(results, 1)

        return fxity, aver_rmsd, nclusts, norm_nclusts, n, matrix 
        
    catch 
        return 0
    end
end
