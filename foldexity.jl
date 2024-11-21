using LinearAlgebra
using Statistics
using Glob


function pdb2matrix(pdb_file::String)
    #pdb_txt = readdlm(pdb_file)
    x,y,z = [], [], []
    for line in eachline(open(pdb_file))
        if startswith(line, "ATOM") && strip(line[13:16]) == "CA"
            push!(x, parse(Float64,strip(line[31:38])))
            push!(y, parse(Float64,strip(line[39:46])))
            push!(z, parse(Float64,strip(line[47:54])))
        end
        if startswith(line, "ENDMDL")
            break
        end
    end
    m = hcat(x,y,z)
    return m
end

function xyz2matrix(input_xyz)
    # Imports .xyz file format as data frame, and converts it to an N x 3 Matrix
    raw_xyz=readdlm(input_xyz)
    natoms=size(raw_xyz)[1]
    just_coords=Array{Float64}(raw_xyz[2:natoms,2:4])
    return just_coords
end

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

function kabsch(m1,m2)

    normalisedP = (translate_to_centroid(m1))
    normalisedQ = (translate_to_centroid(m2))

    xcov = cross_covariance_matrix(normalisedP,normalisedQ)
    orot = optimal_rotation_matrix(xcov)

    num_atoms = size(normalisedP)[1]
    rotated = zeros(Float64,num_atoms,3)

    for i=1:num_atoms
        rotated[i,:] = orot*normalisedP[i,:]
    end

    RMSD_value = norm(rotated-normalisedQ)
    return RMSD_value
end

function matrix2fragments(matrix, wsize=5)
    msize = size(matrix)[1]-wsize
    megam = [matrix[i:i+wsize,:] for i in 1:msize]
    return megam
end



function kabsh_matrix(megax)    
    try        
        n = size(megax)[1]  # Change this to the desired size
        matrix = zeros(Float64, n, n)

        for i = 1:n # Fill the upper triangle
            for j = i+1:n  # Ensure j >= i for the upper triangle
                matrix[i, j] = kabsch(megax[i], megax[j])
            end
        end
        return sum(matrix)/(n*n)
    catch 
        return 0
    end
end



function fxdir(dirpath, printon = false)
    i = 1
    data = []
    for (root, dirs, files) in walkdir(dirpath)
            for file in files
                    try
                        filepath = joinpath(root, file)
                        megax = matrix2fragments(pdb2matrix(filepath), 4)
                        nres = size(megax)[1]
                        fxity = kabsh_matrix(megax)
                        push!(data, [i,filepath,fxity,nres])
                        if printon
                        println("$i    $filepath   $fxity  $nres")
                        end
                    catch e
                        push!(data, [i,filepath, 0, 0])
                        if printon
                                println("$i   $filepath  0  0")
                        end
                    end
                    i+=1
            
            end
    end
    return(data)
end




pdbdir  = ARGS[1]

data = fxdir(pdbdir, true)

