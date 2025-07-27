using Printf
using Distances


#readpdb from file
mutable struct PDBdata
    ndx::Vector{Int} # The sequential index of the atoms 
    index::Vector{Int} # The sequential index of the atoms in the pdb file
    atomname::Vector{String}
    resname::Vector{String}
    chain::Vector{String}
    resid::Vector{Int32} # Number of residue as written in PDB file
    x::Vector{Float32}
    y::Vector{Float32}
    z::Vector{Float32}
end


function readpdb_backbone(pdb_file::String)
    if !isfile(pdb_file)
        error("File not found: $pdb_file")
    end

    pdb = PDBdata(Int[], Int[], String[], String[], String[], Int32[], Float32[], Float32[], Float32[])
    i = 1

    open(pdb_file) do file
        for line in eachline(file)
            if startswith(line, "ATOM")
                atom_type = strip(line[13:16])
                if atom_type in ["CA", "C", "N"]
                    try
                        push!(pdb.ndx, i)
                        push!(pdb.index, parse(Int, strip(line[7:11])))
                        push!(pdb.atomname, atom_type)
                        push!(pdb.resname, strip(line[18:20]))
                        push!(pdb.chain, strip(line[22:22]))
                        push!(pdb.resid, parse(Int32, strip(line[23:26])))
                        push!(pdb.x, parse(Float32, strip(line[31:38])))
                        push!(pdb.y, parse(Float32, strip(line[39:46])))
                        push!(pdb.z, parse(Float32, strip(line[47:54])))
                        i += 1
                    catch e
                        println("Error parsing line: $line\n$e")
                    end
                end
            end
            if startswith(line, "ENDMDL")
                break
            end
        end
    end

    return pdb
end


function readpdb_calpha(pdb_file::String)
    if !isfile(pdb_file)
        error("File not found: $pdb_file")
    end

    pdb = PDBdata(Int[], Int[], String[], String[], String[], Int32[], Float32[], Float32[], Float32[])
    i = 1

    open(pdb_file) do file
        for line in eachline(file)
            if startswith(line, "ATOM")
                atom_type = strip(line[13:16])
                if atom_type == "CA"
                    try
                        push!(pdb.ndx, i)
                        push!(pdb.index, parse(Int, strip(line[7:11])))
                        push!(pdb.atomname, atom_type)
                        push!(pdb.resname, strip(line[18:20]))
                        push!(pdb.chain, strip(line[22:22]))
                        push!(pdb.resid, parse(Int32, strip(line[23:26])))
                        push!(pdb.x, parse(Float32, strip(line[31:38])))
                        push!(pdb.y, parse(Float32, strip(line[39:46])))
                        push!(pdb.z, parse(Float32, strip(line[47:54])))
                        i += 1
                    catch e
                        println("Error parsing line: $line\n$e")
                    end
                end
            end
            if startswith(line, "ENDMDL")
                break
            end
        end
    end

    return pdb
end



function cpptraj(parm, trajin,  b, e, offset, outpath, keep_log::Bool=false)
    
    if isdir(outpath)
        rm(outpath, recursive=true, force=true)
    end
    
    mkpath(outpath)

    output = joinpath(outpath, "frame")

    cpptraj_input = """
    parm $parm 
    trajin $trajin $b $e $offset
    autoimage
    rms fit @CA
    trajout $output pdb nobox multi
    go
    """

    open(`../bin/cpptraj`, "w", stdin) do io
        write(io, cpptraj_input)
    end

    for file in readdir(outpath, join=true)
        frame = split(basename(file),".")[2]
        mv(file, "$outpath/$frame.pdb")
    end

    if !keep_log
        rm("cpptraj.log", force=true)
    end



end


#write a pdb file
function writepdb(pdb, pdbpath="output.pdb")

    function align_name(name)
        name = strip(name)
        length(name) == 1 && return " $(name)  "
        length(name) == 2 && return " $(name) "
        length(name) == 3 && return " $(name)"
        return name
    end

    function align_resname(resname)
        resname = strip(resname)
        length(resname) == 1 && return "  $(resname) "
        length(resname) == 2 && return " $(resname)  "
        length(resname) == 3 && return " $(resname) "
        return resname
    end
    
    occup = 1.00
    beta = 1.00
    model = 1
    segname = "PROT"
    open(pdbpath, "w") do f 
        for i in 1:size(pdb.x)[1] 
            atomline = @sprintf(
                "%-6s%5i%1s%4s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", 
                "ATOM",                         # 1 -  6        Record name   "ATOM  "
                pdb.index[i],                   # 7 - 11        Integer       serial       Atom  serial number.
                " ",                            #
                align_name(pdb.atomname[i]),    #13 - 16        Atom          name         Atom name.
                align_resname(pdb.resname[i]),  #18 - 20        Residue name  resName      Residue name.
                pdb.chain[i],                   #22             Character     chainID      Chain identifier.
                pdb.resid[i],                   #23 - 26        Integer       resSeq       Residue sequence number.
                "    ",                         #27             AChar         iCode        Code for insertion of residues.
                pdb.x[i],                       #31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
                pdb.y[i],                       #39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
                pdb.z[i],                       #47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
                occup,                          #55 - 60        Real(6.2)     occupancy    Occupancy.
                beta,                           #61 - 66        Real(6.2)     tempFactor   Temperature  factor.
                "      ",                       #73 - 76        String        segname      Segment identifier, left-justified (not default of PDB)
                segname,                        #77 - 78        LString(2)    element      Element symbol string, right-justified.
                "  ",                           #79 - 80        LString(2)    charge       Charge  on the atom.
                ) 
            write(f, "$atomline", "\n") 
        end
    end
end


######=====pdb2matrix2pdb====######
function pdb2pdbmatrix(pdb)
    
    pdbmatrix = hcat(pdb.ndx,
                pdb.index, 
                pdb.atomname, 
                pdb.resname, 
                pdb.chain, 
                pdb.resid, 
                pdb.x, 
                pdb.y, 
                pdb.z, 
                )

    return pdbmatrix
end

function pdb2xyz_new(pdb::PDBdata)

    function add_virtual_site(A,B,C)
        # Vector from A to C
        AC = C - A
        # Vector from A to B
        AB = B - A
        # Project vector AB onto AC to get foot of perpendicular D
        t = dot(AB, AC) / dot(AC, AC)
        D = A + t * AC  # Foot of the perpendicular from B to line AC
        # Vector from D to B
        DB = B - D
        # Extend from D in the direction of B by twice the length of the height
        V = D + 4 * DB
        return V
    end

    N_ndx = pdb.atomname .== "N"
    CA_ndx = pdb.atomname .== "CA"
    C_ndx = pdb.atomname .== "C"    
    CA_indeces = pdb.index[CA_ndx]
    N_xyz = hcat(pdb.x[N_ndx], pdb.y[N_ndx], pdb.z[N_ndx])   
    CA_xyz = hcat(pdb.x[CA_ndx], pdb.y[CA_ndx], pdb.z[CA_ndx])   
    C_xyz = hcat(pdb.x[C_ndx], pdb.y[C_ndx], pdb.z[C_ndx])   

    V_xyz = add_virtual_site(N_xyz, CA_xyz, C_xyz)

    return CA_indeces, N_xyz, CA_xyz, C_xyz, V_xyz
end


function pdb2xyz(pdb)
    return hcat(pdb.x, pdb.y, pdb.z)    
end

function pdb2fasta(pdb)

    one2three = Dict('C'=> "CYS", 'D'=> "ASP", 'S'=> "SER", 'Q'=> "GLN", 'K'=> "LYS",
                 'I'=> "ILE", 'P'=> "PRO", 'T'=> "THR", 'F'=> "PHE", 'N'=> "ASN", 
                 'G'=> "GLY", 'H'=> "HIS", 'L'=> "LEU", 'R'=> "ARG", 'W'=> "TRP", 
                 'A'=> "ALA", 'V'=> "VAL", 'E'=> "GLU", 'Y'=> "TYR", 'M'=> "MET")


    three2one = Dict("CYS"=> 'C', "ASP"=> 'D', "SER"=> 'S', "GLN"=> 'Q', "LYS"=> 'K',
                 "ILE"=> 'I', "PRO"=> 'P', "THR"=> 'T', "PHE"=> 'F', "ASN"=> 'N', 
                 "GLY"=> 'G', "HIS"=> 'H', "LEU"=> 'L', "ARG"=> 'R', "TRP"=> 'W', 
                 "ALA"=> 'A', "VAL"=> 'V', "GLU"=> 'E', "TYR"=> 'Y', "MET"=> 'M')

    
    resname_sequence = pdb.resname[pdb.atomname .== "CA"]

    fasta = [three2one[RES] for RES=resname_sequence]

    return fasta 
end

function pdbmatrix2pdb(matrix)
    
    pdb = PDBdata(Int[], Int[], String[], String[], String[], Int32[], Float32[], Float32[], Float32[])

    pdb.ndx = matrix[:,1]
    pdb.index = matrix[:,2]
    pdb.atomname = matrix[:,3]
    pdb.resname = matrix[:,4]
    pdb.chain = matrix[:,5]
    pdb.resid = matrix[:,6]
    pdb.x = matrix[:,7]
    pdb.y = matrix[:,8]
    pdb.z = matrix[:,9]
    
    return pdb
end

function missing_residues(pdb)::Bool
    resid = unique(pdb.resid)
    return !all((diff(resid) .== 1))  
end

######=====pdb2matrix2pdb====######




######=====matrix=fragmentation====######


function distancematrix(xyzcoords::Matrix, min_seq_dist::Int = 0)::Matrix
    # make a distance matrix from xyz coordinates, 
    # the N sequential neighbors can be excluded by min_seq_dist
    matirx_lengnt = size(xyzcoords, 1)
 
    distmatrix = pairwise(Euclidean(), xyzcoords, dims=1)
    #matrixmax = maximum(distmatrix)
    for i in 1:matirx_lengnt-min_seq_dist
     distmatrix[i:i+min_seq_dist, i:i+min_seq_dist,] .= 1e9 #matrixmax
    end

    return distmatrix
end


function knn_form_distance_matrix(D::Matrix, k::Int=10)::Matrix
    # returns k nearest neighbors from a distance matrix
    return mapslices(x -> partialsortperm(x, 1:k), D, dims=2) 
end



function coords2knn(CA_xyz::Matrix=none, wordsize::Int=6, min_seq_dist::Int=0, N_xyz::Matrix=none, C_xyz::Matrix=none, V_xyz::Matrix=none)::Vector{Matrix}
    # returns coordinates corresponding to knn indexes
    nxyz = size(CA_xyz, 1)
    distmatrix = distancematrix(CA_xyz, min_seq_dist)
    neighbor_list_index = knn_form_distance_matrix(distmatrix, wordsize)
    #add a virtual site


    #knnfragments = [CA_xyz[push!(neighbor_list_index[i,:], i),:] for i in 1:nxyz]
    knnfragments = [vcat(N_xyz[neighbor_list_index[i,:],:], 
                        CA_xyz[neighbor_list_index[i,:],:],
                        C_xyz[neighbor_list_index[i,:],:],
                        V_xyz[neighbor_list_index[i,:],:],
                        ) for i in 1:nxyz]

    return knnfragments
end



function coords2kmers(CA_xyz::Matrix=none, wordsize=4, offset=0, N_xyz::Matrix=none, C_xyz::Matrix=none, V_xyz::Matrix=none)::Vector{Matrix}
    #split matrix into fragments
    
    backbone_length = 4 # N, CA, C, V
    offset = offset
    wsize = wordsize * backbone_length # 4 backbone residue fragment contains 16 atoms (4 atoms for each residue: N, CA, C, V).
    msize = size(CA_xyz)[1]
    
    CA_xyz = vcat([m[i, :]' for i in 1:size(CA_xyz, 1) for m in (N_xyz, CA_xyz, C_xyz, V_xyz)]...)


    kmerfragments = [CA_xyz[i:i-1+wsize,:] for i=1:msize-wsize+1]

    return kmerfragments
end


function split2kmers(seq, k::Int)
    seqlen = length(seq)
    @assert k > 0 && k <= seqlen / 2
    return [seq[i:i+k-1] for i in 1:seqlen-k+1]
end


function pdb2seqxyz(pdbpath::String)
    #returns a matrix with the sequence in the 1th column and xyz in 2-4
    pdb = readpdb_calpha(pdbpath)
    seq = pdb2fasta(pdb)
    
    xyz = pdb2xyz(pdb)
    
    return hcat(seq, xyz)
end


######=====matrix=fragmentation====######

function vtor(xyzmatrix::Matrix{T}) where {T}

    @assert size(xyzmatrix) == (4, 3)

    p1 = xyzmatrix[1,1:3]
    p2 = xyzmatrix[2,1:3]
    p3 = xyzmatrix[3,1:3]
    p4 = xyzmatrix[4,1:3]

    b1 = p2 .- p1
    b2 = p3 .- p2
    b3 = p4 .- p3

    # Normalize the vectors
    norm_b1 = norm(b1)
    norm_b2 = norm(b2)
    norm_b3 = norm(b3)

    norm_vec1 = b1 ./ norm_b1
    norm_vec2 = b2 ./ norm_b2
    norm_vec3 = b3 ./ norm_b3

    # Compute normals
    n1 = cross(norm_vec1, norm_vec2)
    n2 = cross(norm_vec2, norm_vec3)

    # Compute the dihidral
    x = dot(n1, n2)
    y = dot(cross(n1, n2), norm_vec2)

    dihidral = atan(y, x) * (180 / π)  # Converted to degrees * (180 / π)

    # Compute angles
    angle1 = 180 - acos(dot(b1, b2) / (norm_b1 * norm_b2)) * (180 / π)
    angle2 = 180 - acos(dot(b2, b3) / (norm_b2 * norm_b3)) * (180 / π)



    return [dihidral, angle1, angle2]

end

function k3angle(xyzmatrix::Matrix{T}) where {T}

    @assert size(xyzmatrix) == (3, 3)

    p1 = xyzmatrix[1,1:3]
    p2 = xyzmatrix[2,1:3]
    p3 = xyzmatrix[3,1:3]

    b1 = p2 .- p1
    b2 = p3 .- p2

    # Normalize the vectors
    norm_b1 = norm(b1)
    norm_b2 = norm(b2)

    # Compute angles
    angle1 = 180 - acos(dot(b1, b2) / (norm_b1 * norm_b2)) * (180 / π)

    return angle1

end


function euclid_dist(p1::Vector, p2::Vector)
    return sqrt(sum((p2-p1).^2))
end
