using Printf
using Distances

one2three = Dict('C'=> "CYS", 'D'=> "ASP", 'S'=> "SER", 'Q'=> "GLN", 'K'=> "LYS",
                 'I'=> "ILE", 'P'=> "PRO", 'T'=> "THR", 'F'=> "PHE", 'N'=> "ASN", 
                 'G'=> "GLY", 'H'=> "HIS", 'L'=> "LEU", 'R'=> "ARG", 'W'=> "TRP", 
                 'A'=> "ALA", 'V'=> "VAL", 'E'=> "GLU", 'Y'=> "TYR", 'M'=> "MET")


three2one = Dict("CYS"=> 'C', "ASP"=> 'D', "SER"=> 'S', "GLN"=> 'Q', "LYS"=> 'K',
                 "ILE"=> 'I', "PRO"=> 'P', "THR"=> 'T', "PHE"=> 'F', "ASN"=> 'N', 
                 "GLY"=> 'G', "HIS"=> 'H', "LEU"=> 'L', "ARG"=> 'R', "TRP"=> 'W', 
                 "ALA"=> 'A', "VAL"=> 'V', "GLU"=> 'E', "TYR"=> 'Y', "MET"=> 'M')

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

    pdb = PDBdata(Int[], Int[], String[], String[], String[], Int32[], Float32[], Float32[], Float32[] )
    i = 1
    for line in eachline(open(pdb_file))
        if startswith(line, "ATOM") 
            if strip(line[13:16]) == "CA" || strip(line[13:16]) == "C" || strip(line[13:16]) == "N"
                push!(pdb.ndx, i) #count from 1 to ...
                push!(pdb.index, parse(Int,strip(line[7:11]))) #index in original pdb
                push!(pdb.atomname, strip(line[13:16]))
                push!(pdb.resname, strip(line[17:21]))
                push!(pdb.chain, strip(line[22:22]))
                push!(pdb.resid, parse(Int,strip(line[23:26])))
                push!(pdb.x, parse(Float32, strip(line[31:38])))
                push!(pdb.y, parse(Float32, strip(line[39:46])))
                push!(pdb.z, parse(Float32, strip(line[47:54])))
                i+=1
            end
        end
        if startswith(line, "ENDMDL")
            break
        end
    end
    return pdb
end

function readpdb_calpha(pdb_file::String)

    pdb = PDBdata(Int[], Int[], String[], String[], String[], Int32[], Float32[], Float32[], Float32[] )
    i = 1
    for line in eachline(open(pdb_file))
        if startswith(line, "ATOM") 
            if strip(line[13:16]) == "CA" 
                push!(pdb.ndx, i) #count from 1 to ...
                push!(pdb.index, parse(Int,strip(line[7:11]))) #index in original pdb
                push!(pdb.atomname, strip(line[13:16]))
                push!(pdb.resname, strip(line[17:21]))
                push!(pdb.chain, strip(line[22:22]))
                push!(pdb.resid, parse(Int,strip(line[23:26])))
                push!(pdb.x, parse(Float32, strip(line[31:38])))
                push!(pdb.y, parse(Float32, strip(line[39:46])))
                push!(pdb.z, parse(Float32, strip(line[47:54])))
                i+=1
            end
        end
        if startswith(line, "ENDMDL")
            break
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

function pdb2xyz(pdb)
    return hcat(pdb.x, pdb.y, pdb.z)    
end

function pdb2fasta(pdb)
    
    resname_sequence = pdb.resname[pdb.atomname .== "CA"]

    fasta = join([three2one[RES] for RES=resname_sequence])

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

function missing_residues(pdb)
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


function matrix_knn(D::Matrix, k::Int=10)::Matrix
    # returns k nearest neighbors from a distance matrix
    return mapslices(x -> partialsortperm(x, 1:k), D, dims=2) 
end


function coords2kmers(matrix, wordsize=4) 
    #split matrix into fragments

    backbone_length = 3
    wsize = wordsize * backbone_length # 4 backbone residue fragment contains 12 atoms (3 atoms for each residue: N, CA, C).
    msize = size(matrix)[1]
    
    fragmentsmatrix = [matrix[i:i-1+wsize,:] for i=1:3:msize-wsize+1]

    return fragmentsmatrix
end


function coords2knn(xyzcoords::Matrix, wordsize::Int=6, min_seq_dist::Int=0)::Vector{Matrix}
    # returns coordinates corresponding to knn indexes
    nxyz = size(xyzcoords, 1)
    distmatrix = distancematrix(xyzcoords, min_seq_dist)
    neighbor_list_index = matrix_knn(distmatrix, wordsize)
    
    knnfragments = [xyzcoords[push!(neighbor_list_index[i,:], i),:] for i in 1:nxyz]

    return knnfragments
end


######=====matrix=fragmentation====######

