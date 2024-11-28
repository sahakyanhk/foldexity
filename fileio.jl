using Printf

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

function readpdb(pdb_file::String)

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


#write a pdb file
function writepdb(pdb, pdbpath="backbone.pdb")

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

function pdb2matrix(pdb)
    return hcat(pdb.x, pdb.y, pdb.z)    
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

#split matrix into fragments
function matrix2fragments(matrix, wordsize=4) 
    
    backbone_length = 3
    wsize = wordsize * backbone_length # 4 backbone residue fragment contains 12 atoms (3 atoms for each residue: N, CA, C).
    msize = size(matrix)[1]
    
    fragmentsmatrix = [matrix[i:i-1+wsize,:] for i=1:3:msize-wsize+1]

    return fragmentsmatrix
end


######=====pdb2matrix2pdb====######

