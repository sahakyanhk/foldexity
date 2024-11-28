include("fileio.jl")
include("kabsch_umeyama.jl")

#calculate fxity for all pdb files in the directory
function fxdir(dirpath, printdata = false)
    i = 1
    data = []
    for (root, dirs, files) in walkdir(dirpath)
            for file in files
                    try     
                            pdbpath = joinpath(root, file)
                            pdb = readpdb(pdbpath)
                            coordmatrix = hcat(pdb.x, pdb.y, pdb.z)
                            megax = matrix2fragments(coordmatrix, 4)
                            nres = size(megax)[1]
                            fxity = kabsh_matrix(megax)
                            push!(data, [i,pdbpath,fxity,nres])
                            if printdata
                            println("$i     $pdbpath        $fxity       $nres")
                            end
                    catch 
                            push!(data, [i,pdbpath,0,0])
                            if printdata
                                   println("$i     $pdbpath       NA")
                            end
                    end
                    i+=1
            
            end
    end
    return(data) #ToDo open and write data to a file
end

#calculate fxity a pdb file
function fxpdb(pdbpath)
    pdb = readpdb(pdbpath)
    coordmatrix = hcat(pdb.x, pdb.y, pdb.z)
    megax = matrix2fragments(coordmatrix, 4)
    nfrags = size(megax)[1]
    fxity = kabsh_matrix(megax)
    return fxity, nfrags
end

if abspath(PROGRAM_FILE) == @__FILE__

    inputpath = ARGS[1]

    if isdir(inputpath)
        data = fxdir(inputpath, true)  
    elseif isfile(inputpath)
        fxity, nfrags = fxpdb(inputpath)  
        println("$inputpath        $fxity       $nfrags")  
    else
        println("Error: The path '$inputpath' is neither a directory nor a file.")
    end

end