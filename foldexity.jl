include("fileio.jl")
include("kabsch_umeyama.jl")

#calculate fxity for a pdb file
function fxpdb(pdbpath)
    println("Starting foldexity ...")
    pdb = readpdb(pdbpath)
    coordmatrix = pdb2matrix(pdb)
    megax = matrix2fragments(coordmatrix, 4)
    nfrags = length(megax)
    fxity, m = fxity_kabsh(megax)
    return fxity, nfrags
end

#calculate fxity for all pdb files in the directory
function fxdir(dirpath, outfile = "fxdata.tsv", fsize=12, printdata = false)
    
    pdbpaths = []
    for (root, _, files) in walkdir(dirpath)
        for file in files
            push!(pdbpaths, joinpath(root, file))
        end
    end
    println("$(length(pdbpaths)) files are collected \nStarting foldexity ...")
    
    i = 1
    data_collector = []

    if isdir(outfile)
        rm(outfile)
    end

    # Thread-safe writing to the output file
    open(outfile, "w") do f 
        write(f, "ndx\tpdbpath\tfxity\tnclusts\taver_rmsd\tnres\n" )
    end

    #start loop with muptithreading
    Threads.@threads for pdbpath in pdbpaths
        try     
            pdb = readpdb(pdbpath)
            coordmatrix = pdb2matrix(pdb)
            megax = matrix2fragments(coordmatrix, fsize)
            nres = length(megax)
            aver_rmsd, nclusts, fxity, nres, matrix = fxity_kabsh(megax)
            data = "$i\t$pdbpath\t$fxity\t$nclusts\t$aver_rmsd\t$nres\n"
            push!(data_collector, data)
        catch 
            push!(data_collector, "$i\t$pdbpath\t\t\t\t\n")

            println("Warrning: no data for $pdbpath")
        end
        i+=1
        if i % 50 == 0 

            output = join(data_collector)

            open(outfile, "a") do f
                write(f, output)
            end

            if printdata
                print(output)
            end

            data_collector = []

        end
    end
    #write the remaining incompleate chunck
    output = join(data_collector)
    open(outfile, "a") do f
        write(f, output)
    end

    if printdata
        print(output)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__

    inputpath = ARGS[1]
    outpath = ARGS[2]
    fsize = parse(Int64, ARGS[3])

    if isdir(inputpath)
        data = fxdir(inputpath, outpath, fsize)  
    elseif isfile(inputpath)
        fxity, nfrags = fxpdb(inputpath)  
        println("$inputpath        $fxity       $nfrags")  
    else
        println("Error: The path '$inputpath' is neither a directory nor a file.")
    end

end

