using ArgParse
using ProgressBars

include("fxio.jl")
include("kabsch_umeyama.jl")

#calculate fxity for a pdb file
function fxpdb(pdbpath::String, ksize::Int = 4, ktype::String = "knn", cutoff::Float32 = 1.0)
    println("Starting foldexity...")
    pdb = readpdb_backbone(pdbpath)
    if missing_residues(pdb)
        println("Warning: $pdbpath probably has missing residues, check the file")
    end
    xyzcoords = pdb2xyz(pdb)
    
    megax = coords2knnfragments(xyzcoords, ksize)
    #megax = coords2fragments(xyzcoords, wsize)
    
    fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix = fxity_kabsh(megax, cutoff)
    return fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix
end

#calculate fxity for all pdb files in the directory
function fxdir(dirpath, outfile = "fxdata.tsv", ksize=4, ktype = "knn", cutoff = 1.0, printdata = false)

    if ktype == "knn" || ktype == "k_nearest_neigbors" 
        fragmeter = coords2knnfragments
    elseif ktype == "seq" || ktype == "sequencial"
        fragmeter = coords2fragments
    else 
        println("Warninin unknown fragmenter, k_nearest_neigbors will be used")
    end
        
    
    pdbpaths = []
    for (root, _, files) in walkdir(dirpath)
        for file in files
            push!(pdbpaths, joinpath(root, file))
        end
    end
    println("$(length(pdbpaths)) files are collected \nStarting foldexity...")
    
    i = 0
    data_collector = []

    if isdir(outfile)
        rm(outfile)
    end

    # Thread-safe writing to the output file
    open(outfile, "w") do f 
        write(f, "ndx\tpdbpath\tfxity\taver_rmsd\tnclusts\tnorm_nclusts\tnfrags\n")
    end

    #start loop with muptithreading
    Threads.@threads for pdbpath in ProgressBar(pdbpaths)
        try     
            pdb = readpdb_backbone(pdbpath)
            if missing_residues(pdb)
                println("Warning: $pdbpath probably has missing residues, skipping")
                continue
            end
            xyzcoords = pdb2xyz(pdb)
#            megax = coords2knnfragments(xyzcoords, ksize)
            megax = coords2fragments(xyzcoords, ksize)
            fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix  = fxity_kabsh(megax, cutoff)
            data = "$i\t$pdbpath\t$fxity\t$aver_rmsd\t$nclusts\t$norm_nclusts\t$nfrags\n"
            push!(data_collector, data)
        catch 
            push!(data_collector, "$i\t$pdbpath\t\t\t\t\t\n")
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

    function parse_commandline()
        s = ArgParseSettings()
        @add_arg_table s begin
            "--inputpath", "-i"
                help = "an option with an argument"
                required = true
            "--outpath", "-o"
                help = "another option with an argument"
                required = true
            "--ksize", "-k"
                help = "kmer size"
                arg_type = Int
                default = 4
            "--ktype", "-t"
                help = "kmer type"
                arg_type = Strung
                default = "knn"
            "--cutoff", "-c"
                arg_type = Float64
                default = 1.0
                help = "rmsd cutoff for clustering"
        end
        return parse_args(ARGS, s)
    end

    args = parse_commandline()

    inputpath = args["inputpath"]
    outpath = args["outpath"]
    ksize = args["ksize"]
    ktype = args["ktype"]
    cutoff = args["cutoff"]

    # inputpath = ARGS[1]
    # outpath = ARGS[2]
    # ksize = parse(Int64, ARGS[3])
    # cutoff = parse(Float64, ARGS[4])
    try
        if isdir(inputpath)
            fxdir(inputpath, outpath, ksize, type, cutoff)  
        elseif isfile(inputpath)
            fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, all_vs_all_matrix = fxpdb(inputpath, ksize, cutoff, ktype)  
            println("$fxity\t$aver_rmsd\t$nclusts\t$norm_nclusts\t$nfrags")  
        else
            println("Error: The path '$inputpath' is neither a directory nor a file.")
        end
    catch err
        println("An error occurred: ", err)
        println(stacktrace(catch_backtrace()))  # Print the full stack trace
    end
end



