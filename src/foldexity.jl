using ArgParse
using ProgressBars

include("fxio.jl")
include("kabsch_umeyama.jl")

#calculate fxity for a pdb file
function fxpdb(pdbpath::String, ksize::Int = 4, kmertype::String = "knn", cutoff = 0.5)
    println("Starting foldexity...")
    
    if kmertype == "knn" || kmertype == "k_nearest_neigbors" 
        fragmeter = coords2knn
        readpdb = readpdb_calpha
    elseif kmertype == "seq" || kmertype == "sequencial"
        fragmeter = coords2kmers
        if ksize < 4 
            readpdb =  readpdb_backbone #readpdb_calpha
        else 
            readpdb = readpdb_backbone
        end
    else 
        println("Warninin unknown fragmenter, k_nearest_neigbors will be used")
    end

    pdb = readpdb_backbone(pdbpath)
    if missing_residues(pdb)
        println("Warning: $pdbpath probably has missing residues, check the file")
    end

    xyzcoords = pdb2xyz(pdb)
    
    megax = fragmeter(xyzcoords, ksize)
    
    foldexity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix = fxity_kabsh(megax, cutoff)
    return foldexity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix
end


#calculate fxity for all pdb files in the directory
function fxdir(dirpath, outfile = "fxdata.tsv", ksize=4, kmertype = "seq", cutoff = 1.0, printdata = false)

    if kmertype == "knn" || kmertype == "k_nearest_neigbors" 
        fragmeter = coords2knn
        readpdb = readpdb_calpha
    elseif kmertype == "seq" || kmertype == "sequencial"
        fragmeter = coords2kmers
        if ksize < 4 
            readpdb = readpdb_backbone
        else 
            readpdb = readpdb_backbone
        end
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
            pdb = readpdb(pdbpath)
            if missing_residues(pdb)
                println("Warning: $pdbpath probably has missing residues, skipping")
                continue
            end
            xyzcoords = pdb2xyz(pdb)
#            megax = coords2knn(xyzcoords, ksize)
            megax = fragmeter(xyzcoords, ksize)
            foldexity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix  = fxity_kabsh(megax, cutoff)
            data = "$i\t$pdbpath\t$foldexity\t$aver_rmsd\t$nclusts\t$norm_nclusts\t$nfrags\n"
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
            "--kmertype", "-t"
                help = "kmer type"
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
    kmertype = args["kmertype"]
    cutoff = args["cutoff"]

    try
        if isdir(inputpath)
            fxdir(inputpath, outpath, ksize, kmertype, cutoff)  
        elseif isfile(inputpath)
            foldexity, aver_rmsd, nclusts, norm_nclusts, nfrags, all_vs_all_matrix = fxpdb(inputpath, ksize, cutoff, kmertype)  
            println("$foldexity\t$aver_rmsd\t$nclusts\t$norm_nclusts\t$nfrags")  
        else
            println("Error: The path '$inputpath' is neither a directory nor a file.")
        end
    catch err
        println("An error occurred: ", err)
        println(stacktrace(catch_backtrace()))  # Print the full stack trace
    end
end



