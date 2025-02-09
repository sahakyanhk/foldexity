using ArgParse
using ProgressBars

include("fileio.jl")
include("kabsch_umeyama.jl")

#calculate fxity for a pdb file
function fxpdb(pdbpath, wsize = 4, cutoff = 1.0)
    println("Starting foldexity...")
    pdb = readpdb(pdbpath)
    if missing_residues(pdb)
        println("Warning: $pdbpath probably has missing residues, check the file")
    end
    coordmatrix = pdb2matrix(pdb)
    megax = matrix2fragments(coordmatrix, wsize)
    
    fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix = fxity_kabsh(megax, cutoff)
    return fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix
end

#calculate fxity for all pdb files in the directory
function fxdir(dirpath, outfile = "fxdata.tsv", ksize=4, cutoff = 1.0, printdata = false)
    
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
        write(f, "ndx\tpdbpath\tfxity\tnorm_fxity\taver_rmsd\tnclusts\tnorm_nclusts\tnfrags\n")
    end

    #start loop with muptithreading
    Threads.@threads for pdbpath in ProgressBar(pdbpaths)
        try     
            pdb = readpdb(pdbpath)
            if missing_residues(pdb)
                println("Warning: $pdbpath probably has missing residues, skipping")
                continue
            end
            coordmatrix = pdb2matrix(pdb)
            megax = matrix2fragments(coordmatrix, ksize)
            fxity, norm_fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix  = fxity_kabsh(megax, cutoff)
            data = "$i\t$pdbpath\t$fxity\t$norm_fxity\t$aver_rmsd\t$nclusts\t$norm_nclusts\t$nfrags\n"
            push!(data_collector, data)
        catch 
            push!(data_collector, "$i\t$pdbpath\t\t\t\t\t\t\n")
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
    cutoff = args["cutoff"]

    # inputpath = ARGS[1]
    # outpath = ARGS[2]
    # ksize = parse(Int64, ARGS[3])
    # cutoff = parse(Float64, ARGS[4])
    try
        if isdir(inputpath)
            fxdir(inputpath, outpath, ksize, cutoff)  
        elseif isfile(inputpath)
            fxity, norm_fxity, aver_rmsd, nclusts, norm_nclusts, nfrags, matrix = fxpdb(inputpath, ksize, cutoff)  
            println("$fxity\n$orm_fxity\t$aver_rmsd\t$nclusts\t$norm_nclusts\t$nfrags")  
        else
            println("Error: The path '$inputpath' is neither a directory nor a file.")
        end
    catch err
        println("An error occurred: ", err)
        println(stacktrace(catch_backtrace()))  # Print the full stack trace
    end
end



