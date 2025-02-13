#foldseek structureto3didescriptor 1hv4.pdb 1hv4_3di.dump
using DataFrames
using CSV

function readfasta(input::String)
    df = DataFrame([[],[]], ["id", "seq"])

    id = ""
    sequence = ""

    for line in eachline(open(input))
        if startswith(line, ">") # is header

            if !isempty(sequence)
                push!(df, [id, sequence])
                sequence=""
            end

            id = string(split(line, ">")[2])

        else    # header
            sequence  = sequence * line           
        end
        
    end
    push!(df, [id, sequence])
    return df
end


function structure2fs3di(input::String, output::String = "tmp", keep_3di::Bool=false)
    #redirect_stdout(devnull)
    run(`../bin/foldseek structureto3didescriptor -v 0 $input $output`)
    df = CSV.read(output, DataFrame, delim="\t", header=["id","seqaa", "seq3di", "coords"])
    if !keep_3di
        rm.([output, "$output.dbtype"], force=true)
    end

    return df
end


function structure2rsmu(input::String, output::String = "tmp", keep_mu::Bool=false)
    #redirect_stdout(devnull)
    run(`../bin/reseek -convert2mu $input -fasta $output  `)
    df = readfasta(output)
    if !keep_mu
        rm.([output, "$output.dbtype"], force=true)
    end

    return df
end

function split2kmers(seq, k)
    return [seq[i:i+k-1] for i in 1:length(seq)-k]
end

function entropy_shannon(kmers, k=1)

    if k > 1
        kmers = split2kmers(kmers, k)
    end

    seqlen = length(kmers)

    counts = Dict{Any, Int}()
    for kmer in kmers
        counts[kmer] = get(counts, kmer, 0) + 1
    end
    
    probabilities = [count / seqlen for count in values(counts)]    
    entropy = -sum([p * log2(p) for p in probabilities])
    #norm_entropy = entropy / length(counts)
    return  entropy #, norm_entropy
end

function entropy_profile(seq, k=12)

    if k > 1
        k_mers = split2kmers(seq, k)
    end

    return entropy_shannon.([kmer for kmer in k_mers])

end





if abspath(PROGRAM_FILE) == @__FILE__

    in = ARGS[1]
    out = ARGS[2]

    df = structure2fs3di(in, out)
    println(df)

end
