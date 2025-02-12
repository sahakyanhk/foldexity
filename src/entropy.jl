#foldseek structureto3didescriptor 1hv4.pdb 1hv4_3di.dump
using DataFrames
using CSV

function structure2fs3di(input::String, output::String = "tmp", keep_3di::Bool=false)
    #redirect_stdout(devnull)
    run(`bin/foldseek structureto3didescriptor -v 0 $input $output`)
    df = CSV.read(output, DataFrame, delim="\t", header=["id","seqaa", "seq3di", "coords"])
    if !keep_3di
        rm.([output, "$output.dbtype"], force=true)
    end

    return df
end

function structure2rsmu(input::String, output::String = "tmp", keep_mu::Bool=false)
    #redirect_stdout(devnull)
    run(`bin/reseek -convert2mu $input -fasta $output  `)
    df = CSV.read($output.mu, DataFrame, delim="\t", header=[```not finished```])
    if !keep_mu
        rm.([output, "$output.dbtype"], force=true)
    end

    return df
end


function shannon(seq, k=1)

    if k > 1
        seq = [seq[i:i+k] for i in 1:length(seq)-k]
    end

    counts = Dict{Any, Int}()
    for aa in seq
        counts[aa] = get(counts, aa, 0) + 1
    end

    seqlen = length(seq)
#    seqlen = 20^k
    
    probabilities = [count / seqlen for count in values(counts)]    
    entropy = -sum([p * log2(p) for p in probabilities])
    norm_entropy = entropy / length(counts)
    return  entropy, norm_entropy
end


if abspath(PROGRAM_FILE) == @__FILE__

    in = ARGS[1]
    out = ARGS[2]

    df = structure2fs3di(in, out)
    println(df)

end
