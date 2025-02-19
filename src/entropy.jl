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


function shuffle_string(s) 
    return String(shuffle(collect(s)))
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
    run(`../bin/reseek -convert2mu $input -fasta $output`)
    df = readfasta(output)
    if !keep_mu
        rm.([output, "$output.dbtype"], force=true)
    end

    return df
end


function structure2dssp(input::String, output::String = "tmp", keep_mu::Bool=false)
    
    path = input
    df = DataFrame([[],[]], ["id", "ss"])

    for pdb in readdir(path)
        pdbname = string(split(pdb, ".")[1])
        try
            dssp_cmd = pipeline(`../bin/dssp $path/$pdb`, `sed -n '/#/,$p'`, `awk '{print substr($0, 17,1)}'`, `tr ' ' 'C'`, `tr -d '\n'`)
            ss = read(dssp_cmd, String)
            if true #check if valid string
                push!(df, [pdbname, ss])
            end
        catch 
            println("warrning")
        end    
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


function lz76(sequence, k)

    if k > 1
        sequence = split2kmers(sequence, k)
    end

    sub_strings = Set()
    n = length(sequence)

    ind = 1
    inc = 0
    while true
        if ind + inc > n
            break
        end
        sub_str = sequence[ind : ind + inc]
        if sub_str in sub_strings
            inc += 1
        else
            push!(sub_strings, sub_str)
            ind += (inc+1)
            inc = 0
        end
    end

    return length(sub_strings)
end
    



if abspath(PROGRAM_FILE) == @__FILE__

    in = ARGS[1]
    out = ARGS[2]

    df = structure2fs3di(in, out)
    println(df)

end
