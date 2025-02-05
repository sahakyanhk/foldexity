#foldseek structureto3didescriptor 1hv4.pdb 1hv4_3di.dump
using DataFrames
using CSV

function structure2fs3di(input, output)
    run(`bin/foldseek structureto3didescriptor $input $output`)
    df = CSV.read(output, DataFrame, delim="\t", header=["id","seqaa", "seq3di", "coords"])
    return df
end


#https://discourse.julialang.org/t/how-to-count-all-unique-character-frequency-in-a-string/19342/3
function countchars(s)  
    res = Dict{Char, Int}()
    for c in s
        res[c] = get(res, c, 0) + 1
    end
    return res
end

# https://www.reddit.com/r/learnpython/comments/g1sdkh/python_programming_challenge_calculating_shannon/

function shannon(seq)
    frequencies = collect(values(countchars(seq)))/length(seq)
    return -sum([f * log2(f) for f in frequencies])
end


if abspath(PROGRAM_FILE) == @__FILE__

    in = ARGS[1]
    out = ARGS[2]

    df = structure2fs3di(in, out)
    println(df)

end
