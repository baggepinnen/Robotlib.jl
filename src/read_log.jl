# This file comtains some useful functions to import logfiles for use in either julia or matlab.
# In julia, the data object from orcalog2mat can be used together with the function getdata
# In matlab, the function csv2mat can be used to convert the textfile to a .mat-file which is very fast to read

using DelimitedFiles

"""
    csv2dict(filename, destination="log.mat"; startline=0, writeNames = true, subsample = 1, separator = ',')

convert a csv table file to a dict
"""
function csv2dict(filename, destination="log.mat"; startline=0, writeNames = true, subsample = 1, separator = ',')

    # Read the csv file into a matrix
    println("Reading file $filename of size $(filesize(filename)/1_000_000) MB")
    df,n = readdlm(filename, separator, header=true, skipstart=startline)
    println("Preprocessing")
    (rows, cols) = size(df)
    # replace NA values with 0
    for  j = 1:cols, i = 1:rows
        ismissing(df[i,j]) && (df[i,j] = 0)
    end

    # Remove strange characters from variable names
    replacechars = [' ', ':', '.', '[', ']']
    stripchars = ['_']
    for (i,symbol) in enumerate(n)
        n[i] = (replace(string(symbol),replacechars => "_"))
        n[i] = (strip(n[i],stripchars))
        writeNames && println(n[i])
    end


    # Write the results to a .mat-file
    written = text = 0
    println("Creating save dict")
    d = Dict{String, Vector{<:Real}}()
    for i in 1:cols
        if !isa(df[1,i],Number)
            text += 1
            continue
        end
        written += 1
        a = replace(convert(Array,df[1:subsample:end,i]), missing=>0)
        d[string(n[i])] = a
    end
    return d
end


"""`findData(pattern, data::Dict)`"""
function findData(pattern, data::AbstractDict)
    inds = falses(length(keys(data)))
    pattern isa Regex || (pattern = Regex(pattern))
    for (i,key) in enumerate(keys(data))
        if match(pattern, key) != nothing
            inds[i] = true
        end
    end
    return inds
end


"""
    `array = getdata(pattern, data::Dict, ds=1; removeNaN = false)`

pattern is a string representation of a regex used to match the data
e.g. \"posRawAbs\"
"""
function getdata(pattern, data::AbstractDict, ds=1; removeNaN = false)
    inds = findData(pattern, data)
    if sum(inds) < 1
        @warn("No data found using $pattern, available headers are")
        foreach(println, keys(data))
        return
    end
    names = collect(keys(data))[inds]
    function lt(x,y) # String comparison function, this soring places 2 in front of 10
        if length(x) < length(y)
            return true
        elseif length(x) > length(y)
            return false
        else
            return x < y
        end
    end
    sort!(names,lt=lt)
    retmat = Matrix{Float64}(undef, ceil(Int,length(data[names[1]])/ds),length(names))
    for (j,key) in enumerate(names)
        for (i,k)= enumerate(1:ds:length(data[key]))
            retmat[i,j] = data[key][k]
        end
    end

    if removeNaN
        retmat[any(isnan(retmat),2),:] = []
    end
    return retmat
end
