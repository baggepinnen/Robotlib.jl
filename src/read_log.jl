# This file comtains some useful functions to import logfiles for use in either julia or matlab.
# In julia, the data object from orcalog2mat can be used together with the function getData
# In matlab, the function csv2mat can be used to convert the textfile to a .mat-file which is very fast to read
using MAT

"""`findData(pattern::AbstractString, data)`"""
function findData(pattern::AbstractString, data)
    indexes = falses(length(keys(data)))
    for (i,key) in enumerate(keys(data))
        if match(Regex(pattern), key) != nothing
            indexes[i] = true
        end
    end
    return indexes
end

"""
`data = orcalog2mat(pathopen, pathsave)`
`pathopen` is the full path to the file
`pathsave` is the full path to the saved file, should typically end with .mat
the returned `data` object is a `Dict{ByteString,Any}`
"""
function orcalog2mat(pathopen, pathsave; separator = ',')
    csv2mat(pathopen, pathsave, writeNames = false, separator = separator)
    data = MAT.matread(pathsave)
    return data
end


"""
`data = readmat(pathopen)`
`pathopen` is the full path to the file
the returned `data` object is a `Dict{ByteString,Any}`
"""
function readmat(pathopen)
    data = MAT.matread(pathopen)
    return data
end

"""
`array = getData(pattern, data, ds=1; removeNaN = false)`
pattern is a string representation of a regex used to match the data
e.g. \"posRawAbs\"
"""
function getData(pattern, data, ds=1; removeNaN = false)
    indexes = findData(pattern, data)
    if sum(indexes) < 1
        @warn("No data found using the string $pattern")
        return
    end
    names = collect(keys(data))[indexes]
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
    retmat = Matrix{Float64}(ceil(Int,length(data[names[1]])/ds),length(names))
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
