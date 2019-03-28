#module Csv2mat
using MAT
using DataFrames
#export csv2mat
a = 13



"""
convert a csv table file to a .mat file for import to Matlab
csv2mat(filename, destination="log.mat"; startline=0, writeNames = true, subsample = 1, separator = ',')
"""
function csv2mat(filename, destination="log.mat"; startline=0, writeNames = true, subsample = 1, separator = ',')

    # Read the csv file into a dataframe
    println("Reading file $filename of size $(filesize(filename)/1_000_000) MB")
    df = readtable(filename, separator = separator, header=true, skipstart=startline)
    println("Preprocessing")
    (rows, cols) = size(df)
    # replace NA values with 0
    for  j = 1:cols, i = 1:rows
        ismissing(df[i,j]) && (df[i,j] = 0)
    end

    # Remove strange characters from variable names
    n = names(df)
    stripchars = [' ', '_', ':']
    for (i,symbol) in enumerate(n)
        n[i] = Symbol(strip(string(symbol),stripchars))
        writeNames && println(n[i])
    end
    names!(df,n) # Update the column names


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
    println("Write the results to $destination")
    @time matwrite(destination, d)
    printstyled("Wrote $((round(Int,ceil(size(df,1)/subsample)),written)) entries to log.mat\n", color=:green)
    printstyled("Skipped $(text) columns with text\n", color=:yellow)

    # Perform requested plotting
    # if doplot
    #     colsToPlot = cols
    #     constants = Dict{Int64,Float64}()
    #     for i = 1:cols
    #         if all(df[:,i] .== df[1,i])
    #             typeof(df[1,i]) <: Number && push!(constants,i,df[1,i])
    #             colsToPlot-= 1
    #         end
    #         if !(typeof(df[1,i])<:Number)
    #             colsToPlot-= 1
    #         end
    #     end
    #     numberOfFigures = round(Int64,ceil(colsToPlot/subplots))
    #     print("Plotting all data in $(numberOfFigures) windows ... ")
    #     for window = 0:(numberOfFigures-1)
    #         plotrows::Int64 = ceil(sqrt(minimum([colsToPlot,subplots])))
    #         plotcols::Int64 = ceil(minimum([colsToPlot,subplots])/plotrows)
    #         doplot && plot(layout=(plotrows,plotcols))
    #         # Write the results to a .mat-file
    #         plotrow = plotcol = 0
    #         skip = 0;
    #
    #         for ii in (window*subplots).+collect(1:minimum([colsToPlot,subplots]))
    #             i = ii+skip
    #             while haskey(constants,i) || !(typeof(df[1,i])<:Number)
    #                 skip += 1
    #                 i = ii+skip
    #             end
    #             try
    #                 plot!(df[:,i], subplot=plotrow*plotcols + plotcol+1, title=replace(string(n[i]),"_"," "))
    #             catch ex
    #                 print_with_color(:red,"Could not plot column: $(string(n[i]))\n")
    #             end
    #             plotcol = (plotcol + 1) % plotcols
    #             if plotcol == 0
    #                 plotrow = (plotrow + 1) % plotrows
    #             end
    #         end
    #         #winston && display(t)
    #         colsToPlot -= plotrows*plotcols
    #     end
    #     println("done")
    #     print_with_color(:yellow,"Constants (not plotted):\n")
    #     for i in keys(constants)
    #         println(string(n[i], ":\t", round(constants[i],4)))
    #     end
    # end

end
