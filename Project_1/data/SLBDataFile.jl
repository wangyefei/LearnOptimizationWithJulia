"""
Return data from SLB 2011
Args:
    None
Returns:
    A dictionary contain all the data, keyword is the name of mineral
"""
function SLB_data_params(name::String="SLB2011.csv")::Dict
    SLB = Dict{}()
    Keys = []
    open(joinpath(@__DIR__,name)) do file
        for line in enumerate(eachline(file))
            lineArray = split(line[2],",")
            if line[1] == 1
                Keys = lineArray
            else
                Key = lineArray[1]
                temp = Dict{}()
                for number in enumerate(lineArray)
                    if number[1] <=2
                        push!(temp,Keys[number[1]] => number[2])
                    else
                        push!(temp,Keys[number[1]] => parse(Float64,number[2]))
                    end
                end
                push!(SLB, Key=>temp)
            end
        end
    end
    return SLB
end

#SLB2011compositiondata.csv
