using DelimitedFiles
using Combinatorics
using Pkg
using Missings
using DataFrames
using StatsBase
using Distributed
using Plots
using Distributions 
using IndexedTables
using StatsPlots 
using Query
using BenchmarkTools
using Dates

# Pkg.add(["DataFrames", "StatsBase", "Distributed", "Plots", "Distributions", "IndexedTables", "StatsPlots", "Query", "BenchmarkTools", "Dates"])

const letters="ABCDEFGHIJ"
const numbers="1234567890"

# Functions borrowed from Tamar Nizharadze (t.nizharadze@dkfz.de)

#=
Arguments:
    barcode::String = a barcode to be checked
Value:
    returns 1 if barcode is legitimate, 0 if a barcode can not be generated in the iFlpscLineage system
=#
function checklegitimacy(barcode::String)
    barcodeisok=0
    #check length
    if length(barcode) in [1,2,3,4,5,6,7,8,9,10]
        println("Barcode length is acceptable.")
        #check characters
        if sum(occursin.(split(barcode,""),"1234567890ABCDEFGHIJ"))==length(barcode)
            println("Barcode only contains existent iFlpscLineage segments.")
            #check duplicates - same segment in same orientation
            if length(unique(split(barcode,"")))==length(barcode)
                println("Barcode only contains each iFlpscLineage segment once.")
                    #check if segments and their inverted copies are both present
                    if  all((occursin.(split("1234567890",""),barcode) + occursin.(split("ABCDEFGHIJ",""),barcode)) .< 2)
                        println("Barcode does not contain any segments in both orientations.")
                        barcodeisok=1
                        println("Yayy, this barcode is OK!.")
                    else
                        println("Barcode Error: At least one iFlpscLineage segment is present in both orientations!")
                    end
            else
                println("Barcode Error: Barcode contains duplicated iFlpscLineage segements!")
            end
        else
            println("Barcode Error: Barcode contains at least one nonexistent iFlpscLineage segment!")
        end
    else
        println("Barcode Error: Barcode length is not acceptable!")
    end
    return barcodeisok
end
#examples
checklegitimacy("12345678") #unacceptable length
checklegitimacy("12345678X") #contains inexistent iFlpscLineage segment
#checklegitimacy("123456787") #contains a segment twice
#checklegitimacy("123456798") #contains segments in wrong positions
#checklegitimacy("12345678G") #contains both orientations of one segment
#checklegitimacy("123456789") #OK barcode
#checklegitimacy("ABC") #OK barcode

#=
Arguments:
    barcode::String             = a barcode to be excised
    FRTsites::Array{Int64,2}   = an array of two integers from 1:10 in increasing order
Value:
    if excision is possible returns an excised barcode string
    if excision is not possible returns "-"
=#
function excise!(barcode::String,FRTsites::Array{Int64,1})
    # check that FRTsites are present in the maximun recombination sites added
    if FRTsites[1] in collect(1:1:11) && FRTsites[2] in collect(1:1:11)
        # sets range between FRTsites that will get excised
        remainingsegmentindices=setdiff(1:length(barcode),FRTsites[1]:(FRTsites[2]-1))
        recombinedbarcode=barcode[remainingsegmentindices]
        return recombinedbarcode
    else
        println("Specified FRT sites are not appropriate for excision.")
        return "-"
    end
end
#examples
#barcode= "ABCDEFGHI"
#FRTsites = [2,4]
#remainingsegmentindices=setdiff(1:length(barcode),FRTsites[1]:(FRTsites[2]-1))

excise!("ABCDEFGHIJ",[2,5])
#excise!("ABCDEFGHI",[2,5])

#=
Arguments:
    barcode::String             = a barcode to be excised
    FRTsites::Array{Int64,2}   = an array of two integers from 1:10 in increasing order
Value:
    if excision is possible returns an excised barcode string
    if excision is not possible or the supplied barcode is not legitimate returns "-"
=#

# Deprecated as the legitimacy step can be performed previously to excision

function checkandexcise!(barcode::String,FRTsites::Array{Int64,1})
    if checklegitimacy(barcode)==1
        if FRTsites[1] in collect(1:1:11) && FRTsites[2] in collect(1:1:11)
            remainingsegmentindices=setdiff(1:length(barcode),FRTsites[1]:(FRTsites[2]-1))
            recombinedbarcode=barcode[remainingsegmentindices]
            return recombinedbarcode
        else
            println("Specified FRT sites are not appropriate for excision.")
            return "-"
        end
    else
        println("A non-legitimate barcode supplied!")
        return "-"
    end
end
#examples
checkandexcise!("ABCDE",[2,3])
#checkandexcise!("ABCDEFGHI",[2,5])
#checkandexcise!("AB1DE",[2,4])

#=
Arguments:
    barcode::String             = a barcode to be excised
    FRTsites::Array{Int64,2}   = an array of two integers from 1:10 in increasing order
Value:
    if inversion is possible returns an inverted barcode string
    if inversion is not possible returns "-"
=#
function invert!(barcode::String,FRTsites::Array{Int64,1})
    # check that FRTsites are present in the maximun recombination sites added
    if FRTsites[1] in collect(1:1:11) && FRTsites[2] in collect(1:1:11)
         # sets range between FRTsites that will get inverted
        segmentforinversion=barcode[FRTsites[1]:(FRTsites[2]-1)]
        invertedsegment="-"^(FRTsites[2]-FRTsites[1])
        # loop to change each character by their letter or number counterpart
        for i=FRTsites[1]:(FRTsites[2]-1)
            # numbers into letters
            if occursin(barcode[i],numbers)
                myindex=findfirst(barcode[i],numbers)
                invertedsegment=replace(invertedsegment,"-"=>letters[myindex],count=1)
            # letters into numbers
            elseif occursin(barcode[i],letters)
                myindex=findfirst(barcode[i],letters)
                invertedsegment=replace(invertedsegment,"-"=>numbers[myindex],count=1)
            end
        end
        # reverses the segment that will be inserted to emulate inversion
        invertedsegment=reverse(invertedsegment)
        # replacing inverted segment by previous segment
        recombinedbarcode=replace(barcode,segmentforinversion=>invertedsegment)
        return recombinedbarcode
    else
        println("Specified FRT sites are not appropriate for inversion.")
        return "-"
    end
end
#examples
invert!("1234567",[2,5])
#invert!("123456789",[2,4])

#=
Arguments:
    barcode::String             = a barcode to be excised
    FRTsites::Array{Int64,2}   = an array of two integers from 1:10 in increasing order
Value:
    if inversion is possible returns an inverted barcode string
    if inversion is not possible or the supplied barcode is not legitimate returns "-"
=#

# Deprecated as the legitimacy step can be performed previously to inversion

function checkandinvert!(barcode::String,FRTsites::Array{Int64,1})
    if checklegitimacy(barcode)==1
        if FRTsites[1] in collect(1:1:11) && FRTsites[2] in collect(1:1:11) && isodd(FRTsites[2]-FRTsites[1])
            segmentforinversion=barcode[FRTsites[1]:(FRTsites[2]-1)]
            invertedsegment="-"^(FRTsites[2]-FRTsites[1])
            for i=FRTsites[1]:(FRTsites[2]-1)
                if occursin(barcode[i],numbers)
                    myindex=findfirst(barcode[i],numbers)
                    invertedsegment=replace(invertedsegment,"-"=>letters[myindex],count=1)
                elseif occursin(barcode[i],letters)
                    myindex=findfirst(barcode[i],letters)
                    invertedsegment=replace(invertedsegment,"-"=>numbers[myindex],count=1)
                end
            end
            invertedsegment=reverse(invertedsegment)
            recombinedbarcode=replace(barcode,segmentforinversion=>invertedsegment)
            return recombinedbarcode
        else
            println("Specified FRT sites are not appropriate for inversion.")
            return "-"
        end
    else
        println("A non-legitimate barcode supplied!")
        return "-"
    end
end
#examples
checkandinvert!("1234567890",[2,5])
#checkandinvert!("123456789",[2,4])
#checkandinvert!("1234565",[2,5])

#=
Arguments:
    barcode::String = a barcode to be recombined once
Value:
    returns an array (a row column) of all possible barcodes reached from the supplied barcode in one recombination
=#

function recombineonce(barcode::String)
    # create all possible recombination combinations 
    FRTpairs=collect(combinations(1:(length(barcode)+1), 2))
    targetbarcodes_1 = Array{String}(undef,1,length(FRTpairs))
    targetbarcodes_2 = Array{String}(undef,1,length(FRTpairs))
    # creat an excision and inversion model for each possible recombination
    for i=1:length(FRTpairs)
            targetbarcodes_1[i]=invert!(barcode,FRTpairs[i])
            targetbarcodes_2[i]=excise!(barcode,FRTpairs[i])
        
    end
    # store them together
    targetbarcodes = [targetbarcodes_1 targetbarcodes_2]
    return targetbarcodes
end
#examples
recombineonce("1234567890")
#recombineonce("ABC")
# collect(combinations(1:(10+1), 2))

#########################################################################

# BARCODE Generation Simulation Code developed by Alvaro Regano (alvaro.regano@gmail.com)

#=
Arguments:
    barcode::String = a barcode to be recombined once
Value:
    returns one value of one possible barcode reached from the supplied barcode in one recombination
=#

function one_recombination(barcode::String)
    # create all possible recombination combinations 
    FRTpairs=collect(combinations(1:(length(barcode)+1), 2))
    # choose random pair from set
    rand_FRT_pair = rand((FRTpairs))
    recomb_barcode = "1234567890"
    # make sure that you end up with one BC left and do not delete the whole array
    if rand_FRT_pair[2]-rand_FRT_pair[1] == length(barcode)
        recomb_barcode=invert!(barcode,rand_FRT_pair)
    else
    # 50% chance of inversion or excision
        if isodd(rand((1,2)))
            recomb_barcode=invert!(barcode,rand_FRT_pair)
        else
            recomb_barcode=excise!(barcode,rand_FRT_pair)
        end
    end
    return recomb_barcode
end

# Example
# one_recombination("123")
one_recombination("1234567890")
# Errors



# Considering only one recombination extra step following 1 length BC recombination. Tamar's approach
#=
Arguments:
    barcode::String = a barcode to be recombined once
    recomb_steps::Int64 = the number of recombination steps to consider. You should reach al possible iFlpscLineage recombinations in 11 steps (JIHGFEDCBA Barcode)
Value:
    returns array where each entry is a string of the BC produced during x recombination step
=#


function barcode_simulation(barcode::String, recomb_steps::Int64)
    # set array with missing values instead of undef to allocate empty entries
    barcode_array = Array{Union{Missing,String}}(missing,1,recomb_steps)
        for i in 1:recomb_steps
            barcode = one_recombination(barcode)
            # check within loop that the length of BC is bigger than 1 if not finish by giving a 50% chance of inversion
            if (length(barcode) == 1 && i != recomb_steps)
                barcode_array[1,i] = barcode
                if isodd(rand((1,2)))
                    barcode_f = one_recombination(barcode)
                    barcode_array[1,i+1] = barcode_f
                    break
                else
                    break
                end
            else
                barcode_array[1,i] = barcode
            end
        end

    barcode_array = coalesce.(barcode_array, "-")
    return barcode_array
end

# Examples
@time begin
    barcode_simulation("123456789", 11)
end


# Letting 1 length BCs invert throughout the rest of the simulation. Maurice's approach

#=
Arguments:
    barcode::String = a barcode to be recombined once
    recomb_steps::Int64 = the number of recombination steps to consider. You should reach al possible iFlpscLineage recombinations in 11 steps (JIHGFEDCBA Barcode)
Value:
    returns array where each entry is a string of the BC produced during x recombination step
=#


function barcode_simulation(barcode::String, recomb_steps::Int64)
    # set array with missing values instead of undef to allocate empty entries
    barcode_array = Array{Union{Missing,String}}(missing,1,recomb_steps)
        for i in 1:recomb_steps
            barcode = one_recombination(barcode)
            barcode_array[1,i] = barcode
        end

    return barcode_array
end

# Examples
@time begin
    barcode_simulation("123456789", 11)
end

#=
Arguments:
    barcode::String = a barcode to be recombined once
    recomb_steps::Int64 = the number of recombination steps to consider. You should reach al possible iFlpscLineage recombinations in 11 steps (JIHGFEDCBA Barcode)
    iterations::Int64 + number of times you want to run the simulation
Value:
    returns matrix where each array is one output of the barcode_simulation function
=#

function n_barcode_simulation(barcode::String, recomb_steps::Int64, iterations::Int64)
    # calculate empty matrix first to already have the object stored
    barcode_array = Array{Union{Missing,String}}(undef,iterations,recomb_steps)
    @distributed for i = 1:iterations
        (i % 10_000_000 == 0) && println("current iteration = $(i), current time = $(now())")
        barcode_set = barcode_simulation(barcode, recomb_steps)
        barcode_array[i, :] = barcode_set
    end
    return barcode_array
end

@time begin
    barcode_array = n_barcode_simulation("1234567890", 11, 10_000)
end

# Alternating among workers for different iterations using @spawnat

@time begin
    barcode_array_1 = @spawnat :any n_barcode_simulation("1234567890", 11, 6126468)
    barcode_array_2 = @spawnat :any n_barcode_simulation("1234567890", 11, 6126468)
    barcode_array = barcode_array_1+barcode_array_2
end


writedlm("simulated_barcode_arrays.txt", barcode_array, ", ")

using SharedArrays

barcode_array

function get_simulated_pathmatrix(barcode_array::SharedArray{Union{Missing,String}})
    bc_indexes = unique(barcode_array)
    m,n = size(barcode_array)
    count_matrix=zeros(Int64, length(bc_indexes),n)
    # count_matrix = Array{Float64}(undef, length(bc_indexes),n)  
    simulated_pathmatrix = hcat(bc_indexes, count_matrix)
    # simulated_pathmatrix[:,1] = bc_indexes
    #BCindmap = Dict(bc_indexes .=> collect(1:1:length(bc_indexes)))
    BCindmap=Dict{String,Int64}(bc_indexes[i] => i for i = 1:length(bc_indexes))
    
    for i = 1:n # @distributed
        bc_matrix_col = barcode_array[:, i]
        dict = countmap(bc_matrix_col)
        id = collect(keys(dict)) # => barcodestring
        freq = collect(values(dict))
        
        for h in 1:length(id)
            target_bc = id[h]
            index = BCindmap[target_bc]
            
            if i == 1
                simulated_pathmatrix[index,2] = freq[h]/sum(freq)
            elseif i == 2
                simulated_pathmatrix[index,3] = freq[h]/sum(freq)
            elseif i == 3
                simulated_pathmatrix[index,4] = freq[h]/sum(freq)
            elseif i == 4
                simulated_pathmatrix[index,5] = freq[h]/sum(freq)
            elseif i == 5
                simulated_pathmatrix[index,6] = freq[h]/sum(freq)
            elseif i == 6
                simulated_pathmatrix[index,7] = freq[h]/sum(freq)
            elseif i == 7
                simulated_pathmatrix[index,8] = freq[h]/sum(freq)
            elseif i == 8
                simulated_pathmatrix[index,9] = freq[h]/sum(freq)
            elseif i == 9
                simulated_pathmatrix[index,10] = freq[h]/sum(freq)
            elseif i == 10
                simulated_pathmatrix[index,11] = freq[h]/sum(freq)
            elseif i == 11
                simulated_pathmatrix[index,12] = freq[h]/sum(freq)
            end

            
        end
            
 
    end
    return simulated_pathmatrix
end

@time begin
    pathmatrix = get_simulated_pathmatrix(barcode_array);
end

# Looking at parallel processing for this function, consider using paralell mapping with pmap()

@time begin
    pathmatrix = pmap(get_simulated_pathmatrix, barcode_array);
end

M = Matrix{Float64}[rand(1000,1000) for i = 1:10];
M
N = pmap(svdvals, M);
N
svdvals(M)
using LinearAlgebra


writedlm("simulated_pathmatrix.txt", pathmatrix, ", ")
pathmatrix


#########################################################################

# Visualization of simulated barcode arrays through a gel band plot

# Get the gel band approximation by establishing an in silico gel band simulation on what is the optimal 
# gel for proper sequencing of iFLpscLineage considering partial recombination


d = Poisson(0.5)

d

rand(d, 10000)

histogram(rand(d, 10000))

x = randn(10^3)
histogram(x)

function barcode_final_length_simulation(barcode::String, recomb_steps::Int64)
    # set array with missing values instead of undef to allocate empty entries
        for i in 1:recomb_steps
            barcode = one_recombination(barcode)
        end
    return length(barcode)
end

# Examples
@time begin
    barcode_final_length_simulation("123456789", 1)
end

barcode_array

#=
Arguments:
    barcode::String = a barcode to be recombined once
    recomb_rate::Int64 = the number of recombination steps per day. To this a Poisson distribution is consider to asses for stochaisticity
    iterations::Int64 + number of times you want to run the simulation
Value:
    returns vector of different barcode lengths produced in the simulation
=#

function n_final_length_barcode_simulation(barcode::String, recomb_rate::Float64, iterations::Int64)
    distr = Poisson(recomb_rate)
    recomb_distr = rand(distr, iterations)
    barcode_length_vec = Array{Int64}(undef, iterations, 1)
    # calculate empty matrix first to already have the object stored
    for i = 1:iterations
        barcode_length = barcode_final_length_simulation(barcode, recomb_distr[i])
        barcode_length_vec[i, 1] = barcode_length
    end
    return barcode_length_vec
end

@time begin
    barcode_length_vec = n_final_length_barcode_simulation("1234567890", 10.0, 1000)
end


typeof(barcode_length_vec)

histogram(barcode_length_vec, normalize=:pdf, 
    bins=0.5:1:maximum(barcode_length_vec)+0.5,
    xticks=1:1:10)

#=
Arguments:
    barcode::String = a barcode to be recombined once
    recomb_rates::Int64 = set of number of recombination steps per day. To this a Poisson distribution is consider to asses for stochaisticity
Value:
    returns array of different barcode lengths produced in the simulation
=#


function array_final_length_barcode_simulation(barcode::String, recomb_rates::Array{Float64}, iterations::Int64)
    barcode_length_array = Array{Int64}(undef, iterations, length(recomb_rates))
    for i in 1:length(recomb_rates)
        barcode_length_array[:, i] = n_final_length_barcode_simulation(barcode, recomb_rates[i], iterations)
    end
    return barcode_length_array
end


recomb_rates = [0.5 1 2 3 4 8 16]

@time begin
    barcode_length_array = array_final_length_barcode_simulation("1234567890", recomb_rates, 1000)
end

# Mapping funcion to get frequencies of BC lengths

barcode_length_array

function gel_distr_matrix(barcode_length_array::Array{Int64, 2})
    
    bc_lengths=unique(barcode_length_array)
    bc_lengths=sort(bc_lengths, rev=true)
    gel_distr = Array{Int64}(undef, length(bc_lengths), size(barcode_length_array, 2))
    for i in 1:size(barcode_length_array, 2)
    map = countmap(barcode_length_array[:,i])
        for n in 1:length(bc_lengths)
            if haskey(map, n) == true
                gel_distr[n,i] = map[n]
            else
                 gel_distr[n,i] = 0
        end 
        end
    end
    gel_distr = gel_distr[end:-1:1, :]
    
    return gel_distr
end



gel_distr = gel_distr_matrix(barcode_length_array)

begin
    gel_distr_norm = Array{Float64}(undef, size(gel_distr, 1), size(gel_distr, 2))
    for i in 1:size(gel_distr,2)
        gel_distr_norm_entry= gel_distr[:,i] ./ maximum(gel_distr[:,i])
        gel_distr_norm[:,i] = gel_distr_norm_entry
    end
    return gel_distr_norm
end
gel_distr_norm
# Plot Maurice code

using Plots
using Plots.PlotMeasures

### plot PCR gel band simulation
# define a function that returns a Plots.Shape
size(gel_distr_norm, 1)*100

plot(rand(10), xticks = (1:10, string.(recomb_rates)))

# Plot Homogeneous band gel plot

begin
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
p= plot(plot_title = "Band sizes of iFlpscLineage libraries in a gel (Homogeneous band width)", xlabel= "Recombination rate", ylabel= "Barcode length band size",grid=false, legend=false, axis=true, xticks = (1:size(gel_distr_norm,2), string.(recomb_rates)), size=((size(gel_distr_norm, 1)*100), (size(gel_distr_norm, 2)*100)), right_margin=10px, left_margin=4px);
    for i in 1:size(gel_distr_norm, 2)
        plot!(rectangle(0.7,0.2,i,10), linecolor=:white, color=:black, fillalpha=gel_distr_norm[1, i]);
        plot!(rectangle(0.7,0.2,i,9), linecolor=:white, color=:black, fillalpha=gel_distr_norm[2, i]);
        plot!(rectangle(0.7,0.2,i,8), linecolor=:white, color=:black, fillalpha=gel_distr_norm[3, i]);
        plot!(rectangle(0.7,0.2,i,7), linecolor=:white, color=:black, fillalpha=gel_distr_norm[4, i]);
        plot!(rectangle(0.7,0.2,i,6), linecolor=:white, color=:black, fillalpha=gel_distr_norm[5, i]);
        plot!(rectangle(0.7,0.2,i,5), linecolor=:white, color=:black, fillalpha=gel_distr_norm[6, i]);
        plot!(rectangle(0.7,0.2,i,4), linecolor=:white, color=:black, fillalpha=gel_distr_norm[7, i]);
        plot!(rectangle(0.7,0.2,i,3), linecolor=:white, color=:black, fillalpha=gel_distr_norm[8, i]);
        plot!(rectangle(0.7,0.2,i,2), linecolor=:white, color=:black, fillalpha=gel_distr_norm[9, i]);
        plot!(rectangle(0.7,0.2,i,1), linecolor=:white, color=:black, fillalpha=gel_distr_norm[10, i])
    end
    return p
end
pwd()
cd("Desktop/PhD/scLineageTracing/Polylox_tables_for_Alvaro/Plots")
savefig("iFlpscLineage_gel_band_Homogeneous_model.pdf")

# Plot Band gel plot considering different band size intensity

gel_band_weights = [10;9;8;7;6;5;4;3;2;1]

gel_distr_weighted = gel_distr.*gel_band_weights

begin
    gel_distr_norm = Array{Float64}(undef, size(gel_distr, 1), size(gel_distr, 2))
    for i in 1:size(gel_distr_weighted,2)
        gel_distr_norm_entry= gel_distr_weighted[:,i] ./ maximum(gel_distr_weighted[:,i])
        gel_distr_norm[:,i] = gel_distr_norm_entry
    end
    return gel_distr_norm
end
gel_distr_norm

begin
    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
    p= plot(plot_title = "Band sizes of iFlpscLineage libraries in a gel (Weighted band width)", xlabel= "Avg. recombination number", ylabel= "Barcode length band size",grid=false, legend=false, axis=true, xticks = (1:size(gel_distr_norm,2), string.(recomb_rates)), size=((size(gel_distr_norm, 1)*100), (size(gel_distr_norm, 2)*100)), right_margin=10px, left_margin=4px);
        for i in 1:size(gel_distr_norm, 2)
            plot!(rectangle(0.7,0.2,i,10), linecolor=:white, color=:black, fillalpha=gel_distr_norm[1, i]);
            plot!(rectangle(0.7,0.2,i,9), linecolor=:white, color=:black, fillalpha=gel_distr_norm[2, i]);
            plot!(rectangle(0.7,0.2,i,8), linecolor=:white, color=:black, fillalpha=gel_distr_norm[3, i]);
            plot!(rectangle(0.7,0.2,i,7), linecolor=:white, color=:black, fillalpha=gel_distr_norm[4, i]);
            plot!(rectangle(0.7,0.2,i,6), linecolor=:white, color=:black, fillalpha=gel_distr_norm[5, i]);
            plot!(rectangle(0.7,0.2,i,5), linecolor=:white, color=:black, fillalpha=gel_distr_norm[6, i]);
            plot!(rectangle(0.7,0.2,i,4), linecolor=:white, color=:black, fillalpha=gel_distr_norm[7, i]);
            plot!(rectangle(0.7,0.2,i,3), linecolor=:white, color=:black, fillalpha=gel_distr_norm[8, i]);
            plot!(rectangle(0.7,0.2,i,2), linecolor=:white, color=:black, fillalpha=gel_distr_norm[9, i]);
            plot!(rectangle(0.7,0.2,i,1), linecolor=:white, color=:black, fillalpha=gel_distr_norm[10, i])
        end
        return p
    end


savefig("iFlpscLineage_gel_band_Weighted_model.pdf")

