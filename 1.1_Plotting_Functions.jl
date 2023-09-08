#This code is in notebook form to make it more readable.
#The same code is available as a Julia script (.jl) file for importing it into other scripts.

#Some of this code requires the Julia package Interact along with the Jupyter extension WebIO. 

using Colors
using ProgressMeter
using Plots
using Interact

# function that plots a graph and indicates the node with maximum closeness centrality (slack bus) by color
function example_plot(n, n0, p, q, r; size=2, reps=100)

    grid = create_grid(n, n0, p, q, r)

    #println("exemplary graph:")

    types = ones(Int,n)
    slack_index = findmax(closeness_centrality(grid.graph))[2]
    types[slack_index] = 2
    colors = [colorant"tomato2", colorant"royalblue"]
    
    solar_map = set_prosumers(n,0,0,slack_index)
    flowtest_p_timeseries = power_randomdraws(solar_map.old, solar_dict_default, consumer_dict_default,
                                              0, reps, slack_index)   
    caps = ceil.(1.75*max_flows(grid,flowtest_p_timeseries, slack_index), digits=16)
    
    
    widths = size*caps/maximum(caps)
    
    loc_x = vertices_loc(grid,1)
    loc_y = vertices_loc(grid,2)
    gplot(grid.graph, loc_x, loc_y, 
          nodefillc=colors[types],
          nodesize=types.+size,# NODESIZE=1,
          edgelinewidth=widths, EDGELINEWIDTH = size, edgestrokec=colorant"black")
end

# function that computes and prints statistics about graphs created by the power grid genaration algorithm
function grid_specs(n, n0, p, q, r; reps=1000)
    
    slacks, lengths, edges, dists, degrees, lcc, gcc, cycles = [], [], [], [], [], [], [], []
    
    for i in 1:reps
        grid = create_grid(n, n0, p, q, r)
        slack_index = findmax(closeness_centrality(grid.graph))[2]
        append!(edges, ne(grid))
        append!(degrees, mean(degree(grid.graph)) )
        for j in nv(grid)
            append!(dists, gdistances(grid, j))
        end
        append!(slacks,gdistances(grid, slack_index))
        append!(lcc, mean(local_clustering_coefficient(grid.graph)) )
        append!(gcc, global_clustering_coefficient(grid.graph) )
        append!(cycles, (simplecyclescount(DiGraph(grid.graph))-ne(grid))/2 )
        append!(lengths, sum(line_lengths(grid)))
    end
    
    println("average graph properties:")
    println()
    println("edges: ",mean(edges)," +- ",std(edges))
    println("geodesic distances: ",mean(dists)," +- ",std(dists))
    println("geodesic distance to slack bus: ",mean(slacks)," +- ",std(slacks))
    println("degree: ", mean(degrees)," +- ",std(degrees))
    println("mean local clustering coefficient: ", mean(lcc)," +- ",std(lcc))
    println("global clustering coefficient: ",mean(gcc)," +- ",std(gcc))
    println("fundamental cycles: ", mean(cycles)," +- ",std(cycles))
    println("total line length: ", mean(lengths)," +- ",std(lengths))
end

# function that concatenates all arrays/lists within a distionary to a single list
function dict_hist(dict)
    
    list = []
    for key in keys(dict)
        for e in dict[key]
            append!(list, e)
        end
    end

    return list
    
end

# struct with two fields
struct combo_struct
    list
    ratio
end

# function to combine dictionaries of consumption + supply time series, and save the combination ratio
function combo_hist(c_dict, s_dict; ratio=1)
    
    list = []
    for c_key in keys(c_dict)
        for s_key in keys(s_dict)
            append!(list, ratio*s_dict[s_key]-c_dict[c_key])
        end
    end

    return combo_struct(sort(list),ratio)
    
end

# function that calculates the maximum span of adverse influence parameters that a type of power grid can cope with
# depending on the number of different grid instances and days of simulation
function rangefinder_test(path; grids_range=1:2, days_range=1:2, reps=10, safety=1.5, test_runs=100)
   
    locs = 1
    ftr = test_runs
    rtr = ftr
    
    ratio_span = Dict()
    P2_span = Dict()
    
    #better use one for loop with parameter tuples
    #@distributed 
    for grids in grids_range
        
        #@distributed 
        for days in days_range
            
            ratio_list = []
            P2_list = []
            
            #@distributed 
            for rep in 1:reps
                append!(ratio_list, ratio_finder("none", 0, 0, 0, 0, 0.99, 1., 10.,
                                                 grids=grids, locs=locs, safety=safety, ftr=ftr, rtr=rtr, days=days))
                append!(P2_list, P2_finder("none", 0, 0, 0, 0, 0.99, 10, 99,
                                           grids=grids, locs=locs, safety=safety, ftr=ftr, rtr=rtr, days=days))
            end
            
            key = (grids, days)
            
            ratio_span[key] = maximum(ratio_list)-minimum(ratio_list)
            P2_span[key] = maximum(P2_list)-minimum(P2_list)
            
        end
    end
    
    xs = collect(grids_range)
    ys = collect(days_range)
    
    function rs(x,y)
        return ratio_span[(x,y)]
    end
    function ps(x,y)
        return P2_span[(x,y)]
    end
    
    title_rest =  " span at $(reps) reps \n for safety=$(safety)
                    and flowtest_runs=$(ftr)"
    filename_rest =  "_span_at_$(reps)_reps_for_safety=$(safety)
                      _and_flowtest_runs=$(ftr)"
    
    h1 = heatmap(xs, ys, rs, title= "ratio"*title_rest, xlabel="grids", ylabel="days")
    h2 = heatmap(xs, ys, ps, title= "P2"*title_rest, xlabel="grids", ylabel="days")
    plot(h1,h2, size=(1000,300))
    savefig(path*"/ratio_and_P2"*filename_rest*".png")
    
end

# improved version of the above function that calculates the maximum span of adverse influence parameters that a type of power grid can cope with
# depending on the number of different grid instances, days of simulation, and safety factor/ number of test runs used to determine the initial grid line capacities
function rangefinder_test_2(;path=plot_path, reps=10,
                            grids_range=1:2, days_range=1:2, safety_range=1.5:0.5:2, test_runs_range=100:100:200)
   
    parameter_list = Dict()
    locs = 1
        
    i = 1
    for grids in grids_range
        for days in days_range
            for safety in safety_range
                for test_runs in test_runs_range
                    for rep in 1:reps
                        parameter_list[i] = (grids,days,safety,test_runs,rep)
                    end
                end
            end
        end
    end
    
    #@distributed
    for p in goodvals(parameter_list)
        grids,days,safety,test_runs,rep = p                
        ftr = test_runs
        rtr = ftr
        
        ratio = ratio_finder("none", 0, 0, 0, 0, 0.99, 1., 10.,
                             grids=grids, locs=locs, safety=safety, ftr=ftr, rtr=rtr, days=days)
        P2 = P2_finder("none", 0, 0, 0, 0, 0.99, 10, 99,
                  grids=grids, locs=locs, safety=safety, ftr=ftr, rtr=rtr, days=days)
        
        CSV.write(path*"/"*filename,Tables.table(tab),delim="\t",header=false)
    end   
end

# function that returns the minimum sustainant value achieved by a type of grid in the absence of adverse influences (prosumers)
# depending on the safety factor and number of test runs for the initial line capacities
function sustainant_test_2d(perf_days, grids, max_testruns, testruns_step, max_safety, safety_step) #todo
    
    delta = Dict()
    
    N = 100
    res = 60
    steps_per_day = Int(ceil(86400/res))
    solar_dict = smoother(solar_dict_default,res)
    consumer_dict = smoother(consumer_dict_default,res)
    solar_map = zeros(N)
    prosumption_ratio = 0
    
    caps = Dict()
    
    for g in 1:grids
        
        grid = create_grid(N, maximum([Int(floor(N*0.15)),1]), 0.15, 0.15, 0.2)
        slack_index = findmax(closeness_centrality(grid.graph))[2]

        test_series = balance_randomdraws(solar_map, solar_dict, consumer_dict, prosumption_ratio, max_testruns+1,
                                          slack_index)
        sus_series = balance_randomdraws(solar_map, solar_dict, consumer_dict, prosumption_ratio,
                                         perf_days*steps_per_day, slack_index)
        
        for safety in 1.:safety_step:max_safety
            for t in 1:testruns_step:max_testruns+1
                key = (t, safety)
                caps[key] = safety*max_any_flows(grid, test_series.balance[1:t,:])

                if haskey(delta,key) == false
                    delta[key] = []
                end
                sus = system_performance_series(grid, caps[key], sus_series, slack_index)
                append!(delta[key],sus.delta)
            end
        end
    end
    
    for key in keys(delta)
        delta[key] = minimum(delta[key])
    end
    
    return delta
end

# function to plot the minimum sustainant (computed by function above) as a heatmap
# values below 0.99 are set to 0
function dict2plot(dict)
    
    all_x = [key[1] for key in keys(dict)]
    all_y = [key[2] for key in keys(dict)]
    
    max_x = maximum(all_x)
    max_y = maximum(all_y)
    
    step_x = min_diff(all_x)
    step_y = round(min_diff(all_y),digits=2)
    
    for key in keys(dict)
        if dict[key] < 0.99
            dict[key] = 0
        end
    end
    
    xs = [t for t in 1:step_x:max_x]
    ys = [s for s in 1:step_y:max_y]
    p(x,y) = dict[(x,y)]
    
    plot(heatmap(xs, ys, p), title = "minimum delta = tau value",
         xlabel="captest runs", ylabel="safety factor")
    
end

# function that calculates the inverse of a quantile (also called quantile function = qf)
function qf(array,limit)
    frac = length(findall(x -> x<limit,array)) / length(array)
    return frac
end

# function that returns keys of a dictionary in a sorted order (rather than random)
function betterkeys(dict)
    return sort([key for key in keys(dict)])
end

# function to generate an estimated sample from a list of quantiles and their probabilities
# the histogram is created with the lowest (=worst) values possible
function worst_from_list(p_values, qs; counts=50*60*24*60)
    
    values = []
    for i in 1:(length(p_values)-1)
        
        low = qs[i]
        reps = counts * (p_values[i+1]-p_values[i])
        append!(values,fill(low,Int(ceil(reps))))
    end
    return values
end

# workaround function to add a super-title to an existing plot or grid of subplots
function superplot(args...; title="title", layout = length(args), size=(1200,800))
    
    y = [1,2,3] 
    diy_title = Plots.scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], Plots.text(title, :center)),
                              axis=false, leg=false, grid=false)

    Plots.plot( diy_title, Plots.plot(args..., layout=layout, size=size), layout=(2, 1))
    #layout=grid(2, 1, heights=[0.1, 0.9]) #maybe 1.5 fixes grid?
end

# function to generate a fake histogram (column plot) from a list of quantiles and their probabilities  
function fake_histogram_data(ps, qs)
    
    xs = []
    ys = []
    
    append!(xs,qs[1])
    append!(ys,ps[1]^2/ps[2])
    
    for i in 2:length(ps)
        
        append!(xs,[qs[i-1],qs[i]])
        append!(ys,[ps[i]-ps[i-1], ps[i]-ps[i-1]])
    end
    
    return (xs,ys)
end    

# function to create a contour plot from 3D data along with another subplot displaying the error
function contour_with_error(x, y, z, z_err, ;z_max=maximum(z), xlabel="x", ylabel="y", zlabel="z", elabel="z error",
                            cmap=:Spectral_11, err_cmap=:inferno, dpi=100,
                            xlog=true, frame_fac=0, sep_cmap=true, clevels=500, cexp=-10) 
    
    cmap = cgrad(cmap, scale=x->-log10(1-x+1*10. ^cexp))
    #err_cmap = cgrad(cmap, scale=x->-log10(1-x+1*10. ^(cexp/50)))
    
    clims = (0, z_max)
    err_clims = :none #(0,1)
    
    if sep_cmap == false
        err_cmap = cmap
        err_clims = clims
    end
    
    if xlog == true
        xscale = :log10
        xspan = log10(x[end]) - log10(x[1])
        xlims = (10^(x[1]-frame_fac*xspan), 10^(x[end]+frame_fac*xspan) )
    else
        xscale = :identity
        xspan = x[end]-x[1]
        xlims = (x[1]-frame_fac*xspan, x[end]+frame_fac*xspan)
    end
    
    yspan = y[end]-y[1]
    ylims = (y[1]-frame_fac*yspan, y[end]+yspan*frame_fac)
    
    vp = contour(x, y, z, xlabel=xlabel, ylabel=ylabel, color=cmap, fill=true, levels=clevels, colorbar=:right,
                 title=zlabel, yticks=y, xtickfontrotation=90, clims=clims, #xlims=xlims, ylims=ylims,
                 right_margin=dpi/100*5*mm, xscale=xscale, xticks=(x,x), framestyle=:box)
    
    ep = contour(x, y, z_err, xlabel=xlabel, ylabel=ylabel, color=err_cmap, fill=true, levels=clevels, colorbar=:right,
                 title=elabel, yticks=y, xtickfontrotation=90, clims=err_clims, #xlims=xlims, ylims=ylims,
                 right_margin=dpi/100*5*mm, xscale=xscale, xticks=(x,x), framestyle=:box)
    
    return plot(vp, ep, layout=@layout[a b], size=(1000,500))
end

# function that returns values from an array even for index values larger than the array length (by cycling back to the first array element)
function mod_array(number,array)
    return array[mod(number,1:length(array))]
end

# function to create two marginal line plots from 3D data
function marginal_plot(x, y, z, z_err; xlab="x", xlabel="x", ylabel="y", zlabel="z",cmap=:viridis, #:rainbow1 :RdYlBu_3 :rainbow_bgyr_35_85_c73_n256
                       xlog=true, zlog=false, frame_fac=0, dpi=100)#0.05)
    
    if xlog == true
        xscale = :log10
        xspan = log10(x[end]) - log10(x[1])
        xlims = ( 10^(x[1]-frame_fac*xspan), 10^(x[end]+frame_fac*xspan) )
    else
        xscale = :identity
        xspan = x[end]-x[1]
        xlims = (x[1]-frame_fac*xspan, x[end]+frame_fac*xspan)
    end
    
    if zlog == true
        zscale = :log10
        zspan = log10(maximum(z)) - log10(minimum(z))
        zlims = ( 10^(minimum(z)-frame_fac*zspan), 10^(maximum(z)+frame_fac*zspan) )
    else
        zscale = :identity
        zspan = maximum(z)-minimum(z)
        zlims = (minimum(z)-frame_fac*zspan, maximum(z)+zspan*frame_fac)
    end
    
    yspan = y[end]-y[1]
    ylims = (y[1]-frame_fac*yspan, y[end]+yspan*frame_fac)
    
    styles =  [:solid, :dash, :dot, :dashdot,]
    
    xp = plot(framestyle=:box)
    for i in 1:length(y)
        line_c = cgrad(cmap)[(i-1)/(length(y)-1)]
        line_s = mod_array(i,styles)
        plot!(x,z[i,:],xlabel=xlabel,ylabel=zlabel,label="$(ylabel) = $(y[i])",leg=:outertop, #xlims=xlims, ylims=zlims,
              linewidth=2,linecolor=line_c, xticks=(x,x), xtickfontrotation=90, xscale=xscale, yscale=zscale,
              ribbon=z_err, fillalpha=0.2, fillcolor=line_c, linestyle=line_s, 
              right_margin=dpi/100*5*mm, top_margin=dpi/100*5*mm)
    end
    
    yp = plot(framestyle=:box)
    for j in 1:length(x)
        line_c = cgrad(cmap)[(j-1)/(length(x)-1)]
        line_s = mod_array(j,styles)
        plot!(y,z[:,j],xlabel=ylabel,ylabel=zlabel,label="$(xlab) = $(x[j])",leg=:outertop, #xlims=ylims, ylims=zlims,
              linewidth=2,linecolor=line_c, xticks=y, xtickfontrotation=90, xscale=:identity,
              ribbon=z_err, fillalpha=0.2, fillcolor=line_c, linestyle=line_s,
              right_margin=dpi/100*15*mm, top_margin=dpi/100*5*mm)
    end
    
    return plot(xp, yp, layout=@layout[a b], size=(1000,600))
    
end

# function to get all element of a list except for one
function find_missing_values(selectable, selected)
    
    vals = []
    for i in selectable
        if (i in selected) == false
            append!(vals,i)
        end
    end
    return vals
end

# function to delete a row of a data frame
function delete_rows(df, name, value)
    df[df[!,name].!=value,:]
end

# function that returns the frequency of the sustainant deviating from the ideal value by more than a threshold
function resilience_freq(mean_weighted_deviations; threshold=1/99/365) #10^-5)
    return length(filter(x -> x<=threshold, mean_weighted_deviations))/length(mean_weighted_deviations)
end

# multiple possible influence density functions

normy = 1#/(9900-4950)/(50-0.995)
decreasing_chi(np, rp) = normy * (100-np) * (10-rp) 

normy2 = 1# /(10-0.1)/99
uniform_chi(np ,rp) = normy2

indicator_chi(np, rp) = 1

# function that calculates the resilience value of a given response strategy from the QMC samples covering the influence space
# the QMC samples are read from their folder and a new .csv file containing the summary is created in the parent directory
function qmc_summary(results_path; uniform_chi=uniform_chi, non_uniform_chi=decreasing_chi)
    
    @manipulate for b = button("Summarize QMC sampling results"; value=0)
        if b > 0
            b = 0

            results_directory = readdir(results_path)
            for allocation_folder in results_directory

                if occursin("csv",allocation_folder) == false

                    allocation_directory = readdir(results_path*"/"*allocation_folder)
                    @showprogress for parameter_folder in allocation_directory

                        if occursin("csv",parameter_folder) == false

                            samples = Dict()        
                            counts = 0
                            non_uni_prob_sum = 0
                            uni_delta_sum, uni_tau_sum, uni_delta_squaresum, uni_tau_squaresum = 0, 0, 0, 0
                            non_uni_delta_sum, non_uni_tau_sum, non_uni_delta_squaresum, non_uni_tau_squaresum = 0, 0, 0, 0

                            folder = readdir(results_path*"/"*allocation_folder*"/"*parameter_folder)
                            for file in folder
                                
                                thing = results_path*"/"*allocation_folder*"/"*parameter_folder*"/"*file
                                #params = CSV.read(thing, DataFrame, transpose=true,
                                                   #limit=1, drop=collect(24:45))
                                devs = CSV.read(thing, DataFrame, skipto=27, header=26)

                                delta_devs = devs."delta mean weighted deviation"
                                tau_devs = devs."tau mean weighted deviation"

                                delta_freq = 1.0 #resilience_freq(delta_devs)
                                tau_freq = resilience_freq(tau_devs)

                                file = string(file)
                                P2 = value_from_filename(file, "P2")
                                ratio = value_from_filename(file, "rp")
                                key = (P2, ratio)
                                
                                counts += 1
                                
                                uni_inf_prob = uniform_chi(P2,ratio)
                                non_uni_inf_prob = non_uniform_chi(P2,ratio)
                                non_uni_prob_sum += non_uni_inf_prob
                                
                                samples[key] = [uni_inf_prob, non_uni_inf_prob, delta_freq, tau_freq]
                                
                                uni_delta_sum += delta_freq * uni_inf_prob
                                uni_tau_sum += tau_freq * uni_inf_prob
                                
                                uni_delta_squaresum += (delta_freq * uni_inf_prob)^2
                                uni_tau_squaresum += (tau_freq * uni_inf_prob)^2
                                
                                non_uni_delta_sum += delta_freq * non_uni_inf_prob
                                non_uni_tau_sum += tau_freq * non_uni_inf_prob
                                
                                non_uni_delta_squaresum += (delta_freq * non_uni_inf_prob)^2
                                non_uni_tau_squaresum += (tau_freq * non_uni_inf_prob)^2
                                
                            end
                            
                            max_P2 = value_from_filename(parameter_folder, "max_P2")
                            max_ratio = value_from_filename(parameter_folder, "max_rp")
                            allocation = string_value_from_filename(parameter_folder, "alloc")
                            line_budget_factor = value_from_filename(parameter_folder, "line_bf")
                            new_lines = value_from_filename(parameter_folder, "new_l")
                            battery_budget_factor = value_from_filename(parameter_folder, "batt_bf")
                            battery_reliance = value_from_filename(parameter_folder, "b_rel")

                            qmc_area = max_P2 * (max_ratio-0.1)
                            
                            
                            #normalizing a posteriori
                            non_uni_delta_sum *= (counts/qmc_area/non_uni_prob_sum)
                            non_uni_tau_sum *= (counts/qmc_area/non_uni_prob_sum)
                            non_uni_delta_squaresum *= (counts/qmc_area/non_uni_prob_sum)^2
                            non_uni_tau_squaresum *= (counts/qmc_area/non_uni_prob_sum)^2
                            
                            
                            
                            uni_delta_integral = uni_delta_sum * qmc_area/counts
                            uni_tau_integral = uni_tau_sum * qmc_area/counts
                            
                            non_uni_delta_integral = non_uni_delta_sum * qmc_area/counts
                            non_uni_tau_integral = non_uni_tau_sum * qmc_area/counts
                            
                            if (uni_delta_squaresum - (1/counts)*uni_delta_sum^2) < 0
                                uni_delta_error = 0
                            else
                                uni_delta_error = qmc_area /sqrt(counts^2-counts) * 
                                                  sqrt( uni_delta_squaresum - (1/counts)*uni_delta_sum^2 )
                            end
                            if (uni_tau_squaresum - (1/counts)*uni_tau_sum^2) < 0
                                uni_tau_error = 0
                            else
                                uni_tau_error = qmc_area /sqrt(counts^2-counts) * 
                                                sqrt( uni_tau_squaresum - (1/counts)*uni_tau_sum^2 )
                            end
                            
                            if (non_uni_delta_squaresum - (1/counts)*non_uni_delta_sum^2) < 0
                                non_uni_delta_error = 0
                            else
                                non_uni_delta_error = qmc_area /sqrt(counts^2-counts) * 
                                                      sqrt( non_uni_delta_squaresum - (1/counts)*non_uni_delta_sum^2 )
                            end
                            if (non_uni_tau_squaresum - (1/counts)*non_uni_tau_sum^2) < 0
                                non_uni_tau_error = 0
                            else
                                non_uni_tau_error = qmc_area /sqrt(counts^2-counts) * 
                                                    sqrt( non_uni_tau_squaresum - (1/counts)*non_uni_tau_sum^2 )
                            end
                            
                            output_file = ""

                            if allocation == "none"
                                output_file = "res_for_alloc=$(allocation)"*
                                                  ",_max_rp=$(max_ratio),_max_P2=$(max_P2)"
                            elseif line_budget_factor == 0
                                output_file = "res_for_alloc=$(allocation)"*
                                                  ",_batt_bf=$(battery_budget_factor),_b_rel=$(battery_reliance)"*
                                                  ",_max_rp=$(max_ratio),_max_P2=$(max_P2)"
                            elseif battery_budget_factor == 0
                                output_file = "res_for_alloc=$(allocation)"*
                                                  ",_line_bf=$(line_budget_factor),_new_l=$(new_lines)"*
                                                  ",_max_rp=$(max_ratio),_max_P2=$(max_P2)"
                            end

                            output_file = replace(output_file,Pair(".","p"))*".csv"

                            parameter_names = ["qmc_area",
                                               "uniform_delta_integral", "uniform_delta_error",
                                               "uniform_tau_integral", "uniform_tau_error",
                                               "non_uniform_delta_integral", "non_uniform_delta_error",
                                               "non_uniform_tau_integral", "non_uniform_tau_error"]
                            parameter_values = [qmc_area,
                                                uni_delta_integral, uni_delta_error,
                                                uni_tau_integral, uni_tau_error,
                                                non_uni_delta_integral, non_uni_delta_error,
                                                non_uni_tau_integral, non_uni_tau_error]
                            empty_column = [" ", " ", " ", " ", " ", " ", " ", " ", " "]
                            list = hcat(parameter_names,parameter_values,
                                        empty_column, empty_column, empty_column, empty_column)
                            emptyrow = [" " " " " " " " " " " "]
                            sample_list = ["P2" "prosumer_ratio" "uniform_influence_probability" "non_uniform_influence_probability" "delta_resilience_assessment" "tau_resilience_assessment"]

                            for key in keys(samples)
                                newline = zeros(1,6)
                                newline[1] = key[1]
                                newline[2] = key[2]
                                newline[3] = samples[key][1]
                                newline[4] = samples[key][2]
                                newline[5] = samples[key][3]
                                newline[6] = samples[key][4]
                                sample_list = vcat(sample_list,newline)
                            end

                            tab = vcat(list,emptyrow,sample_list)
                            CSV.write(results_path*"/"*allocation_folder*"/"*output_file, Tables.table(tab),
                                      delim="\t", header=false)
                        end
                    end
                end
            end
        end
    end
end

# struct to contain mean and error of a result
struct error_struct
    mean
    error
end

# function that calculates the resilience curve of a given response strategy from its resilience value at different parameter instances
# the resilience values are the read from the files created by the above function and a new .csv file containing the curve summary is created in the parent directory
function strategy_summary(results_path)
    
    @manipulate for b = button("Summarize allocation results"; value=0)
        if b > 0
            b = 0

            results_directory = readdir(results_path)
            filter!(e->!occursin("csv",e), results_directory)
            @showprogress for allocation_folder in results_directory

                uni_delta_res, uni_tau_res = Dict(), Dict()
                non_uni_delta_res, non_uni_tau_res = Dict(), Dict()
                variables = []

                allocation_directory = readdir(results_path*"/"*allocation_folder)
                filter!(e->occursin("csv",e), allocation_directory)
                for file in allocation_directory

                    key = []

                    allocation = string_value_from_filename(file, "alloc")
                    line_budget_factor = value_from_filename(file, "line_bf")
                    new_lines = value_from_filename(file, "new_l")
                    battery_budget_factor = value_from_filename(file, "batt_bf")
                    battery_reliance = value_from_filename(file, "b_rel")

                    if allocation == "none"
                        key = []
                        variables = []
                    elseif line_budget_factor == 0
                        key = [battery_budget_factor, battery_reliance]
                        variables = ["battery_budget_factor", "battery_reliance"]
                    elseif battery_budget_factor == 0
                        key = [line_budget_factor, new_lines]
                        variables = ["line_budget_factor", "new_lines"]
                    #else
                        #key = (battery_budget_factor, battery_reliance, line_budget_factor, new_lines)
                        #variables = ["battery_budget_factor", "battery_reliance", "line_budget_factor", "new_lines"]
                    end
                    
                    thing = results_path*"/"*allocation_folder*"/"*file
                    qmc = CSV.read(thing, DataFrame, transpose=true, limit=1, drop=collect(11:268))

                    uni_delta_res[key] = error_struct(qmc."uniform_delta_integral"[1],qmc."uniform_delta_error"[1])
                    uni_tau_res[key] = error_struct(qmc."uniform_tau_integral"[1],qmc."uniform_tau_error"[1])
                    non_uni_delta_res[key] = error_struct(qmc."non_uniform_delta_integral"[1],qmc."non_uniform_delta_error"[1])
                    non_uni_tau_res[key] = error_struct(qmc."non_uniform_tau_integral"[1],qmc."non_uniform_tau_error"[1])
                end

                output_file = "summary_for_"*allocation_folder*".csv"

                lines = length(uni_delta_res)+1
                columns = length(variables)+8

                tab = Array{Any}(undef, lines, columns)
                tab[1,1:length(variables)] = variables
                tab[1,end-7:end] = ["uniform_delta_resilience","uniform_delta_error",
                                    "uniform_tau_resilience","uniform_tau_error",
                                    "non_uniform_delta_resilience","non_uniform_delta_error",
                                    "non_uniform_tau_resilience","non_uniform_tau_error"]

                i = 2
                for key in betterkeys(uni_delta_res)
                    tab[i,1:length(variables)] = key
                    tab[i,end-7] = uni_delta_res[key].mean
                    tab[i,end-6] = uni_delta_res[key].error
                    tab[i,end-5] = uni_tau_res[key].mean
                    tab[i,end-4] = uni_tau_res[key].error
                    tab[i,end-3] = non_uni_delta_res[key].mean
                    tab[i,end-2] = non_uni_delta_res[key].error
                    tab[i,end-1] = non_uni_tau_res[key].mean
                    tab[i,end] = non_uni_tau_res[key].error
                    i += 1
                end

                CSV.write(results_path*"/"*output_file,Tables.table(tab),delim="\t",header=false)
            end
        end
    end
end

# function to generate the resilience value histogram of the status quo (non-adapted) power grid and save it as a .csv file
function status_quo_summary(path)
    
    @manipulate for b = button("Summarize status quo results"; value=0)
        if b > 0
            b = 0
    
            coordinate_list = readdir(path)

            @showprogress for coordinate in coordinate_l^ist

                file_list = readdir(path*"/"*coordinate)
                filter!(e->!occursin("frac",e), file_list)

                worst_dist = []

                for file in file_list
                    
                    thing = path*"/"*coordinate*"/"*file
                    params = CSV.read(thing, DataFrame, transpose=true,
                                      limit=1, drop=collect(24:45))
                    quants = CSV.read(thing, DataFrame, skipto=26, header=25)

                    s = params."line_safety"[1]
                    ftr = params."flowtest_runs"[1]

                    ps = quants."p-value"
                    qs = quants."delta p-quantile"

                    dist = worst_from_list(ps, qs)
                    append!(worst_dist, dist)
                end

                frac = qf(worst_dist, 1)
                CSV.write(path*"/"*coordinate*"/frac.csv",Tables.table([frac frac; frac frac]),delim="\t",header=false)       
            end
        end
    end
end

# function to select a directory interactively
function choose_output_path(array)
    @manipulate for path in array
        global output_path = path
        println(path)
    end
end

# function to select a baseline value interactively
# the baseline could be e.g. the resilience of the non-adapted power grid
function choose_baseline(array)
    @manipulate for tuple in array
        global baseline = tuple
        println(tuple)
    end
end

# function to interactively show and export scatter plots of the QMC samples estimating the resilience basins
function qmc_visualizer(path; sizzle=(1200,800),fontfac=1, dpi=100, cmap=:inferno) ###change cmap?
    
    IJulia.clear_output()
    
    directory = readdir(path)
    filter!(e->!occursin("csv",e), directory)
    
    @manipulate for allocation = dropdown( directory, label = "allocation"),
                    #influence_density in ["uniform", "non_uniform"],
                    sustainant in ["tau", "delta"]
        
        folder = readdir(path*"/"*allocation)
        filter!(e->occursin("csv",e), folder)
        
        LBF, NL, BBF, BR = [], [], [], []
        
        for filename in folder
            if occursin("res_for_",filename) == true
                append!(LBF,value_from_filename(filename, "line_bf"))
                append!(NL,value_from_filename(filename, "new_l"))
                append!(BBF,value_from_filename(filename, "batt_bf"))
                append!(BR,value_from_filename(filename, "b_rel"))
            end
        end
        
        LBF = sort(unique(LBF))
        NL = sort(unique(NL))
        BBF = sort(unique(BBF))
        BR = sort(unique(BR))
        
        @manipulate for line_budget_factor in LBF, new_lines in NL,
                        battery_budget_factor in BBF, battery_reliance in BR,
                        cbtoggle = toggle(label="plot colorbar", value=0),
                        b = button("Save as .PNG"; value=0),
                        size_h = spinbox(label="horizontal plot size"; value=1220),
                        size_v = spinbox(label="vertical plot size"; value=1000)
            
            
            plot_filename = ""
            for filename in folder
                if (occursin("res_for_",filename) == true &&
                    value_from_filename(filename, "line_bf") == line_budget_factor &&
                    value_from_filename(filename, "new_l") == new_lines &&
                    value_from_filename(filename, "batt_bf") == battery_budget_factor &&
                    value_from_filename(filename, "b_rel") == battery_reliance )
                    
                    plot_filename = filename
                end
            end
            
            thing = path*"/"*allocation*"/"*plot_filename
            #values = CSV.read(thing, DataFrame, transpose=true, limit=1, drop=collect(10:267))
            samples = CSV.read(thing, DataFrame, skipto=12, header=11)
            
            #res_col = "$(influence_density)_$(sustainant)_integral"
            #err_col = "$(influence_density)_$(sustainant)_error"
            
            #res = round( values[:,res_col][1], digits=2)
            #err = ceil( values[:,err_col][1], digits=2)

            plot_filename = plot_filename[1:end-4]

            max_P2 = value_from_filename(plot_filename, "max_P2")
            max_ratio = value_from_filename(plot_filename, "max_rp")
            
            #x_delta_stable, y_delta_stable = [], []
            #x_delta_unstable, y_delta_unstable = [], []
            #x_tau_stable, y_tau_stable = [], []
            #x_tau_unstable, y_tau_unstable = [], []
            rps, nps, freqs = [], [], []
            
            number = size(samples)[1]
            for L in 1:number
                row = samples[L,:]
                append!(rps,row."prosumer_ratio"[1])
                append!(nps,row."P2"[1])
                if sustainant == "delta"
                    append!(freqs,row."delta_resilience_assessment"[1])
                elseif sustainant == "tau"
                    append!(freqs,row."tau_resilience_assessment"[1])
                end
            end

            col = [(cgrad(cmap)[(freqs[i]/1)]) for i in 1:length(freqs)]
            
            col_func(rp,np) = freqs[ intersect( findall(e->e==rp, rps), findall(e->e==np, nps) )[1] ]
            
            cbool = :none
            xticks=:auto
            if cbtoggle == 1
                cbool = :right
                #xticks = [1]
            end
            
            p = plot(size=sizzle, dpi=dpi,
                     xlabel="\$r_p\$",
                     #ylim=(0,ceil(max_P2, digits=-1)), xlim=(0,ceil(max_ratio,digits=1)),
                     ylabel="\$n_p\$", legend=false, ytickfontrotation=90,
                     labelfontsize=11*fontfac, legendfontsize=8*fontfac, tickfontsize=8*fontfac, framestyle=:box)
            scatter!(rps, nps, markersize=25, colorbar=cbool, #color=col,
                     marker_z = col_func, cmap =:inferno, colorbar_title="\$\\alpha\$",
                     colorbar_tickfontsize=8*fontfac, colorbar_titlefontsize=11*fontfac, clims=(0,1), xticks=xticks,
                     right_margin=5*mm, markerstrokecolor=:grey)
                
                
            if b > 0
                savefig("Plots/$(sustainant)-based_QMC_samples"*plot_filename[4:end])
                b = 0
            end
            
            p
        end
    end
end

# function to convert filename strings into readable text (without underscores) and format mathematical symbols in LaTeX
function texify(string)
    if string == "line_budget_factor"
        return "\$\\phi_L\$"
    elseif string == "new_lines"
        return "\$\\varepsilon\$" #"number of new lines \$\\lambda\$"
        
    elseif string == "battery_budget_factor"
        return "\$\\phi_B\$"
    elseif string == "battery_reliance"
        return "\$\\lambda\$"
        
    elseif string == "line_adaptation_uniform"
        return "uniform line adaptation"
    elseif string == "battery_adaptation_uniform"
        return "uniform battery adaptation"
        
    elseif string == "line_adaptation_closeness"
        return "closeness-based line adaptation"
    elseif string == "battery_adaptation_closeness"
        return "closeness-based battery adaptation"
        
    elseif string == "line_adaptation_risktest"
        return "flowtest-based line adaptation"
    elseif string == "battery_adaptation_risktest"
        return "flowtest-based battery adaptation"
    
    end
end

# function to interactively show and export line plots of the resilience curves showing the parameter dependency of the response strategies
function strategy_visualizer(path)
    
    directory = readdir(path)
    filter!(e->occursin("csv",e), directory)
    
    better_directory = copy(directory)
    for i in 1:length(directory)
        better_directory[i] = directory[i][13:end-4]
    end
    
    @manipulate for allocation = dropdown( better_directory, label = "allocation"),
                    diagram in ["colormap", "budget plot", "secondary plot"],
                    density in ["probability (non-uniform)", "possibility (uniform)"]
        
        thing = path*"/summary_for_"*allocation*".csv"
        table = CSV.read(thing, DataFrame, header=true)
        
        if allocation == "no_adaptation"
            println("uniform tau resilience: $(table."uniform_tau_resilience"[1]) +- $(table."uniform_tau_error"[1])")
            println()
            println("non-uniform tau resilience: $(table."non_uniform_tau_resilience"[1]) +- $(table."non_uniform_tau_error"[1])")
        else
            xlab = names(table)[1]
            ylab = names(table)[2]
            
            x = unique(table[:,xlab])
            y = unique(table[:,ylab])
            
            @manipulate for xchoice = toggles(x, value=x, label="filter values for $(xlab)"),
                            ychoice = toggles(y, value=y, label="filter values for $(ylab)")

                sort!(xchoice)
                sort!(ychoice)

                LX = length(xchoice)
                LY = length(ychoice)

                x_hidden = find_missing_values(x, xchoice)
                y_hidden = find_missing_values(y, ychoice)

                newtable = copy(table)

                for i in x_hidden
                    newtable = delete_rows(newtable, xlab, i)
                end
                for j in y_hidden
                    newtable = delete_rows(newtable, ylab, j)
                end

                if density == "possibility (uniform)"
                    delta = reshape(newtable."uniform_delta_resilience", (LY,LX))
                    delta_error = reshape(newtable."uniform_delta_error", (LY,LX))
                    tau = reshape(newtable."uniform_tau_resilience", (LY,LX))
                    tau_error = reshape(newtable."uniform_delta_error", (LY,LX))
                    R_max = 980.1
                elseif density == "probability (non-uniform)"
                    delta = reshape(newtable."non_uniform_delta_resilience", (LY,LX))
                    delta_error = reshape(newtable."non_uniform_delta_error", (LY,LX))
                    tau = reshape(newtable."non_uniform_tau_resilience", (LY,LX))
                    tau_error = reshape(newtable."non_uniform_delta_error", (LY,LX))
                    R_max = 1
                end                
                
                R = []
                R_err = []
                R_label =""

                R = tau
                R_err = tau_error
                R_label = "\$R\$"          

                xlabel = texify(xlab)
                ylabel = texify(ylab)
                
                if diagram == "colormap"
                    @manipulate for size_h = spinbox(label="horizontal plot size"; value=1220),
                                    size_v = spinbox(label="vertical plot size"; value=1000),
                                    budget_scale in ["log", "linear"],
                                    secondary_scale in ["linear", "log"],
                                    font_fac = spinbox(label="fontsize factor"; value=2.),
                                    clevels = spinbox(label="color levels"; value=20)
                        
                        #R normalization according to max in x
                        for i in 1:size(R)[2]
                            R[:,i] /= maximum(R[:,i])
                        end
                        
                        if budget_scale == "log"
                            xscale = :log10
                            xlabel = texify(xlab)*" (logarithmic scale)"
                        else
                            xscale = :identity
                            xlabel = texify(xlab)
                        end
                        
                        if secondary_scale == "log"
                            yscale = :log10
                            ylabel = texify(ylab)*" (logarithmic scale)"
                        else
                            yscale = :identity
                            ylabel = texify(ylab)
                        end
                        
                        contour(x, y, R, xlabel=xlabel, ylabel=ylabel, color=:inferno, fill=true,
                                colorbar=:right, xticks=(x,x), yticks=(y,y), xtickfontrotation=90,
                                xscale=xscale, yscale=yscale, colorbar_title="\$R/R_{opt(\\phi_L)}\$",
                                #clims=(0, R_max),
                                framestyle=:box, size=(size_h, size_v),
                                labelfontsize=11*font_fac, tickfontsize=8*font_fac,
                                colorbar_tickfontsize=8*font_fac, colorbar_titlefontsize=11*font_fac,
                                levels=clevels)

                    end
                    
                elseif diagram == "budget plot"
                    @manipulate for size_h = spinbox(label="horizontal plot size"; value=1000),
                                    size_v = spinbox(label="vertical plot size"; value=1000),
                                    budget_scale in ["log", "linear"],
                                    R_scale in ["log", "linear"],
                                    cmap in [:viridis, :thermometer, :coolwarm],
                                    legend_loc in [:topleft, :bottomright],
                                    font_fac = spinbox(label="fontsize factor"; value=2.),
                                    alph = spinbox(label="ribbon alpha"; value=0.1)
                        
                        if budget_scale == "log"
                            xscale = :log10
                            xlabel = texify(xlab)*" (logarithmic scale)"
                        else
                            xscale = :identity
                            xlabel = texify(xlab)
                        end
                        
                        if R_scale == "log"
                            Rscale = :log10
                            R_label = "\$R\$ (logarithmic scale)"
                        else
                            Rscale = :identity
                            R_label = "\$R\$"
                        end
                        
                        styles =  [:solid, :dash, :dot, :dashdot,]
    
                        xp = plot(framestyle=:box)
                        for i in 1:length(y)
                            line_c = cgrad(cmap)[(i-1)/(length(y)-1)]
                            line_s = mod_array(i,styles)
                            plot!(x, R[i,:], xlabel=xlabel, ylabel=R_label, label="$(ylabel) = $(y[i])",
                                  leg=legend_loc, linewidth=2, linecolor=line_c, xticks=(x,x),
                                  xtickfontrotation=90, xscale=xscale, yscale=Rscale, ribbon=R_err,
                                  fillalpha=alph, fillcolor=line_c, linestyle=line_s, size=(size_h, size_v),
                                  labelfontsize=11*font_fac, legendfontsize=7*font_fac, tickfontsize=8*font_fac)
                        end
                        for i in 1:length(y)
                            line_c = cgrad(cmap)[(i-1)/(length(y)-1)]
                            line_s = mod_array(i,styles)
                            plot!(x, R[i,:], label=false, linewidth=2, linecolor=line_c, linestyle=line_s)
                        end
                        xp
                    end                   
                    
                elseif diagram == "secondary plot"
                    @manipulate for size_h = spinbox(label="horizontal plot size"; value=1000),
                                    size_v = spinbox(label="vertical plot size"; value=1000),
                                    secondary_scale in ["linear", "log"],
                                    R_scale in ["linear", "log"],
                                    cmap in [:redsblues, :vikO, :thermometer, :viridis, :coolwarm],
                                    legend_loc in [:topleft, :bottomright],
                                    font_fac = spinbox(label="fontsize factor"; value=2.),
                                    alph = spinbox(label="ribbon alpha"; value=0.1)
                                                    
                    if secondary_scale == "log"
                            yscale = :log10
                            ylabel = texify(ylab)*" (logarithmic scale)"
                        else
                            yscale = :identity
                            ylabel = texify(ylab)
                        end
                        
                        if R_scale == "log"
                            Rscale = :log10
                            R_label = "\$R\$ (logarithmic scale)"
                        else
                            Rscale = :identity
                            R_label = "\$R\$"
                        end
                        
                        styles =  [:solid, :dash, :dot, :dashdot,]
                    
                        yp = plot(framestyle=:box)
                        for j in 1:length(x)
                            line_c = cgrad(cmap)[(j-1)/(length(x)-1)]
                            line_s = mod_array(j,styles)
                            plot!(y, R[:,j],xlabel=ylabel,ylabel=R_label, label="$(xlabel) = $(x[j])", 
                                  leg=legend_loc, linewidth=2,linecolor=line_c, xticks=(y,y), xtickfontrotation=90,
                                  xscale=yscale, yscale=Rscale, ribbon=R_err, fillalpha=alph,
                                  fillcolor=line_c, linestyle=line_s, size=(size_h, size_v),
                                  labelfontsize=11*font_fac, legendfontsize=7*font_fac, tickfontsize=8*font_fac)
                        end
                        for j in 1:length(x)
                            line_c = cgrad(cmap)[(j-1)/(length(x)-1)]
                            line_s = mod_array(j,styles)
                            plot!(y, R[:,j], label=false, linewidth=2, linecolor=line_c, linestyle=line_s)
                        end
                        yp
                        
                    end     
                end    
            end
        end
    end
end

# deprecated version of the above plotting function
function strategy_visualizer_old(path; subtract=[0,0], exp=-1, sizzle=(1000,1333), dpi=100, fontfac=1)
    
    IJulia.clear_output()
    
    directory = readdir(path)
    filter!(e->occursin("csv",e), directory)
    
    better_directory = copy(directory)
    for i in 1:length(directory)
        better_directory[i] = directory[i][13:end-4]
    end
    
    @manipulate for allocation = dropdown( better_directory, label = "allocation")
        
        thing = path*"/summary_for_"*allocation*".csv"
        table = CSV.read(thing, DataFrame, header=true)
        
        if allocation == "no_adaptation"
            println("uniform delta resilience: $(table."uniform_delta_resilience"[1]) +- $(table."uniform_delta_error"[1])")
            println("uniform tau resilience: $(table."uniform_tau_resilience"[1]) +- $(table."uniform_tau_error"[1])")
            println()
            println("non-uniform delta resilience: $(table."non_uniform_delta_resilience"[1]) +- $(table."non_uniform_delta_error"[1])")
            println("non-uniform tau resilience: $(table."non_uniform_tau_resilience"[1]) +- $(table."non_uniform_tau_error"[1])")
        else
            xlab = names(table)[1]
            ylab = names(table)[2]
            
            x = unique(table[:,xlab])
            y = unique(table[:,ylab])

            @manipulate for sustainant in ["tau", "delta"],
                            mode in ["probability (non-uniform)", "possibility (uniform)"],
                            error_mode in ["absolute", "relative"],
                            resilience_colorscale_exp = dropdown(-20:1: -1, label="resilience_colorscale_exp",
                                                                 value=exp),
                            resilience_colorlevels = spinbox(100:100:1200, label="resilience_colorlevels",value=500),
                            budget_scale in ["logarithmic", "linear"],
                            R_scale in ["logarithmic", "linear"],
                            xchoice = toggles(x, value=x, label="filter values for $(xlab)"),
                            ychoice = toggles(y, value=y, label="filter values for $(ylab)"),
                            margs = toggle(label="plot marginals", value=1),
                            b = button("Save as .PNG"; value=0)

                sort!(xchoice)
                sort!(ychoice)

                LX = length(xchoice)
                LY = length(ychoice)

                x_hidden = find_missing_values(x, xchoice)
                y_hidden = find_missing_values(y, ychoice)

                newtable = copy(table)

                for i in x_hidden
                    newtable = delete_rows(newtable, xlab, i)
                end
                for j in y_hidden
                    newtable = delete_rows(newtable, ylab, j)
                end

                if mode == "possibility (uniform)"
                    delta = reshape(newtable."uniform_delta_resilience", (LY,LX))
                    delta_error = reshape(newtable."uniform_delta_error", (LY,LX))
                    tau = reshape(newtable."uniform_tau_resilience", (LY,LX))
                    tau_error = reshape(newtable."uniform_delta_error", (LY,LX))
                    z_max = 980.1
                elseif mode == "probability (non-uniform)"
                    delta = reshape(newtable."non_uniform_delta_resilience", (LY,LX))
                    delta_error = reshape(newtable."non_uniform_delta_error", (LY,LX))
                    tau = reshape(newtable."non_uniform_tau_resilience", (LY,LX))
                    tau_error = reshape(newtable."non_uniform_delta_error", (LY,LX))
                    z_max = 1
                end

                z = []
                title = ""
                plot_name = ""

                if sustainant == "delta"
                    z = delta .-subtract[1]
                    z_err = delta_error
                    sus = "\$\\delta\$"
                    R = "\$R_\\delta\$"
                elseif sustainant == "tau"
                    z = tau .-subtract[2]
                    z_err = tau_error
                    sus = "\$\\tau\$"
                    R = "\$R_\\tau\$"
                end                
                
                marg_err = copy(z_err)
                if error_mode == "relative"
                    z_err = z_err./z
                    sep_cmap = true
                    elabel = "relative error of "*R
                elseif error_mode == "absolute"
                    sep_cmap = true
                    elabel = "absolute error \$\\sigma'_{R_\\$sustainant}\$"
                end

                zlabel = "\$R\$"
                plot_name = "$(sustainant)-based_resilience_for_$(texify(allocation))"

                xlog = false
                zlog=false
                xlabel = texify(xlab)
                leg_xlabel = xlabel
                ylabel = texify(ylab)
                if budget_scale == "logarithmic"
                    xlog = true
                    xlabel = xlabel*" (log scale)"
                end
                if R_scale == "logarithmic"
                    zlog = true
                    zlabel = zlabel*" (log scale)"
                end

                if minimum(size(z)) < 2
                    println("not enough parameter combos. Minimum is 2x2")
                else                    
                    cp = contour_with_error(xchoice, ychoice, z, z_err; z_max=z_max, xlabel=xlabel, ylabel=ylabel,
                                            zlabel=zlabel, elabel=elabel, xlog=xlog, sep_cmap=sep_cmap, dpi=dpi,
                                            clevels=resilience_colorlevels, cexp=resilience_colorscale_exp)
                    if margs == 0

                        p = cp

                    else

                        mp = marginal_plot(xchoice, ychoice, z, marg_err; dpi=dpi,
                                           xlab=leg_xlabel, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                                           xlog=xlog, zlog=zlog)

                        p = plot(cp, mp, layout=@layout[a; b], size=sizzle, dpi=dpi,
                                 labelfontsize=11*fontfac, legendfontsize=8*fontfac, tickfontsize=8*fontfac,
                                 titlefontsite=14*fontfac)

                    end

                    if b > 0
                        savefig("Plots/"*plot_name)
                        b = 0
                    end

                    p
                end
            end
        end
    end
end

# function that performs a glm fit on the resilience curve of a given response strategy
function glm_fit(path; subtract=[0,0])
    
    directory = readdir(path)
    filter!(e->occursin("csv",e), directory)
    
    better_directory = copy(directory)
    for i in 1:length(directory)
        better_directory[i] = directory[i][13:end-4]
    end
    
    @manipulate for allocation = dropdown( better_directory, label = "allocation"),
                    sustainant in ["tau", "delta"]
                    #mode
                    #errors?
        
        thing = path*"/summary_for_"*allocation*".csv"
        table = CSV.read(thing, DataFrame, header=true)
        
        table."delta_resilience" .-=subtract[1]
        table."tau_resilience" .-=subtract[2]
            
        if sustainant == "delta"
            newtable = rename(table, (names(table)[3] => "res" ))
        elseif sustainant == "tau"
            newtable = rename(table, (names(table)[4] => "res" ))      
        end
        
        newtable = rename(newtable, (names(newtable)[1] => "x" ), ( names(newtable)[2] => "y" ))
        
        @manipulate for link_function = dropdown( [IdentityLink(), CauchitLink(), CloglogLink(), InverseLink(),
                                          InverseSquareLink(), LogitLink(), LogLink(),
                                          ProbitLink(), SqrtLink()] , label="link function"),
                        family = dropdown( [Normal(), Bernoulli(), Binomial(), Gamma(), Poisson()] , label="family"),
                        b = button("import fit formula"; value=0)
            if b > 0
                b = 0
                
                fit = glm(form, newtable, family, link_function)
                
                println(fit)
                #println("R-squared: ", r2(fit, variant=:devianceratio))
            end
        end        
    end
end

# function that performs a nonlinear least squares fit on the resilience curve of a given response strategy
function lsq_fit(path; subtract=[0,0])
    
    directory = readdir(path)
    filter!(e->occursin("csv",e), directory)
    
    better_directory = copy(directory)
    for i in 1:length(directory)
        better_directory[i] = directory[i][13:end-4]
    end
    
    @manipulate for allocation = dropdown( better_directory, label = "allocation"),
                    sustainant in ["tau", "delta"]
                    #mode
                    #errors?
        
        thing = path*"/summary_for_"*allocation*".csv"
        table = CSV.read(thing, DataFrame, header=true)
        
        table."delta_capacity" .-=subtract[1]
        table."tau_capacity" .-=subtract[2]
        
        xlabel = names(table)[1]
        ylabel = names(table)[2]
        
        xy = hcat(table[:,xlabel], table[:,ylabel])
        z = []
        
        if sustainant == "delta"
            z = table."delta_capacity"
        elseif sustainant == "tau"
            z = table."tau_capacity"
        end
        
        @manipulate for cb = toggle(label="use parameter bounds", value=1),
                        b = button("apply fit formula"; value=0)
            if b > 0
                b = 0
                
                if cb == 0
                    fit = curve_fit(model, xy, z, Float64.(p0))
                else
                    fit = curve_fit(model, xy, z, Float64.(p0),
                                    lower=Float64.(lower_bounds), upper = Float64.(upper_bounds))
                end

                params = DataFrame(parameter_value = fit.param, standard_deviation = stderror(fit))
                println(sustainant*" resilience for "*allocation*":")
                println(params)
                println()
            end            
        end        
    end
end

# function to compare the two differently defined resilience measures (delta and tau) in a scatter plot
# the response strategies and their parameters are encoded as shape, color, and size
function resilience_scatter(path; subtract=[0,0], cmap=:bilbao, min_size=3, max_size=10,
                            mode="capacity", scale=:log10, sizzle=(2000,1000), fontfac=1, lw=1, dpi=100,
                            lens=false, lens_x1=25, lens_x2=900, lens_y1=18, lens_y2=900, reg=true)
    
    IJulia.clear_output()
    
    directory = readdir(path)
    filter!(e->occursin("csv",e), directory)
    
    better_directory = copy(directory)
    for i in 1:length(directory)
        better_directory[i] = directory[i][13:end-4]
    end
    
    allocation_dict = Dict()

    for allocation in better_directory
        
        thing = path*"/summary_for_"*allocation*".csv"
        table = CSV.read(thing, DataFrame, header=true)
        
        table."delta_capacity" .-=subtract[1]
        table."tau_capacity" .-=subtract[2]
        
        allocation_dict[allocation] = table
    end
    
    shapes = [:square, :utriangle, :diamond, :dtriangle, :circle, :star5]
    strats = better_directory
    
    dlabel = "\$\\delta\$-based resilience $(mode) \$R_\\delta\$"
    tlabel = "\$\\tau\$-based resilience $(mode) \$R_\\tau\$"
    if scale != :identity
        dlabel = "\$\\delta\$-based resilience $(mode) \$R_\\delta\$ (logarithmic scale)"
        tlabel = "\$\\tau\$-based resilience $(mode) \$R_\\tau\$ (logarithmic scale)"
    end

    p = plot(scale=scale, xlabel=dlabel, ylabel=tlabel, legend=:bottomright, size=sizzle, framestyle=:box,
             dpi=dpi, legendfontsize=8*fontfac, guidefontsize=11*fontfac, tickfontsize=8*fontfac, gridlinewidth=lw)

    low, high = 980.1, 0

    j = 0
    for allocation in sort(strats)

        j+=1
        allocation_index = j
        table = allocation_dict[allocation]

        deltas = table."delta_capacity"
        taus = table."tau_capacity"

        if mode == "probability"
            deltas = deltas /= 980.1
            taus = taus /= 980.1
        end

        if allocation != "no_adaptation"

            xlabel = names(table)[1]
            ylabel = names(table)[2]

            x = table[:,xlabel]
            y = table[:,ylabel]
            
            perm = sortperm(x, rev=true)
            x = x[perm]
            y = y[perm]
            deltas = deltas[perm]
            taus = taus[perm]
            
            if scale != :identity
                L = length(deltas)
                for i in 0:L-1
                    if (deltas[L-i] <= 0) || (taus[L-i] <= 0)
                        deleteat!(deltas, L-i)
                        deleteat!(taus, L-i)
                        deleteat!(x, L-i)
                        deleteat!(y, L-i)
                    end
                end
            end
            
            low = minimum([low, minimum(deltas), minimum(taus)])
            high = maximum([high, maximum(deltas), maximum(taus)])
            
            y_max = maximum(y)
            col = [(cgrad(cmap)[(y[i]/y_max)]) for i in 1:length(y)]

            log_x = log10.(x)
            x_max = maximum(log_x)
            x_min = minimum(log_x)
            x_span = x_max - x_min               
            size_diff = max_size-min_size
            size = [(min_size + size_diff*( (log_x[i]-x_min)/x_span )) for i in 1:length(x)]

            scatter!(deltas, taus, label=texify(allocation), markersize=size, markercolor=col, markerstrokecolor=:grey,
                     markerstrokewidth=1, markershape=shapes[allocation_index])
        else
            low = minimum([low, minimum(deltas), minimum(taus)])
            high = maximum([high, maximum(deltas), maximum(taus)])
            scatter!(deltas, taus, label=allocation, markersize=min_size, markercolor=:black, markerstrokecolor=:grey,
                     markerstrokewidth=1, markershape=:x)
        end 
    end
    
     plot!(smooth=reg)
    
    duo = [low, high]
    plot!(duo, duo, linecolor=:grey, linealpha=0.5, linestyle=:dot, label="identity line", linewidth=lw)
    
    if lens == true
        lens!([lens_x1,lens_x2], [lens_y1,lens_y2], inset = (1, bbox(0.05, 0.05, 0.48, 0.49)),
              scale=scale, gridlinewidth=lw, tickfontsize=8*fontfac)
    end
      
    p
end

# function that visualizes the fraction of suboptimal resilience values that occur without adaptaion
# depending on the safety factor and number of flowtest runs which determine the initial power line capacities
function status_quo_visualizer(path)
    
    c = readdir(path)
    coordinate_list = [c[3:4];c[1:2];c[7:8];c[5:6];c[11:12];c[9:10,];c[15:16];c[13:14]] 
    
    safeties = []
    flowtest_runs = []
    fail_fracs = []
    
    @showprogress for coordinate in coordinate_list
        
        
        thing = path*"/"*coordinate*"/frac.csv"
        table = CSV.read(thing, DataFrame, header=false, limit=1, drop=[2])
        frac = table[1,1]
        
        append!(fail_fracs,frac)
        
        s = value_from_filename(coordinate, "safety")
        ftr = value_from_filename(coordinate, "flowtest_runs")
        
        append!(safeties, s)
        append!(flowtest_runs, ftr)
    end

    safeties = unique(safeties)
    flowtest_runs = unique(flowtest_runs)

    LS = length(safeties)
    LF = length(flowtest_runs)

    fail_fracs = reshape(fail_fracs, (LF,LS))       
    
    @manipulate for scale in ["linear", "log"],
                    b = button("Save as .PNG"; value=0)
        
        if scale == "linear"
            zlabel="fraction of sustainant values < 1"
            title = "Status quo test: \n safety vs. flowtest_runs"
            plot_name = "Status_quo_test_-_safety_vs_flowtest_runs"
            z = fail_fracs
        elseif scale == "log"
            zlabel="log (fraction of sustainant values < 1)"
            title = "Status quo test: \n safety vs. flowtest_runs \n (logarithmic)"
            plot_name = "Status_quo_test_-_safety_vs_flowtest_runs_(logarithmic)"
            z = log10.(fail_fracs)
        end
        
        p = contour_marginals(safeties, flowtest_runs, z;
                              xlabel="line_safety", ylabel="flowtest_runs", zlabel=zlabel ,
                              title=title)

        if b > 0
            savefig("Plots/"*plot_name)
            b = 0
        end

        p
    end
end
