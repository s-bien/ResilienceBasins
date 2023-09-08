This code is in Julia script (.jl) form for importing it into other scripts.
The same code is available as a Jupyter Notebook (.ipynb) file for better readability.

using ProgressMeter
using Distributed

# function that divides a time series into chunks and saves them as a dictionary
function data2dict(data,days)
    dict = Dict()
    L = Int(floor(length(data)/days))
    for d in 1:days
        dict[d] = data[(d-1)*L+1:d*L]
    end
    return dict
end

# function that returns keys of a dictionary in a sorted order (rather than random)
function goodkeys(dict)
    list = []
    for key in keys(dict)
        append!(list,key)
    end
    return sort(list)
end

# function that returns the values of a dictionary in the corresponding order
function goodvals(dict)
    list = []
    for key in goodkeys(dict)
        append!(list,dict[key])
    end
    return list
end

# function that resamples a time series array to a new resolution
function merger(array,resolution)
    new_array = []
    bin = 0
    count = 0
    for i in array
        bin += i/resolution
        count += 1
        if count == resolution
            append!(new_array,bin)
            bin = 0
            count = 0
        end
    end
    return new_array
end

# function that applies the resampling above to a whole dictionary of arrays
function smoother(dict,resolution)
    new_dict = Dict()
    for k in keys(dict)
        new_dict[k] = merger(dict[k],resolution)
    end
    return new_dict
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

# function to create a power grid-like graph
function create_grid(n, n0, p, q, r)
    s = 0.
    u = 1.
    return generate_graph( RandomPowerGrid(n, n0, p, q, r, s, u) )
end

# function to create a copy of a graph object (pointers will overwrite the original instance)
function graphcopy(graph)
    n = nv(graph)
    a = adjacency_matrix(graph)
    copy = SimpleGraph(n)
    for i in 1:n
        for j in i+1:n
            if a[i,j] !== 0
                add_edge!(copy,i,j)
            end
        end
    end
    return copy
end

# function to copy an embedded graph object
function embedded_graphcopy(eg)

    loc_x = vertices_loc(eg,1)
    loc_y = vertices_loc(eg,2)
    loc_vector = [[loc_x[i],loc_y[i]+1] for i in 1:length(loc_x)]
    return EmbeddedGraph(graphcopy(eg), loc_vector)
    
end

# function to get the edge object of a certain index within a graph
function index2edge(graph, index)
    i=1
    for e in edges(graph)
        if i == index
            return e
        else
            i+=1
        end
    end
end

# function to get the index number of a certain edge object
function edge2index(graph, src, dst)
    a=adjacency_matrix(graph)
    n = nv(graph)
    counter = 0
    for i in 1:src
        for j in i+1:n
            if a[i,j] !== 0
                counter += 1
                if i==src && j==dst
                    return counter
                end
            end
        end
    end
end

# struct with two fields for old and new states of a power grid
struct old_new_struct
    old
    new
end

# function that randomly allocates the prosumer identity to a given number of nodes
# two stages are returned, first with P1 and second with P2 prosumers
function set_prosumers(N,P1,P2,S)
    
    draw = sample(1:N-1,P2,replace=false)
    for i in 1:P2
        if draw[i] >= S
            draw[i] += 1
        end
    end
    vec1 = zeros(Int,N)
    vec2 = zeros(Int,N)
    
    for i in 1:P2
        
        if i <= P1
            vec1[draw[i]] = 1
        end
        vec2[draw[i]] = 1
             
    end
        
    return old_new_struct(vec1,vec2)
end

# function that selects a value from a random timeseries in a dictionary at a given time
function randomdraw(dict, time)
    day = rand(keys(dict))
    index = Int( ceil( time*length(dict[day]) ))
    return dict[day][index]
end

# struct with two fields for power demand and supply values
struct power_struct
    demand
    supply
end

# function to generate random step-wise time series of demand and supply for each node of a given power grid
# these are only used to calculate initial line capacities
function power_randomdraws(solar_map, solar_dict, consumer_dict, prosumer_ratio, timesteps, slack_index)
    
    NV = length(solar_map)
    consumption, production = zeros(timesteps,NV), zeros(timesteps,NV)
    
    for t in 1:timesteps
        
        time = rand()
        
        for n in 1:NV
            if n != slack_index
                consumption[t,n] = randomdraw(consumer_dict,time)
            end
            if solar_map[n] == 1
                production[t,n] = prosumer_ratio*randomdraw(solar_dict,time)
            end
        end
        
        #consumption[t,slack_index] = -sum(consumption[t,:])
        #production[t,slack_index] = -sum(production[t,:])
    end
    
    return power_struct(consumption,production)
end

# function to generate random chunk-wise time series of demand and supply for each node of a given power grid
# these are used in the actual simulation
function power_timeseries(solar_map, solar_dict, consumer_dict, prosumer_ratio, days, time_resolution, slack_index)
    
    NV = length(solar_map)
    timesteps = Int(ceil(24*60/time_resolution))
    
    L = days*timesteps
    consumption, production = zeros(L,NV), zeros(L,NV)
    
    for n in 1:NV
        for d in 1:days
            
            if n != slack_index
                c_draw = sample(1:consumption_days)
                consumption[ (d-1)*timesteps+1 : d*timesteps,n] = consumer_dict[c_draw]
                
                if solar_map[n] == 1
                    s_draw = sample(1:solar_days)
                    production[ (d-1)*timesteps+1 : d*timesteps,n] = prosumer_ratio*solar_dict[s_draw]
                end
            end
        end
    end
    
    #for t in 1:L
    #    consumption[t,slack_index] = -sum(consumption[t,:])
    #    production[t,slack_index] = -sum(production[t,:])
    #end
        
    return power_struct(consumption,production)
end

# function to compute the linear power flows on a given graph with given power injections
function flow(g, p)
    b = incidence_matrix(g, Float64, oriented=true)
    f = lsqr(b, p)
    @assert all(isfinite.(f))
    return -floor.(f, digits=16)
end

# function to extract the maximum magnitude of flows occuring on each edge given a time series of injections
function max_flows(grid, p_timeseries, slack_index)
    NE = ne(grid)
    len = size(p_timeseries.demand)[1]
    fseries = zeros(len,NE)
    
    injection = p_timeseries.supply-p_timeseries.demand
    
    for i in 1:len
        injection[i,slack_index] = -sum(injection[i,:])
        fseries[i,:] = abs.(flow(grid,injection[i,:]))
    end
    
    max = zeros(NE)
    for j in 1:NE
        max[j] = maximum(fseries[:,j])
    end
    
    return floor.(max, digits=16)
end

# function that corrects a given edge budget allocation
# such that the minimum non-zero allocation is large enough to prevent numerical errors
function zero_corrector(weights, nll)
    
    new_weights = floor.( copy(weights), digits=16)
    
    if sum(new_weights) == 0
        new_weights = ones(length(weights))
    else
        min_weight = minimum( filter(e->e!=0, weights) )
        @assert min_weight > 0

        if length(nll) > 0
            for e in nll
                if new_weights[e] == 0
                    new_weights[e] = min_weight
                end
            end
        end
    end
    
    new_weights = floor.( new_weights, digits=16)
    @assert any(new_weights.>0)
    return new_weights
end

# function to compute euclidean lenghts of edges in an embedded graph
function line_lengths(grid)
    
    lengths = zeros(ne(grid))
    
    for e in 1:ne(grid)
        edge = index2edge(grid,e)
        i = src(edge)
        j = dst(edge)
        lengths[e] = euclidean(grid.vertexpos[i], grid.vertexpos[j])
    end
    
    lengths = floor.(lengths, digits=16)
    @assert all(lengths.>0)
    return lengths
end

# function to generate edge weights proportional to endge length
# this corresponds to uniform capacity upgrades
function edge_weights_uniform(grid)
   
    weights = ones(ne(grid))
    #weight corrector unnecessary
    
    for e in 1:length(weights)
        edge = index2edge(grid,e)
        i = src(edge)
        j = dst(edge)
        weights[e] *= euclidean(grid.vertexpos[i], grid.vertexpos[j])
    end
    
    weights /= sum(weights)
    return floor.(weights,digits=16)
end

# function to compute pairwise graph distances between nodes of a graph
function gdist_matrix(grid)
    
    NV = nv(grid)
    m = zeros(NV, NV)
    
    for i in 1:NV
        m[i,:] = gdistances(grid,i)
    end
    
    return m
end

# function to generate edge weights proportional to their length times the sum of their closenesses to all prosumer nodes
# closeness is here defined as the distance raised to a negative exponent
function edge_weights_by_closeness(grid, map, exponent)
   
    gdm = gdist_matrix(grid)
    
    weights = zeros(ne(grid))
    for e in 1:length(weights)
        edge = index2edge(grid,e)
        i = src(edge)
        j = dst(edge)
        
        dist_list = [minimum([gdm[i,v],gdm[j,v]])+1 for v in 1:nv(grid)]
        weights[e] = sum(map .* (dist_list.^exponent) ) * euclidean(grid.vertexpos[i], grid.vertexpos[j])
    end
    
    #weights = zero_corrector(weights, [])
    weights /= sum(weights)
    return floor.(weights,digits=16)
end

# function that computes the exceedings of power flows over line capacities
# based on a gives power injection snapshot
function count_violations(grid, caps, power)
    
    diff = abs.(flow(grid,power.supply-power.demand)) - caps
    
    for i in 1:length(diff)
        if diff[i] < 0
            diff[i] = 0
        end
    end
    
    return floor.(diff, digits=16)
end

# function that computes the maximum occuring exceeding of power flows on each line
function max_violations(grid, caps, p_timeseries)
    
    max_vio = zeros(length(caps))
    timesteps = size(p_timeseries.demand)[1]
    
    for t in 1:timesteps
        
        vio = count_violations(grid, caps, power_struct(p_timeseries.demand[t,:],p_timeseries.supply[t,:]))
        for e in 1:length(caps)
            
            max_vio[e] = maximum([max_vio[e],vio[e]])
        end
    end
        
    return floor.(max_vio, digits=16)
end

# function that computes the edge weights based on line length and maximum power exceedings
function edge_weights_by_violation(grid, caps, p_timeseries, nll)
   
    violations = max_violations(grid, caps, p_timeseries)
    
    weights = zero_corrector(violations, nll)
    
    for e in 1:length(weights)
        edge = index2edge(grid,e)
        i = src(edge)
        j = dst(edge)
        weights[e] *= euclidean(grid.vertexpos[i], grid.vertexpos[j])
    end
    
    weights /= sum(weights)
    return floor.(weights,digits=16)
end

# function that computes node weights based on the sum of exceeding-based edge weights of its incident edges
function node_weights_by_violation(grid, caps, p_timeseries)
    
    max_vio = max_violations(grid, caps, p_timeseries)    
    NV = nv(grid)
    weights = zeros(NV)
    
    for i in 1:length(caps)
        e = index2edge(grid,i)
        s,d = src(e), dst(e)
        weights[s] += max_vio[i]
        weights[d] += max_vio[i]
    end
    
    weights = zero_corrector(weights,[])
    
    weights /= sum(weights)
    return floor.(weights,digits=16)
end

# function that finds the pair of non-neighboring nodes with the largest redundancy effect if connected by a new edge
function max_redundancy_pair(graph, r)
    
    path_length = zeros(nv(graph), nv(graph))
    air_dist = zeros(nv(graph), nv(graph))
    
    for i in 1:nv(graph)
        path_length[i,:] = dijkstra_shortest_paths(graph,i).dists
        air_dist[i,:] = map(j -> euclidean(graph.vertexpos[i], graph.vertexpos[j]), 1:nv(graph))
    end
    
    rho = ((path_length .+ air_dist) .^ r) ./ air_dist
    for i in 1:nv(graph)
        rho[i,i] = 0
    end
    
    arg = argmax(rho)
    while has_edge(graph,arg[1],arg[2]) == true
        rho[arg[1],arg[2]] = 0
        rho[arg[2],arg[1]] = 0
        arg = argmax(rho)
    end
    
    return sort([arg[1],arg[2]])
end

# struct for returning an altered version of a grid
struct newline_struct
    grid
    caps
end

# function that adds a single line to an existing grid in the most redundant location
function add_line(graph, caps, r)
    
    pair = max_redundancy_pair(graph,r)
    new_grid = embedded_graphcopy(graph)
    add_edge!(new_grid, pair[1], pair[2])
    
    
    index = edge2index(new_grid, pair[1], pair[2])
    new_caps = zeros(length(caps)+1)
    for i in 1:length(caps)
        if i < index
            new_caps[i] = caps[i]
        else
            new_caps[i+1] = caps[i]
        end
    end
    
    return newline_struct(new_grid, new_caps)
end

# struct to keep track of battery states
struct battery_struct
    demand
    supply
    battery_content #at end of time step
    battery_use
    battery_emergency_use
    injection
    waste
end

# function that calculates the power injections of a grid with batteries
# based on current demand and supply, previous battery charging levels, and the battery reliance parameter
# "emergency_mode" is a deprecated argument
function batterize(grid, capacities, demand_supply, battery_map, content_old, reliance, slack_index, emergency_mode)
    
    tolerance = 0 #10^-16 #?
    
    g = graphcopy(grid.graph)
    caps = copy(capacities)
    
    demand = copy(demand_supply.demand)
    supply = copy(demand_supply.supply)
    demand[slack_index] = 0
    supply[slack_index] = 0
    
    NV = length(demand)
    
    battery_content = zeros(NV) #at end of time step
    battery_use = zeros(NV)
    battery_emergency_use = zeros(NV)
    injection = zeros(NV)
    waste = zeros(NV)
     
    for node in 1:NV
        if battery_map[node] > 0+tolerance
            if supply[node]-demand[node] > 0+tolerance
                charge_capacity = (battery_map[node]-content_old[node])
                battery_use[node] = minimum([charge_capacity, reliance*abs(supply[node]-demand[node])])
            elseif supply[node]-demand[node] < 0-tolerance
                charge_capacity = content_old[node]
                battery_use[node] = - minimum([charge_capacity, reliance*abs(supply[node]-demand[node])])
            end
        end
    end
    
    injection = supply - demand - battery_use
    injection[slack_index] = -sum(injection)
    
    diff = caps - abs.(flow(g, injection))
    
    while minimum(diff) < 0-tolerance
        
        indices = findall(x -> x<0-tolerance, diff)
        sort!(indices,rev=true)
        for i in indices
            e = index2edge(g,i)
            rem_edge!(g, e)
            deleteat!(caps,i)
        end
        
        cc = connected_components(g)            
            
        injection[slack_index] = 0
        for c in cc
            subgraph_balance = 0
            for node in c
                subgraph_balance += injection[node]
            end
            subgraph_balance = round(subgraph_balance, digits=14)
            
            if slack_index in c
                injection[slack_index] = -subgraph_balance
                
            else
                if subgraph_balance > 0+tolerance
                    
                    subgraph_battery_vacancy = 0
                                        
                    if emergency_mode == true
                        for node in c
                        subgraph_battery_vacancy += (battery_map[node] - content_old[node] 
                                                     - battery_use[node] - battery_emergency_use[node] )
                        end
                        subgraph_battery_vacancy = round(subgraph_battery_vacancy, digits=14)
                        
                        for node in c
                            if subgraph_battery_vacancy > 0
                                charge = (battery_map[node] - content_old[node]
                                          - battery_use[node] - battery_emergency_use[node] )*
                                          minimum([1, subgraph_balance/subgraph_battery_vacancy])
                            else
                                charge = 0
                            end
                            battery_emergency_use[node] += charge
                            injection[node] -= charge
                        end
                    end
                    
                    subgraph_infeed = 0
                    for node in c
                        if injection[node] > 0+tolerance
                            subgraph_infeed += injection[node]
                        end
                    end
                    subgraph_infeed = round(subgraph_infeed, digits=14)
                    
                    if subgraph_infeed > 0
                        wasting_factor = minimum(
                                         [1, maximum([0, (subgraph_balance-subgraph_battery_vacancy)/subgraph_infeed ]) ])
                    else
                        wasting_factor = 0
                    end
                    
                    @assert wasting_factor <= 1+10^-14 "wasting factor > 1: w=$(wasting_factor), J=$(subgraph_balance), V=$(subgraph_battery_vacancy), J+=$(subgraph_infeed)"
                    
                    for node in c
                        if injection[node] > 0+tolerance
                            waste[node] += injection[node] * wasting_factor
                            injection[node] -= injection[node] * wasting_factor
                        end
                    end
                        
                elseif subgraph_balance < 0-tolerance
                    
                    if emergency_mode == true
                        
                        subgraph_battery_content = 0
                        for node in c
                        subgraph_battery_content += (content_old[node] +
                                                     battery_use[node] + battery_emergency_use[node] )
                        end
                        subgraph_battery_content = round(subgraph_battery_content, digits=14)
                        
                        if subgraph_battery_content >= abs(subgraph_balance)
                        
                            for node in c
                                if subgraph_battery_content > 0
                                    charge = (content_old[node] + battery_use[node] + battery_emergency_use[node] )*
                                              minimum([1, (abs(subgraph_balance)/subgraph_battery_content) ])
                                else
                                    charge = 0
                                end
                                battery_emergency_use[node] -= charge
                                injection[node] += charge
                            end
                        end
                    else
                        for node in c
                            
                            if emergency_mode == true
                                
                                if injection[node] > 0+tolerance
                                    charge_capacity = (battery_map[node] - content_old[node] -
                                                       battery_use[node] - battery_emergency_use[node])
                                    charge = minimum([charge_capacity, injection[node]])
                                    battery_emergency_use[node] += charge
                                    injection[node] -= charge
                                    
                                elseif injection[node] < 0-tolerance
                                    charge_capacity = content_old[node] +
                                                      battery_use[node] + battery_emergency_use[node]
                                    charge = minimum([charge_capacity, abs(injection[node])])
                                    battery_emergency_use[node] -= charge
                                    injection[node] += charge
                                end
                            end
                            
                            for node in c
                                waste[node] += maximum([0, injection[node]])
                                injection[node] = 0
                            end
                        end
                    end
                end
            end
        end
        
        if ne(g) == 0
            break
        else
            diff = caps - abs.(flow(g,injection))
        end
    end
    
    battery_content = content_old + battery_use + battery_emergency_use
    
    for node in 1:NV
        if (supply[node] > demand[node])
            waste[node] = minimum([waste[node],supply[node]-demand[node]])
        end
    end
    
    power = battery_struct(demand, supply, battery_content, battery_use, battery_emergency_use,
                           injection, waste)
    
    return power
end

# deprecated function for calculating battery behavior in "emergency_mode"
function batterize_old_em(grid, capacities, demand_supply, battery_map, content_old, reliance, slack_index, emergency_mode)
    
    tolerance = 0 #10^-16 #?
    
    sigma_solar = 1.9908019172377822
    sigma_consumption = 0.6017518938044338
    
    g = graphcopy(grid.graph)
    caps = copy(capacities)
    
    demand = copy(demand_supply.demand)
    supply = copy(demand_supply.supply)
    demand[slack_index] = 0
    supply[slack_index] = 0
    
    NV = length(demand)
    
    battery_content = zeros(NV) #at end of time step
    battery_use = zeros(NV)
    battery_emergency_use = zeros(NV)
    injection = zeros(NV)
    waste = zeros(NV)
     
    for node in 1:NV
        if battery_map[node] > 0+tolerance
            if supply[node]-demand[node] > 0+tolerance
                charge_capacity = (battery_map[node]-content_old[node])
                battery_use[node] = minimum([charge_capacity, reliance*abs(supply[node]-demand[node])])
            elseif supply[node]-demand[node] < 0-tolerance
                charge_capacity = content_old[node]
                battery_use[node] = - minimum([charge_capacity, reliance*abs(supply[node]-demand[node])])
            end
        end
    end
    
    injection = supply - demand - battery_use
    injection[slack_index] = -sum(injection)
    
    diff = caps - abs.(flow(g, injection))
    
    while minimum(diff) < 0-tolerance
        
        indices = findall(x -> x<0-tolerance, diff)
        sort!(indices,rev=true)
        for i in indices
            e = index2edge(g,i)
            rem_edge!(g, e)
            deleteat!(caps,i)
        end
        
        cc = connected_components(g)            
            
        injection[slack_index] = 0
        for c in cc
            subgraph_balance = 0
            for node in c
                subgraph_balance += injection[node]
            end
            subgraph_balance = round(subgraph_balance, digits=14)
            
            if slack_index in c
                injection[slack_index] = -subgraph_balance
                
            else
                if subgraph_balance > 0+tolerance
                    
                    subgraph_battery_vacancy = 0
                                        
                    if emergency_mode == true
                        for node in c
                        subgraph_battery_vacancy += (battery_map[node] - content_old[node] 
                                                     - battery_use[node] - battery_emergency_use[node] )
                        end
                        subgraph_battery_vacancy = round(subgraph_battery_vacancy, digits=14)
                        
                        for node in c
                            if subgraph_battery_vacancy > 0
                                charge = (battery_map[node] - content_old[node]
                                          - battery_use[node] - battery_emergency_use[node] )*
                                          minimum([1, subgraph_balance/subgraph_battery_vacancy])
                            else
                                charge = 0
                            end
                            battery_emergency_use[node] += charge
                            injection[node] -= charge
                        end
                    end
                    
                    subgraph_infeed = 0
                    for node in c
                        if injection[node] > 0+tolerance
                            subgraph_infeed += injection[node]
                        end
                    end
                    subgraph_infeed = round(subgraph_infeed, digits=14)
                    
                    if subgraph_infeed > 0
                        wasting_factor = minimum(
                                         [1, maximum([0, (subgraph_balance-subgraph_battery_vacancy)/subgraph_infeed ]) ])
                    else
                        wasting_factor = 0
                    end
                    
                    @assert wasting_factor <= 1+10^-14 "wasting factor > 1: w=$(wasting_factor), J=$(subgraph_balance), V=$(subgraph_battery_vacancy), J+=$(subgraph_infeed)"
                    
                    for node in c
                        if injection[node] > 0+tolerance
                            waste[node] += injection[node] * wasting_factor
                            injection[node] -= injection[node] * wasting_factor
                        end
                    end
                        
                elseif subgraph_balance < 0-tolerance
                    
                    if emergency_mode == true
                        
                        subgraph_battery_content = 0
                        for node in c
                        subgraph_battery_content += (content_old[node] +
                                                     battery_use[node] + battery_emergency_use[node] )
                        end
                        subgraph_battery_content = round(subgraph_battery_content, digits=14)
                        
                        if subgraph_battery_content >= abs(subgraph_balance)
                        
                            for node in c
                                if subgraph_battery_content > 0
                                    charge = (content_old[node] + battery_use[node] + battery_emergency_use[node] )*
                                              minimum([1, (abs(subgraph_balance)/subgraph_battery_content) ])
                                else
                                    charge = 0
                                end
                                battery_emergency_use[node] -= charge
                                injection[node] += charge
                            end
                        end
                    else
                        for node in c
                            
                            if emergency_mode == true
                                
                                if injection[node] > 0+tolerance
                                    charge_capacity = (battery_map[node] - content_old[node] -
                                                       battery_use[node] - battery_emergency_use[node])
                                    charge = minimum([charge_capacity, injection[node]])
                                    battery_emergency_use[node] += charge
                                    injection[node] -= charge
                                    
                                elseif injection[node] < 0-tolerance
                                    charge_capacity = content_old[node] +
                                                      battery_use[node] + battery_emergency_use[node]
                                    charge = minimum([charge_capacity, abs(injection[node])])
                                    battery_emergency_use[node] -= charge
                                    injection[node] += charge
                                end
                            end
                            
                            for node in c
                                waste[node] += maximum([0, injection[node]])
                                injection[node] = 0
                            end
                        end
                    end
                end
            end
        end
        
        if ne(g) == 0
            break
        else
            diff = caps - abs.(flow(g,injection))
        end
    end
    
    battery_content = content_old + battery_use + battery_emergency_use
    
    for node in 1:NV
        if (supply[node] > demand[node])
            waste[node] = minimum([waste[node],supply[node]-demand[node]])
        end
    end
    
    power = battery_struct(demand, supply, battery_content, battery_use, battery_emergency_use,
                           injection, waste)
    
    return power
end

# function to compute the power injection and battery level time series
# from a time series of supply and demand, intitial battery charge, and the reliance parameter
function batterize_series(grid, caps, demand_supply_series, battery_map, initial_content,
                          reliance, slack_index, emergency_mode)
       
    timesteps, NV = size(demand_supply_series.demand)
    
    demand = demand_supply_series.demand
    supply = demand_supply_series.supply
    battery_content = zeros(timesteps, NV)
    battery_use = zeros(timesteps, NV)
    battery_emergency_use = zeros(timesteps, NV)
    injection = zeros(timesteps, NV)
    waste = zeros(timesteps, NV)
    
    content_0 = battery_map .* initial_content

    batt_0 = batterize(grid, caps, power_struct(demand[1,:],supply[1,:]), battery_map, content_0,
                       reliance, slack_index, emergency_mode)
    battery_content[1,:] = batt_0.battery_content
    battery_use[1,:] = batt_0.battery_use
    battery_emergency_use[1,:] = batt_0.battery_emergency_use
    injection[1,:] = batt_0.injection
    waste[1,:] = batt_0.waste

    for t in 2:timesteps
        batt = batterize(grid, caps, power_struct(demand[t,:],supply[t,:]), battery_map, battery_content[t-1,:],
                         reliance, slack_index, emergency_mode)
        battery_content[t,:] = batt.battery_content
        battery_use[t,:] = batt.battery_use
        battery_emergency_use[t,:] = batt.battery_emergency_use
        injection[t,:] = batt.injection
        waste[t,:] = batt.waste
    end
    
    return battery_struct(demand, supply, battery_content, battery_use, battery_emergency_use,
                          injection, waste)   
end

# deprecated version of the above battery struct
struct battery_struct_old
    charge #resulting charge at end of timestep
    demand
    supply #supply visible to the grid
end

# deprecated version of the above battery function
function batterize_old(battery_map, charge, reliance, power, slack_index)
    
    NV = length(charge)
    new_charge, visible_supply = zeros(NV), zeros(NV)
    
    for i in 1:NV
        
        if battery_map[i] > 0

            if power.supply[i]-power.demand[i] > 0
                charge_budget = reliance*(battery_map[i]-charge[i])
                diff = minimum([charge_budget, abs(power.supply[i]-power.demand[i])])
                visible_supply[i] = power.supply[i] - diff
                new_charge[i] = charge[i] + diff
            elseif power.supply[i]-power.demand[i] < 0
                charge_budget = reliance*charge[i]
                diff = minimum([charge_budget, abs(power.supply[i]-power.demand[i]) ])
                visible_supply[i] = power.supply[i] + diff
                new_charge[i] = charge[i] - diff
            elseif power.supply[i]-power.demand[i] == 0
                visible_supply[i] = power.supply[i]
                new_charge[i] = charge[i]
            end
            
        else 
            if i != slack_index
                visible_supply[i] = power.supply[i]
            end
        end
    end

    visible_supply[slack_index] = -sum(visible_supply)
    
    return battery_struct(new_charge, power.demand, visible_supply)
end

# deprecated version of the above battery function
function batterize_series_old(battery_map, initial_content, reliance, power_series, slack_index)
    
    p = power_struct(copy(power_series.demand),copy(power_series.supply))
    timesteps, NV = size(p.demand)
    
    charge = zeros(timesteps, NV)
    charge_0 = battery_map .* initial_content

    batt_0 = batterize(battery_map, charge_0, reliance, power_struct(p.demand[1,:],p.supply[1,:]), slack_index)
    charge[1,:] = batt_0.charge
    p.supply[1,:] = batt_0.supply

    for t in 2:timesteps
        batt = batterize(battery_map, charge[t-1,:], reliance, power_struct(p.demand[t,:],p.supply[t,:]), slack_index)
        charge[t,:] = batt.charge
        p.supply[t,:] = batt.supply
    end

    return battery_struct(charge, p.demand, p.supply)    
end

# struct that contains two different sustainant values
# (only one of the two was used in the final article)
struct sustainant_struct
    delta
    tau
end

# function that computes two different sustainant values at a single point in time based on the current power injections
# (only the one called tau was used in the final article)
function sustainants(power, slack_index)
    
    tolerance = 10^-14 #0
    
    power.demand[slack_index] = 0
    power.supply[slack_index] = 0
    power.battery_content[slack_index] = 0
    power.battery_use[slack_index] = 0
    power.battery_emergency_use[slack_index] = 0
    power.injection[slack_index] = 0
    power.waste[slack_index] = 0
    
    lack = -power.supply +power.demand +power.battery_use +power.battery_emergency_use +power.injection +power.waste
    
    served_demand = power.demand - lack #power.demand + mismatch.*Int.(mismatch.<0)
    
    if sum(power.demand) != 0
        delta = sum(served_demand)/sum(power.demand)
    else
        delta = 1
    end
    
    if -tolerance <= delta <0
        delta = 0
    elseif 1 < delta <= 1+tolerance
        delta = 1
    end
    
    @assert delta >= 0 "negative delta: $(delta)"
    @assert delta <= 1 "delta above 1: $(delta)"
    
    damage = lack + power.waste
    
    #aka served transmission, served injection?
    served_mismatch = abs.(power.supply-power.demand) - abs.(damage) # damage ohne abs?
    
    if sum(abs.(power.supply-power.demand)) != 0
        tau = sum(served_mismatch)/sum(abs.(power.supply-power.demand))
    else
        tau = 1
    end
    
    if -tolerance <= tau <0
        tau = 0
    elseif 1 < tau <= 1+tolerance
        tau = 1
    end
    
    @assert tau >= 0 "negative tau: $(tau)"
    @assert tau <= 1 "tau above 1: $(tau)"
    
    return sustainant_struct(delta, tau)
end

# function that computes the sustainant time series from the power injection time series
function sustainants_series(power_series, slack_index)
    
    delta, tau = [], []
    
    for i in 1:size(power_series.demand)[1]
        
        d = power_series.demand[i,:]
        s = power_series.supply[i,:]
        bc = power_series.battery_content[i,:]
        bu = power_series.battery_use[i,:]
        beu = power_series.battery_emergency_use[i,:]
        j = power_series.injection[i,:]
        w = power_series.waste[i,:]
        p = battery_struct(d,s,bc,bu,beu,j,w)

        sus = sustainants(p, slack_index)
        append!(delta,sus.delta)
        append!(tau,sus.tau)
    end
    
    return sustainant_struct(delta,tau)
end

# deprecated version of the above sustainant function
function sustainant_old(grid, capacities, power, slack_index)
    
    tolerance = 0 #10^-16
    
    g = graphcopy(grid.graph)
    caps = copy(capacities)
    
    d = copy(power.demand)
    d[slack_index] = 0

    b = copy(power.supply-power.demand)
    b[slack_index] = 0

    total_demand = sum(d)
    total_balance = sum(abs.(b))

    d_realized = copy(power.demand)
    b_realized = copy(power.supply-power.demand)

    diff = caps - abs.(flow(g,b_realized))
    
    while minimum(diff) < 0-tolerance
        
        indices = findall(x-> x<0-tolerance, diff)
        sort!(indices,rev=true)
        for i in indices
            e = index2edge(g,i)
            rem_edge!(g, e)
            deleteat!(caps,i)
        end
        
        cc = connected_components(g)
        
        b_realized[slack_index] = 0
        for c in cc
            subgraph_balance = 0
            for i in c
                subgraph_balance += b_realized[i] 
            end
            
            if slack_index in c
                b_realized[slack_index] = -subgraph_balance
            
            else
                if subgraph_balance < 0-tolerance
                    for i in c
                        b_realized[i] = 0
                        d_realized[i] = 0 #better: d_realized[i] = power.supply[i]
                    end
                    
                elseif subgraph_balance > 0+tolerance
                    subgraph_surplus = 0
                    for i in c
                        if b_realized[i] > 0+tolerance
                            subgraph_surplus += b_realized[i]
                        end
                    end

                    shrinking_factor = (subgraph_surplus-subgraph_balance)/subgraph_surplus

                    for i in c
                        if b_realized[i] > 0+tolerance
                            b_realized[i] *= shrinking_factor
                        end
                    end
                end
            end
        end
        
        if ne(g) == 0
            break
        else
            diff = caps - abs.(flow(g,b_realized))
        end
    end
    
    @assert all(isfinite.(d))
    @assert all(isfinite.(b))
    @assert all(isfinite.(d_realized))
    @assert all(isfinite.(b_realized))

    d_realized[slack_index] = 0
    b_realized[slack_index] = 0
    total_demand_realized = sum(d_realized)
    total_balance_realized = sum(abs.(b_realized))
    
    delta = 0.
    tau = 0.

    if total_demand == 0
        delta = 1
    else
        delta = total_demand_realized / total_demand
    end

    if total_balance == 0
        tau = 1
    else
        tau = total_balance_realized / total_balance
    end    

    return sustainant_struct(delta,tau)
end

# deprecated version of the above sustainant function
function sustainant_series_old(grid, capacities, power_series, slack_index)
    
    delta, tau = [], []
    
    for i in 1:size(power_series.demand)[1]
        
        d = power_series.demand[i,:]
        s = power_series.supply[i,:]
        p = power_struct(d,s)

        sus = sustainant(grid, capacities, p, slack_index)
        append!(delta,sus.delta)
        append!(tau,sus.tau)
    end
    
    return sustainant_struct(delta,tau)
end

# struct to contain an upper and lower limit
struct lower_upper_struct
    lower
    upper
end

# function that calculates the upper and lower limit of injections that are allowed by the battery charging algorithm
# higher and lower injections will be partially absorbed by the battery
function limit_finder(rp, np, reliance; solar_dict=solar_dict, consumer_dict=consumer_dict)
    
    j_pop = combo_hist(consumer_dict, solar_dict; ratio=rp).list
    j_mean = rp-1
    
    goal_frac = minimum([1, reliance*99*10/rp/np]) #rp dependence is approximate and bad at low rp
    
    max_cum_ex = sum(filter(x -> x>0, j_pop.-j_mean))
    
    tolerance = 10^-5
    
    lim = j_mean
    rim = j_pop[end]
    upper_limit = j_pop[end]
    old_frac = -Inf
    frac = 0
    while abs(frac-old_frac) > tolerance && abs(frac-goal_frac) > tolerance
        
        upper_limit = mean([lim, rim])
        cum_ex = sum(filter(x -> x>0, j_pop.-upper_limit))
        old_frac = copy(frac)
        frac = cum_ex/max_cum_ex
        
        if frac < goal_frac
            rim = copy(upper_limit)
        elseif frac > goal_frac
            lim = copy(upper_limit)
        end
    end
    
    uppercut_j_pop = replace(x -> x>upper_limit ? upper_limit : x, j_pop)
    
    lim = j_pop[1]
    rim = j_mean   
    lower_limit = -j_pop[1]
    old_man = -1000
    man = -j_pop[1]
    while abs(man-old_man) > tolerance && abs(man-j_mean) > tolerance
        
        lower_limit = mean([lim, rim])
        old_man = copy(man)
        man = mean(replace(x -> x<lower_limit ? lower_limit : x, uppercut_j_pop))
        
        if man < j_mean
            lim = copy(lower_limit)
        elseif man > j_mean
            rim = copy(lower_limit)
        end
    end
    
    println("rp: ", rp)
    println("np: ", np)
    println("reliance: ", reliance)
    println("j_mean: ", j_mean)
    println("goal_frac: ", goal_frac) 
    println("upper_limit: ", upper_limit)
    println("lower_limit: ", lower_limit)
    println()
    
    return lower_upper_struct(lower_limit, upper_limit)
end

# struct to contain two sunstainant values and the battery content distribution
struct triple_struct
    delta
    tau
    battery_content
end

# function that saves memory by performing the battery algorithm and the sustainant calculation at the same time
# this avoids sending large arrays of injections and battery content between functions
function batt_sus(grid, capacities, demand, supply, battery_map, content_old, reliance, slack_index,
                  emergency_mode, limits)
    
    tolerance = 0 #10^-16 #?
    
    g = graphcopy(grid.graph)
    caps = copy(capacities)
    
    demand[slack_index] = 0
    supply[slack_index] = 0
    
    NV = length(demand)
    
    battery_content = zeros(NV) #at end of time step
    battery_use = zeros(NV)
    battery_emergency_use = zeros(NV)
    injection = zeros(NV)
    waste = zeros(NV)
     
    for node in 1:NV
        if battery_map[node] > 0+tolerance
            if supply[node]-demand[node] > limits.upper
                charge_capacity = (battery_map[node]-content_old[node])
                charge_goal = supply[node]-demand[node] - limits.upper
                battery_use[node] = minimum([charge_capacity, charge_goal])
                
            elseif supply[node]-demand[node] < limits.lower
                charge_capacity = content_old[node]
                charge_goal = limits.lower - (supply[node]-demand[node])
                battery_use[node] = -minimum([charge_capacity, charge_goal])
            end
        end
    end
    
    injection_0 = supply - demand - battery_use
    injection = supply - demand - battery_use
    injection[slack_index] = -sum(injection)
    
    diff = caps - abs.(flow(g, injection))
    
    while minimum(diff) < 0-tolerance
        
        indices = findall(x -> x<0-tolerance, diff)
        sort!(indices,rev=true)
        for i in indices
            e = index2edge(g,i)
            rem_edge!(g, e)
            deleteat!(caps,i)
        end
        
        cc = connected_components(g)            
            
        injection[slack_index] = 0
        for c in cc
            subgraph_balance = 0
            for node in c
                subgraph_balance += injection[node]
            end
            subgraph_balance = round(subgraph_balance, digits=14)
            
            if slack_index in c
                injection[slack_index] = -subgraph_balance
                
            else
                if subgraph_balance > 0+tolerance
                    
                    subgraph_battery_vacancy = 0
                                        
                    if emergency_mode == true
                        for node in c
                        subgraph_battery_vacancy += (battery_map[node] - content_old[node] 
                                                     - battery_use[node] - battery_emergency_use[node] )
                        end
                        subgraph_battery_vacancy = round(subgraph_battery_vacancy, digits=14)
                        
                        for node in c
                            if subgraph_battery_vacancy > 0
                                charge = (battery_map[node] - content_old[node]
                                          - battery_use[node] - battery_emergency_use[node] )*
                                          minimum([1, subgraph_balance/subgraph_battery_vacancy])
                            else
                                charge = 0
                            end
                            battery_emergency_use[node] += charge
                            injection[node] -= charge
                        end
                    end
                    
                    subgraph_infeed = 0
                    for node in c
                        if injection[node] > 0+tolerance
                            subgraph_infeed += injection[node]
                        end
                    end
                    subgraph_infeed = round(subgraph_infeed, digits=14)
                    
                    if subgraph_infeed > 0
                        wasting_factor = minimum(
                                         [1, maximum([0, (subgraph_balance-subgraph_battery_vacancy)/subgraph_infeed ]) ])
                    else
                        wasting_factor = 0
                    end
                    
                    @assert wasting_factor <= 1+10^-14 "wasting factor > 1: w=$(wasting_factor), "*
                            "J=$(subgraph_balance), V=$(subgraph_battery_vacancy), J+=$(subgraph_infeed)"
                    
                    for node in c
                        if injection[node] > 0+tolerance
                            waste[node] += injection[node] * wasting_factor
                            injection[node] -= injection[node] * wasting_factor
                        end
                    end
                        
                elseif subgraph_balance < 0-tolerance
                    
                    if emergency_mode == true
                        
                        subgraph_battery_content = 0
                        for node in c
                        subgraph_battery_content += (content_old[node] +
                                                     battery_use[node] + battery_emergency_use[node] )
                        end
                        subgraph_battery_content = round(subgraph_battery_content, digits=14)
                        
                        if subgraph_battery_content >= abs(subgraph_balance)
                        
                            for node in c
                                if subgraph_battery_content > 0
                                    charge = (content_old[node] + battery_use[node] + battery_emergency_use[node] )*
                                              minimum([1, (abs(subgraph_balance)/subgraph_battery_content) ])
                                else
                                    charge = 0
                                end
                                battery_emergency_use[node] -= charge
                                injection[node] += charge
                            end
                        end
                    else
                        for node in c
                            
                            if emergency_mode == true
                                
                                if injection[node] > 0+tolerance
                                    charge_capacity = (battery_map[node] - content_old[node] -
                                                       battery_use[node] - battery_emergency_use[node])
                                    charge = minimum([charge_capacity, injection[node]])
                                    battery_emergency_use[node] += charge
                                    injection[node] -= charge
                                    
                                elseif injection[node] < 0-tolerance
                                    charge_capacity = content_old[node] +
                                                      battery_use[node] + battery_emergency_use[node]
                                    charge = minimum([charge_capacity, abs(injection[node])])
                                    battery_emergency_use[node] -= charge
                                    injection[node] += charge
                                end
                            end
                            
                            for node in c
                                waste[node] += maximum([0, injection[node]])
                                injection[node] = 0
                            end
                        end
                    end
                end
            end
        end
        
        if ne(g) == 0
            break
        else
            diff = caps - abs.(flow(g,injection))
        end
    end
    
    battery_content = content_old + battery_use + battery_emergency_use
    
    demand[slack_index] = 0
    supply[slack_index] = 0
    battery_content[slack_index] = 0
    battery_use[slack_index] = 0
    battery_emergency_use[slack_index] = 0
    injection[slack_index] = 0
    waste[slack_index] = 0
    
    for node in 1:NV
        if (supply[node] > demand[node])
            waste[node] = minimum([waste[node],supply[node]-demand[node]])
        end
    end
    
    power = battery_struct(demand, supply, battery_content, battery_use, battery_emergency_use,
                           injection, waste)
    
    tolerance = 10^-14 #0
    
    lack = -power.supply +power.demand +power.battery_use +power.battery_emergency_use +power.injection +power.waste
    
    served_demand = power.demand - lack #power.demand + mismatch.*Int.(mismatch.<0)
    
    if sum(power.demand) != 0
        delta = 0 #sum(served_demand)/sum(power.demand)
    else
        delta = 1
    end
    
    if -tolerance <= delta <0
        delta = 0
    elseif 1 < delta <= 1+tolerance
        delta = 1
    end
    
    @assert delta >= 0 "negative delta: $(delta)"
    @assert delta <= 1 "delta above 1: $(delta)"
    
    damage = abs.(lack) + abs.(power.waste)
    
    #aka served transmission, served injection?
    served_mismatch = abs.(injection_0) - damage # damage ohne abs?
    
    if sum(abs.(injection_0)) != 0
        tau = sum(served_mismatch)/sum(abs.(injection_0))
    else
        tau = 1
    end
    
    if -tolerance <= tau <0
        tau = 0
    elseif 1 < tau <= 1+tolerance
        tau = 1
    end
    
    @assert tau >= 0 "negative tau: $(tau)"
    @assert tau <= 1 "tau above 1: $(tau)"
    
    return triple_struct(delta, tau, power.battery_content)
end

# a deprecated (and probably erroneous) version of the above function
function batt_sus_old(grid, capacities, demand, supply, battery_map, content_old, reliance, slack_index, emergency_mode)
    
    tolerance = 0 #10^-16 #?
    
    g = graphcopy(grid.graph)
    caps = copy(capacities)
    
    demand[slack_index] = 0
    supply[slack_index] = 0
    
    NV = length(demand)
    
    battery_content = zeros(NV) #at end of time step
    battery_use = zeros(NV)
    battery_emergency_use = zeros(NV)
    injection = zeros(NV)
    waste = zeros(NV)
     
    for node in 1:NV
        if battery_map[node] > 0+tolerance
            if supply[node]-demand[node] > 0+tolerance
                charge_capacity = (battery_map[node]-content_old[node])
                battery_use[node] = minimum([charge_capacity, reliance*abs(supply[node]-demand[node])])
            elseif supply[node]-demand[node] < 0-tolerance
                charge_capacity = content_old[node]
                battery_use[node] = - minimum([charge_capacity, reliance*abs(supply[node]-demand[node])])
            end
        end
    end
    
    injection = supply - demand - battery_use
    injection[slack_index] = -sum(injection)
    
    diff = caps - abs.(flow(g, injection))
    
    while minimum(diff) < 0-tolerance
        
        indices = findall(x -> x<0-tolerance, diff)
        sort!(indices,rev=true)
        for i in indices
            e = index2edge(g,i)
            rem_edge!(g, e)
            deleteat!(caps,i)
        end
        
        cc = connected_components(g)            
            
        injection[slack_index] = 0
        for c in cc
            subgraph_balance = 0
            for node in c
                subgraph_balance += injection[node]
            end
            subgraph_balance = round(subgraph_balance, digits=14)
            
            if slack_index in c
                injection[slack_index] = -subgraph_balance
                
            else
                if subgraph_balance > 0+tolerance
                    
                    subgraph_battery_vacancy = 0
                                        
                    if emergency_mode == true
                        for node in c
                        subgraph_battery_vacancy += (battery_map[node] - content_old[node] 
                                                     - battery_use[node] - battery_emergency_use[node] )
                        end
                        subgraph_battery_vacancy = round(subgraph_battery_vacancy, digits=14)
                        
                        for node in c
                            if subgraph_battery_vacancy > 0
                                charge = (battery_map[node] - content_old[node]
                                          - battery_use[node] - battery_emergency_use[node] )*
                                          minimum([1, subgraph_balance/subgraph_battery_vacancy])
                            else
                                charge = 0
                            end
                            battery_emergency_use[node] += charge
                            injection[node] -= charge
                        end
                    end
                    
                    subgraph_infeed = 0
                    for node in c
                        if injection[node] > 0+tolerance
                            subgraph_infeed += injection[node]
                        end
                    end
                    subgraph_infeed = round(subgraph_infeed, digits=14)
                    
                    if subgraph_infeed > 0
                        wasting_factor = minimum(
                                         [1, maximum([0, (subgraph_balance-subgraph_battery_vacancy)/subgraph_infeed ]) ])
                    else
                        wasting_factor = 0
                    end
                    
                    @assert wasting_factor <= 1+10^-14 "wasting factor > 1: w=$(wasting_factor), J=$(subgraph_balance), V=$(subgraph_battery_vacancy), J+=$(subgraph_infeed)"
                    
                    for node in c
                        if injection[node] > 0+tolerance
                            waste[node] += injection[node] * wasting_factor
                            injection[node] -= injection[node] * wasting_factor
                        end
                    end
                        
                elseif subgraph_balance < 0-tolerance
                    
                    if emergency_mode == true
                        
                        subgraph_battery_content = 0
                        for node in c
                        subgraph_battery_content += (content_old[node] +
                                                     battery_use[node] + battery_emergency_use[node] )
                        end
                        subgraph_battery_content = round(subgraph_battery_content, digits=14)
                        
                        if subgraph_battery_content >= abs(subgraph_balance)
                        
                            for node in c
                                if subgraph_battery_content > 0
                                    charge = (content_old[node] + battery_use[node] + battery_emergency_use[node] )*
                                              minimum([1, (abs(subgraph_balance)/subgraph_battery_content) ])
                                else
                                    charge = 0
                                end
                                battery_emergency_use[node] -= charge
                                injection[node] += charge
                            end
                        end
                    else
                        for node in c
                            
                            if emergency_mode == true
                                
                                if injection[node] > 0+tolerance
                                    charge_capacity = (battery_map[node] - content_old[node] -
                                                       battery_use[node] - battery_emergency_use[node])
                                    charge = minimum([charge_capacity, injection[node]])
                                    battery_emergency_use[node] += charge
                                    injection[node] -= charge
                                    
                                elseif injection[node] < 0-tolerance
                                    charge_capacity = content_old[node] +
                                                      battery_use[node] + battery_emergency_use[node]
                                    charge = minimum([charge_capacity, abs(injection[node])])
                                    battery_emergency_use[node] -= charge
                                    injection[node] += charge
                                end
                            end
                            
                            for node in c
                                waste[node] += maximum([0, injection[node]])
                                injection[node] = 0
                            end
                        end
                    end
                end
            end
        end
        
        if ne(g) == 0
            break
        else
            diff = caps - abs.(flow(g,injection))
        end
    end
    
    battery_content = content_old + battery_use + battery_emergency_use
    
    demand[slack_index] = 0
    supply[slack_index] = 0
    battery_content[slack_index] = 0
    battery_use[slack_index] = 0
    battery_emergency_use[slack_index] = 0
    injection[slack_index] = 0
    waste[slack_index] = 0
    
    for node in 1:NV
        if (supply[node] > demand[node])
            waste[node] = minimum([waste[node],supply[node]-demand[node]])
        end
    end
    
    power = battery_struct(demand, supply, battery_content, battery_use, battery_emergency_use,
                           injection, waste)
    
    tolerance = 10^-14 #0
    
    lack = -power.supply +power.demand +power.battery_use +power.battery_emergency_use +power.injection +power.waste
    
    served_demand = power.demand - lack #power.demand + mismatch.*Int.(mismatch.<0)
    
    if sum(power.demand) != 0
        delta = sum(served_demand)/sum(power.demand)
    else
        delta = 1
    end
    
    if -tolerance <= delta <0
        delta = 0
    elseif 1 < delta <= 1+tolerance
        delta = 1
    end
    
    @assert delta >= 0 "negative delta: $(delta)"
    @assert delta <= 1 "delta above 1: $(delta)"
    
    damage = lack + power.waste
    
    #aka served transmission, served injection?
    served_mismatch = abs.(power.supply-power.demand) - abs.(damage) # damage ohne abs?
    
    if sum(abs.(power.supply-power.demand)) != 0
        tau = sum(served_mismatch)/sum(abs.(power.supply-power.demand))
    else
        tau = 1
    end
    
    if -tolerance <= tau <0
        tau = 0
    elseif 1 < tau <= 1+tolerance
        tau = 1
    end
    
    @assert tau >= 0 "negative tau: $(tau)"
    @assert tau <= 1 "tau above 1: $(tau)"
    
    return triple_struct(delta, tau, power.battery_content)
end

# function that calculates the battery behavior and sustainant time series at the same time
function batt_sus_series(grid, caps, demand, supply, battery_map, initial_content,
                         reliance, slack_index, emergency_mode, limits)
    
    timesteps, NV = size(demand)
    
    delta = zeros(timesteps)
    tau = zeros(timesteps)
    
    content_0 = battery_map .* initial_content

    batt_0 = batt_sus(grid, caps, demand[1,:], supply[1,:], battery_map, content_0,
                       reliance, slack_index, emergency_mode, limits)
    content = batt_0.battery_content
    delta[1] = batt_0.delta
    tau[1] = batt_0.tau

    for t in 2:timesteps
        batt = batt_sus(grid, caps, demand[t,:], supply[t,:], battery_map, content,
                         reliance, slack_index, emergency_mode, limits)
        content = batt.battery_content
        delta[t] = batt.delta
        tau[t] = batt.tau
    end
    
    return sustainant_struct(delta,tau)     
end

# function to quickly extract a range of quantiles from a series of values
function quantile_list(series, p_values=[0, 10^-10, 10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2,
                    10^-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    
    #@assert all(series.==1.0) == false
    return [quantile(series,p) for p in p_values]
end

# function that calculates the inverse of a quantile (also called quantile function = qf)
function qf(series; threshold=1)
    return length(filter(x -> x<threshold, series))/length(series)
end

# function that performs a complete power grid simulation and resilience analysis
# for a single adverse influence and a single response strategy
function adaptation_result(P2, prosumer_ratio,
                           allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance;
                           grids = default_grids, locations = default_locs, 
                           line_safety = default_safety, flowtest_runs = default_flowtest_runs,
                           risktest_runs = default_risktest_runs, performance_days = default_days,
                           N = 100, P1 = 0, 
                           time_resolution = 1,
                           initial_content = 0.5, digits = 16,
                           n0 = maximum([Int(floor(N*0.15)),1]), p = 0.15 ,q = 0.15, r = 0.2,
                           w = x->x, emergency_mode=true)

    solar_dict = smoother(solar_dict_default,time_resolution)
    consumer_dict = smoother(consumer_dict_default,time_resolution)
    steps_per_day = Int(ceil(24*60/time_resolution))
    
    Delta_P = P2-P1
    
    results = sustainant_struct(zeros(grids*locations), zeros(grids*locations))
    
    limits = lower_upper_struct(-Inf, Inf)
    if battery_reliance > 0
        limits = limit_finder(prosumer_ratio, P2, battery_reliance;solar_dict=solar_dict, consumer_dict=consumer_dict)
    end
    
    for g in 1:grids

        grid = create_grid(N, n0, p, q, r)
        NV = N
        NE = ne(grid)
        slack_index = findmax(closeness_centrality(grid.graph))[2]
        
        caps = []
        line_budget_initial = 0
        
        for L in 1:locations
            
            solar_map = set_prosumers(NV,P1,P2,slack_index)
            if P2 > 0
                @assert sum(solar_map.new) > 0
            end
            
            if P1 == 0 && L > 1
                caps = caps
            else
                flowtest_p_timeseries = power_randomdraws(solar_map.old, solar_dict, consumer_dict,
                                                       prosumer_ratio, flowtest_runs, slack_index)
                caps = floor.(line_safety*max_flows(grid,flowtest_p_timeseries, slack_index), digits=digits)
                caps = zero_corrector(caps, collect(1:NE))
                @assert all(caps.>0) "$caps"
                line_budget_initial = sum(caps .* line_lengths(grid))
            end
            
            
            line_budget = line_budget_factor * line_budget_initial            
            
            new_grid = embedded_graphcopy(grid)
            new_caps = copy(caps)
            
            if line_budget == 0

                #no grid change

            elseif line_budget > 0
                
                if allocation == "uniform"
                    NE += new_lines
                    for L in 1:new_lines
                        nl = add_line(new_grid, new_caps, r)
                        new_grid = nl.grid
                        new_caps = nl.caps
                    end
                end
                
                @assert NE == ne(new_grid)
                
                new_line_locs = findall(e -> e == 0, new_caps)

                if allocation == "uniform"
                    
                    budget_weights = edge_weights_uniform(new_grid)
                    while sum(budget_weights) > 1
                        budget_weights /= sum(budget_weights)
                    end
                    @assert sum(budget_weights) <= 1
                    cap_boosts = budget_weights * line_budget ./ line_lengths(new_grid)
                    new_caps = new_caps + cap_boosts

                elseif allocation == "closeness"
                    
                    exponent = new_lines
                    budget_weights = edge_weights_by_closeness(new_grid, solar_map.new, exponent)
                    while sum(budget_weights) > 1
                        budget_weights /= sum(budget_weights)
                    end
                    @assert sum(budget_weights) <= 1
                    cap_boosts = budget_weights * line_budget ./ line_lengths(new_grid)
                    new_caps = new_caps + cap_boosts

                elseif allocation == "risktest"

                    risktest_p_timeseries = power_randomdraws(solar_map.new, solar_dict, consumer_dict,
                                                           prosumer_ratio, risktest_runs, slack_index)
                    budget_weights = edge_weights_by_violation(new_grid, new_caps, risktest_p_timeseries, new_line_locs)
                    while sum(budget_weights) > 1
                        budget_weights /= sum(budget_weights)
                    end
                    @assert sum(budget_weights) <= 1
                    cap_boosts = budget_weights * line_budget ./ line_lengths(new_grid)
                    new_caps = new_caps + cap_boosts
                    
                end
            end
            
            @assert all(isfinite.(new_caps))
            
            new_caps = floor.(new_caps, digits=digits)
            @assert all(new_caps.>0) "$new_caps"
            
            battery_budget = (N-1)*24*60*battery_budget_factor
            
            if battery_budget == 0

                battery_map = 0*solar_map.new

            elseif battery_budget > 0
                
                if allocation == "uniform"

                    battery_size = battery_budget / P2
                    battery_map = battery_size * solar_map.new

                elseif allocation == "closeness"
                    
                    weight_vector = node_weights_by_closeness(grid, solar_map.new)
                    while sum(weight_vector) > 1
                        weight_vector /= sum(weight_vector)
                    end
                    @assert sum(weight_vector) <= 1 "sum=$(sum(weight_vector))"
                    battery_map = battery_budget * weight_vector .* solar_map.new

                elseif allocation == "risktest"

                    risktest_p_timeseries = power_randomdraws(solar_map.new, solar_dict, consumer_dict,
                                                           prosumer_ratio, risktest_runs, slack_index)
                    weight_vector = node_weights_by_violation(grid, caps, risktest_p_timeseries)
                    while sum(weight_vector) > 1
                        weight_vector /= sum(weight_vector)
                    end
                    @assert sum(weight_vector) <= 1 "sum=$(sum(weight_vector))"
                    battery_map = battery_budget * weight_vector .* solar_map.new

                end
                
                battery_map = floor.(battery_map, digits=digits)
            end
            
            demand_supply_series = power_timeseries(solar_map.new, solar_dict, consumer_dict,
                                                   prosumer_ratio, performance_days, time_resolution, slack_index)
                
            sus = batt_sus_series(new_grid, new_caps, demand_supply_series.demand, demand_supply_series.supply,
                                  battery_map, initial_content, battery_reliance, slack_index, emergency_mode,
                                  limits)
            
            # failure fraction
            #results.delta[(g-1)*locations+L] = qf(sus.delta)
            #results.tau[(g-1)*locations+L] = qf(sus.tau)
            
            # mean, which is equivalent to integral
            #results.delta[(g-1)*locations+L] = mean(sus.delta)
            #results.tau[(g-1)*locations+L] = mean(sus.tau)
            
            # mean of weighted deviations from optimum
            results.delta[(g-1)*locations+L] = mean( w.( 1 .- sus.delta))
            results.tau[(g-1)*locations+L] = mean( w.( 1 .- sus.tau))            
        end            
    end  
    
    return results
    
end

# function that performs the above resilience analysis calculation and saves the results as a .csv file
function adaptation_result_file(P2, prosumer_ratio,
                                allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                                path, w, emergency_mode;
                                grids = default_grids, locations = default_locs, 
                                line_safety = default_safety, flowtest_runs = default_flowtest_runs,
                                risktest_runs = default_risktest_runs, performance_days = default_days,
                                N = 100, P1 = 0, 
                                time_resolution = 1,
                                initial_content = 0.5, digits = 16,
                                n0 = maximum([Int(floor(N*0.15)),1]), p = 0.15 ,q = 0.15, r = 0.2,
                                filename_override="")
    
    filename = ""
    
    if P2-P1 == 0
        prosumer_ratio = 0
    end
    
    if prosumer_ratio == 0
        P2 = copy(P1)
    end
    
    if battery_reliance == 0
        battery_budget_factor = 0
    end
    
    if battery_budget_factor == 0
        battery_reliance = 0
    end
    
    if line_budget_factor == 0
        new_lines = 0
    end
    
    if line_budget_factor + battery_budget_factor == 0
        allocation = "none"
    end
    
    if allocation == "none"
        line_budget_factor = 0
        battery_budget_factor = 0
        filename = "sus_for_P2=$(P2),_rp=$(round(prosumer_ratio,digits=4)),_alloc=$(allocation)"
    elseif line_budget_factor == 0
        filename = "sus_for_P2=$(P2),_rp=$(round(prosumer_ratio,digits=4))"*
                   ",_alloc=$(allocation),_em=$(emergency_mode)"*
                   ",_batt_bf=$(battery_budget_factor),_b_rel=$(battery_reliance)"
    elseif battery_budget_factor == 0
        filename = "sus_for_P2=$(P2),_rp=$(round(prosumer_ratio,digits=4))"*
                   ",_alloc=$(allocation)"*
                   ",_line_bf=$(line_budget_factor),_new_l=$(new_lines)"
    end
    filename = replace(filename,Pair(".","p"))*".csv"
    
    if filename_override !== ""
        filename = filename_override*".csv"
    end
    
    folder = readdir(path)
    if (filename in folder) == false
        
        parameter_names = ["N", "line_safety", "P1",
                           "N0", "p", "q", "r",
                           "P2", "Delta_P", "prosumer_ratio",
                           "allocation", "line_budget_factor", "new_lines", "battery_budget_factor", "battery_reliance",
                           "emergency_mode", "initial_content", "grids", "locations", "flowtest_runs", "risktest_runs",
                           "performance_days", "time_resolution", "digits"]
        Delta_P = P2-P1
        parameter_values = [N, line_safety, P1,
                            n0, p, q, r,
                            P2, Delta_P, prosumer_ratio,
                            allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                            emergency_mode, initial_content, grids, locations, flowtest_runs, risktest_runs,
                            performance_days, time_resolution, digits]
        emptycolumn = fill(" ",length(parameter_names))
        parameters = hcat(parameter_names, parameter_values, emptycolumn)
        emptyrow = [" " " " " "]
        
        header = ["instance" "delta mean weighted deviation" "tau mean weighted deviation"] #
        
        results = adaptation_result(P2, prosumer_ratio,
                                allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                                grids = grids, locations = locations, 
                                line_safety = line_safety, flowtest_runs = flowtest_runs, risktest_runs = risktest_runs,
                                performance_days = performance_days,
                                N = N, P1 = P1, 
                                time_resolution = time_resolution,
                                initial_content = initial_content, digits = digits,
                                n0 = n0, p = p, q = q, r = r, w=w, emergency_mode=emergency_mode)
        
        instances = collect(1:length(results.delta))
        list = hcat(instances, results.delta, results.tau)
        
        tab = vcat(parameters,emptyrow,header,list)

        CSV.write(path*"/"*filename,Tables.table(tab),delim="\t",header=false)
        
    else
        println("file $(filename) already exists in folder $(path) . skipping computation")
    end
end

# function that checks whether the frequency of deviations (no matter the magnitude) from the ideal sustainant value (=1) are rare enough
function stability_check(sus_series)
    
    #if quantile(sus_series,10^-5) == 1
    if qf(sus_series) <= 10^-5 #
        return true
    else
        return false
    end
end    

# function that returns the frequency of the sustainant deviating from the ideal value by more than a threshold
function resilience_freq(mean_weighted_deviations; threshold=10^-5)
    return length(filter(x -> x<=threshold, mean_weighted_deviations))/length(mean_weighted_deviations)
end

# function that determines the maximum prosumer ratio that a given response strategy can cope with
# the ratio is determined by interval halving up to a set resolution
function ratio_finder(allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                      ratio_res, max_ratio, w, emergency_mode)
    
    P2 = 1 #2 for closeness
    r_left = 0.
    r_right = copy(max_ratio)
    ratio = copy(r_right)
    while (r_right-r_left) > ratio_res
        
        results = adaptation_result(P2, ratio,
                               allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                               w=w, emergency_mode=emergency_mode)
        
        if resilience_freq(results.delta) == 0 && resilience_freq(results.tau) == 0 ###only one of the sustainants?
            r_right = ratio
        else
            r_left = ratio
        end
        ratio = 0.5*(r_left + r_right)
    end
    
    return r_right
end

# function that determines the maximum prosumer count that a given response strategy can cope with
# the count is determined by interval halving up to a set resolution
function P2_finder(allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                   P2_res, max_P2, w, emergency_mode)
    
    ratio = 0.5
    p_low = 1
    p_high = copy(max_P2)
    p = copy(p_high)
    while (p_high-p_low) > P2_res
        
        results = adaptation_result(p, ratio,
                               allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                               w=w, emergency_mode=emergency_mode)
        
        if resilience_freq(results.delta) == 0 && resilience_freq(results.tau) == 0 ###only one of the sustainants?
            p_high = p
        else
            p_low = p
        end
        p = Int(round(0.5*(p_low + p_high)))
    end
    
    return p_high
end

# function that extracts a numeric value from a filename (or any string) in the form of "...,_variable=value,_..."
function value_from_filename(filename, variable, splitter=",_")
    chunks = string.(split(filename,splitter))
    found = false
    for c in chunks
        if occursin(variable,c)
            value = replace(string.(split(c,"="))[2],Pair("p","."))
            found = true
            return Meta.parse(value)
        end
    end
    if found == false
        return 0
    end
end

# function that extracts a string value from a filename (or any string) in the form of "...,_variable=value,_..."
function string_value_from_filename(filename, variable, splitter=",_")
    chunks = string.(split(filename,splitter))
    found = false
    for c in chunks
        if occursin(variable,c)
            value = replace(string.(split(c,"="))[2],Pair("p","."))
            found = true
            return value
        end
    end
    if found == false
        return ""
    end
end

# function that creates subfolders with combinations of response strategies and their parameters
function qmc_folder_prep(allocation, line_budget_factor_range, new_lines_range,
                         battery_budget_factor_range, battery_reliance_range,
                         results_path, weight_function, emergency_mode;
                         ratio_res = 1.5, P2_res = 15, max_P2 = 99, max_ratio = 10.)
    
    allocation_folder = ""
    if allocation == "none"
        allocation_folder = "no_adaptation"
    elseif line_budget_factor_range == 0:0 && battery_budget_factor_range == 0:0
        allocation_folder = "no_adaptation"
    elseif line_budget_factor_range == 0:0 && battery_budget_factor_range != 0:0
        if emergency_mode == true
            allocation_folder = "emergency_battery_adaptation_"*allocation
        else
            allocation_folder = "simple_battery_adaptation_"*allocation
        end
    elseif line_budget_factor_range != 0:0 && battery_budget_factor_range == 0:0
        allocation_folder = "line_adaptation_"*allocation
    end
    
    combos = [(axis, LBF, NL, BBF, BR) for axis in ["rp","np"],
                                           LBF in line_budget_factor_range,
                                           NL in new_lines_range,
                                           BBF in battery_budget_factor_range,
                                           BR in battery_reliance_range]
                            
    function main_call(combo)
        (axis, LBF, NL, BBF, BR) = combo
        if axis == "rp"
            return ratio_finder(allocation, LBF, NL, BBF, BR, ratio_res, max_ratio, weight_function, emergency_mode)
        elseif axis == "np"
            return P2_finder(allocation, LBF, NL, BBF, BR, P2_res, max_P2, weight_function, emergency_mode)
        end
    end

    #cutoff = pmap(main_call, combos)

    for i in 1:2:length(combos)
        max_r, max_p = 10, 99 #cutoff[i], cutoff[i+1]
        (axis, LBF, NL, BBF, BR) = combos[i]

        parameter_folder = ""

        if allocation == "none"
            parameter_folder = "sus_samples_for_alloc=$(allocation)"*
                          ",_max_P2=$(max_p),_max_rp=$(max_r)"
        elseif LBF == 0
            parameter_folder = "sus_samples_for_alloc=$(allocation)"*",_em=$(emergency_mode)"*
                          ",_batt_bf=$(BBF),_b_rel=$(BR)"*
                          ",_max_P2=$(max_p),_max_rp=$(max_r)"
        elseif BBF == 0
            parameter_folder = "sus_samples_for_alloc=$(allocation)"*
                          ",_line_bf=$(LBF),_new_l=$(NL)"*
                          ",_max_P2=$(max_p),_max_rp=$(max_r)"
        end

        parameter_folder = replace(parameter_folder,Pair(".","p"))

        parameter_folder_path = results_path*"/"*allocation_folder*"/"*parameter_folder
        mkpath(parameter_folder_path)
    end
    
    return allocation_folder
end

# function that generates the path to a subfolder given the corresponding parameters
function qmc_folder_pointer(allocation, line_budget_factor_range, new_lines_range,
                         battery_budget_factor_range, battery_reliance_range,
                         results_path, weight_function, emergency_mode;
                         ratio_res = 1.5, P2_res = 15, max_P2 = 99, max_ratio = 10.)
    
    allocation_folder = ""
    if allocation == "none"
        allocation_folder = "no_adaptation"
    elseif line_budget_factor_range == 0:0 && battery_budget_factor_range == 0:0
        allocation_folder = "no_adaptation"
    elseif line_budget_factor_range == 0:0 && battery_budget_factor_range != 0:0
        if emergency_mode == true
            allocation_folder = "emergency_battery_adaptation_"*allocation
        else
            allocation_folder = "simple_battery_adaptation_"*allocation
        end
    elseif line_budget_factor_range != 0:0 && battery_budget_factor_range == 0:0
        allocation_folder = "line_adaptation_"*allocation
    end
        
    return allocation_folder
end

# function that performs the QMC sampling of the resilience analysis
# it is parallelized to optimze computation time
# the argument "offset" allows to expand on an existing QMC sample
function qmc_sampler(allocation_path, weight_function; monte_carlo_runs = mc_runs, offset = 0, skip_errors = true)
    
    allocation_directory = readdir(allocation_path)
    filter!(e->!occursin("csv",e), allocation_directory)
    
    for parameter_folder in allocation_directory
        
        max_P2 = value_from_filename(parameter_folder, "max_P2")
        max_ratio = value_from_filename(parameter_folder, "max_rp")
        allocation = string_value_from_filename(parameter_folder, "alloc")
        line_budget_factor = value_from_filename(parameter_folder, "line_bf")
        new_lines = value_from_filename(parameter_folder, "new_l")
        battery_budget_factor = value_from_filename(parameter_folder, "batt_bf")
        battery_reliance = value_from_filename(parameter_folder, "b_rel")
        emergency_mode = value_from_filename(parameter_folder, "em")
        if emergency_mode == 0; emergency_mode = false; end
        
        @assert max_ratio > 0.1
        @assert max_P2 > 1
        
        samples = Array{Any}(QuasiMonteCarlo.sample(monte_carlo_runs+offset,
                                                    [0.1, 0.5 + 10^-16], [max_ratio, max_P2 + 0.5 - 10^-16],
                                                    SobolSample()))[:,offset+1:end]
        samples[2,:] = Int.(round.(samples[2,:]))
        
        main_call(sample) = adaptation_result_file(sample[2], sample[1],
                allocation, line_budget_factor, new_lines, battery_budget_factor, battery_reliance,
                allocation_path * "/" * parameter_folder, weight_function, emergency_mode)
        
        if skip_errors == true
            pmap(main_call, [samples[:,i] for i in 1:monte_carlo_runs]; on_error = identity,
                 retry_delays = ExponentialBackOff(n = 3))
        else
            pmap(main_call, [samples[:,i] for i in 1:monte_carlo_runs]; on_error = throw)
        end
    end
end
