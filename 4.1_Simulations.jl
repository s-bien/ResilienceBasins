# this script is used to compute the power grid performance with a fixed response strategy, but a range of parameters.
# it is optimized to run on a SLURM-based high performance cluster

# checkpoint print commands help identify errors
println("checkpoint 1")

# load packages to manage cluster
using Distributed, SlurmClusterManager
println("checkpoint 2")

# automatically add processes on all assigned CPUs
addprocs(SlurmManager())     # comment out if not on cluster
println("checkpoint 3")

# get a response from all CPUs
@everywhere println("Hello from worker $(myid()) @ $(gethostname()) !")
println("checkpoint 4")

# load package manager on all CPUs
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()
# @everywhere Pkg.resolve()
println("checkpoint 5")

# load data and function script on all CPUs
@everywhere begin
    consumption_path = "Data/consumption"
    solar_path = "Data/solar"
    functions_path = "1.0_Functions.jl"
    output_path = "output"
end
println("checkpoint 6")

# load packages on all CPUs
@everywhere begin
    using SyntheticNetworks
    using LightGraphs
    using StatsBase
    using EmbeddedGraphs
    using DelimitedFiles
    using CSV
    using IterativeSolvers
    using Tables
    using Distances
    using QuasiMonteCarlo
    #using ProgressMeter
    #using Colors
end
println("checkpoint 7")

# read and process data on all CPUs
@everywhere include(functions_path);
begin
    folder = readdir(consumption_path)
    consumption_days = length(folder)
    consumption_data = []

    for file in folder
        data = open(readdlm,consumption_path*"/"*file)
        stringy = string.(data[:,2])
        floaty = parse.(Float64,stringy)
        append!(consumption_data,floaty)
    end

    consumption_data /= mean(consumption_data);

    consumer_dict_raw = data2dict(consumption_data,consumption_days)

    consumer_dict_default = smoother(consumer_dict_raw,30);

    folder = readdir(solar_path)
    solar_days = length(folder)
    solar_data = []

    for file in folder
        data = open(readdlm,solar_path*"/"*file)
        stringy = string.(data[:,2])
        floaty = parse.(Float64,stringy)
        append!(solar_data,floaty)
    end

    solar_data /= mean(solar_data);

    solar_dict_default = data2dict(solar_data,solar_days);
end
println("checkpoint 8")

# set data as global variables on all CPUs
@allocated @everywhere begin
    solar_dict_default = $solar_dict_default
    consumer_dict_default = $consumer_dict_default
    consumption_days = $consumption_days
    solar_days = $solar_days
end
println("checkpoint 9")

# select the response strategy an the ranges for the associated parameters
# as well number of QMC samples & iterations
# setting TEST_MODE=true reduces the number of computations to quickly check for errors
@everywhere TEST_MODE = true
@everywhere skip_errors = false

@everywhere include(functions_path);
@everywhere begin
    
    allocation = "uniform" #fixed
    failure_weight_function = x -> x #fixed
    
    typ = "battery" # none / line / battery
    battery_emergency_mode = false
    
    if TEST_MODE == true
        mc_runs = 2
        
        default_grids = 1
        default_days = 1
        default_locs = 1

        default_safety = 1.75
        default_flowtest_runs = 1000
        default_risktest_runs = 0

        line_budget_factor_range = [0]
        new_lines_range = [0]

        battery_budget_factor_range = [1000, 100000]
        battery_reliance_range = [0.01, 0.1]
    else
        mc_runs = 256
        
        default_grids = 50
        default_days = 60
        default_locs = 1

        default_safety = 1.75
        default_flowtest_runs = 1000
        default_risktest_runs = 0

        line_budget_factor_range = [0.01, 1, 100]
        new_lines_range = [0, 10, 100]

        battery_budget_factor_range = [1000]
        battery_reliance_range = [0.1]
    end
    
    if allocation == "none" || typ == "none"
        line_budget_factor_range = [0]
        new_lines_range = [0]
        battery_budget_factor_range = [0]
        battery_reliance_range = [0]
    elseif typ == "battery"
        line_budget_factor_range = [0]
        new_lines_range = [0]
    elseif typ == "line"
        battery_budget_factor_range = [0]
        battery_reliance_range = [0]
    end
end

# set the output directory globally
dir = qmc_folder_prep(allocation, line_budget_factor_range, new_lines_range,
                      battery_budget_factor_range, battery_reliance_range,
                      output_path, failure_weight_function, battery_emergency_mode)

@everywhere dir = $dir
@everywhere println("prepared folders in $dir")

# main simulation
qmc_sampler(output_path*"/"*dir, failure_weight_function;
            monte_carlo_runs = mc_runs, offset = 0, skip_errors=skip_errors);
