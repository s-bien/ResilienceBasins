# this script is used to compute the power grid performance without any adaptation (the "status quo")
# it is optimized to run on a SLURM-based high performance cluster but can be tested locally on any CPU

# load packages to manage cluster
using Distributed
using ClusterManagers

# is the script being run on the cluster?
const SLURM = false

N_cpus = try
    parse(Int, ARGS[1])
catch
    3
end

# 1 main process and $N_worker sub processes
# on the cluster add all processes
SLURM ? addprocs(SlurmManager(N_cpus)) : addprocs(N_cpus - 1);

# load package manager on all CPUs
@everywhere using Pkg
@everywhere Pkg.activate(".")
# @everywhere Pkg.instantiate()
# @everywhere Pkg.resolve()

# load data and function script on all CPUs
@everywhere begin
    consumption_path = "Data/consumption"
    solar_path = "Data/solar"
    functions_path = "1.0_Functions.jl"
    output_path = "output"
end

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
    using ProgressMeter
    using Colors
end

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

# set data as global variables on all CPUs
@everywhere begin
    solar_dict_default = $solar_dict_default
    consumer_dict_default = $consumer_dict_default
    consumption_days = $consumption_days
    solar_days = $solar_days
end

# set the simulation parameter ranges and number of iterations
# setting TEST_MODE=true reduces the number of computations to quickly check for errors
@everywhere TEST_MODE = true
@everywhere skip_errors = false

@everywhere include(functions_path);
@everywhere begin
    
    allocation = "uniform" #fixed
    failure_weight_function = x -> x #fixed
    default_risktest_runs = 0 #fixed
    
    typ = "none" # none / line / battery
    battery_emergency_mode = false
    
    if TEST_MODE == true
        safety_range = [2.]
        flowtest_runs_range = [1000]
        reps = 1
        
        default_grids = 10
        default_days = 6
        default_locs = 1
    else
        safety_range = [1.5, 1.75, 2.]
        flowtest_runs_range = [500,750,1000]
        reps = 1
        
        default_grids = 50
        default_days = 60
        default_locs = 1
    end
end

# main simulation
@everywhere samples = [(s, ftr, rep) for s in safety_range, ftr in flowtest_runs_range, rep in 1:reps]

@everywhere main_call(sample) = ( mkpath(output_path * "/status_quo_test/safety=$(sample[1]),_flowtest_runs=$(sample[2])") ;
                                  adaptation_result_file(0, 0, typ, 0, 0, 0, 0,
                                  output_path * "/status_quo_test/safety=$(sample[1]),_flowtest_runs=$(sample[2])",
                                  failure_weight_function, battery_emergency_mode;
                                  line_safety = sample[1], flowtest_runs = sample[2],
                                  filename_override="rep $(sample[3])" ) )

pmap(main_call, samples; on_error = throw);


