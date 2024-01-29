module Whale

    using CSV, DataFrames, DelimitedFiles, FixedEffectModels, LinearAlgebra, Plots, Optim, Random, Parameters, StatsBase, SplitApplyCombine, GLM, Impute, StatsPlots, Binscatters, Combinatorics, SparseArrays, TickTock, Kronecker, Distributions, CategoricalArrays, Distances, BenchmarkTools, DataStructures, Distributed, SpecialFunctions, Base.Threads, Measures, ColorSchemes, ForwardDiff, Statistics, Parameters

    include("0_mutable_struct.jl")
    include("1_build_up.jl")
    include("2_model_estimation.jl")
    include("3_model_simulation.jl")
    include("get_aggregate_variables.jl")
    include("get_demand_parameters.jl")
    include("get_entrants_state.jl")
    include("get_initial_expected_industry_state.jl")
    include("get_production_parameters.jl")
    include("get_state_space.jl")
    
end