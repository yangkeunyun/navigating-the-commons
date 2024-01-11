#=
# this script installs all required packages to run this script
using Pkg
println("Installation of required packages to run this script")
Pkg.update()
package_list = ["CSV", "DataFrames", "DelimitedFiles", "FixedEffectModels", "LinearAlgebra", 
                "Plots", "Optim", "Random", "Parameters", "StatsBase", "SplitApplyCombine", "GLM",
                "Impute", "StatsPlots", "Binscatters", "Combinatorics", "SparseArrays", "TickTock", "Kronecker",
                "Distributions", "CategoricalArrays", "Distances", "BenchmarkTools", "DataStructures", "Distributed", 
                "SpecialFunctions", "ForwardDiff", "Measures"]
for package in package_list
    if ! in(package,keys(Pkg.installed())) Pkg.add(package) end
end
@warn("Packages installed!")
=#
# Additional packages
using Pkg

using CSV, DataFrames, DelimitedFiles   # loading data
using FixedEffectModels # create lagged variables in DataFrame
using LinearAlgebra     # basic math
using Plots             # draw plots
using Optim # minimization of functions
using Random            # generate random numbers
using Parameters        # decorates a type definition to allow default values and a keyword constructor
using StatsBase
using SplitApplyCombine # split string variables into multiple columns
using GLM
using Impute
using StatsPlots
using Binscatters
using Combinatorics
using SparseArrays
using TickTock
using Kronecker
using Distributions
using CategoricalArrays # find quantiles
using Distances
using BenchmarkTools
using DataStructures
using Distributed
using SpecialFunctions
using Base.Threads
using Measures
using ColorSchemes
using ForwardDiff
#using Optimization, OptimizationOptimJL, ForwardDiff # In Linux server, I had trouble with installing 'OptimizationOptimJL'. So I did using Pkg
                                                     #                                                                                Pkg.activate("my_project")
                                                     #                                                                                Pkg.add("OptimizationOptimJL")
                                                     #                                                                       and then using Pkg
                                                     #                                                                                Pkg.update("OptimizationOptimJL")

                                                     
# Data structure
abstract type Data end                      # General type for any data
@with_kw mutable struct PF_data <: Data
    y::Array{Float64}                       # output quantity
    L1y::Array{Float64}
    k::Array{Float64}                       # predetermined input, e.g. capital
    L1k::Array{Float64}
    l::Array{Float64}                       # flexible inputs, e.g. labor, materials. fX[1] is a static variable optimized by firm every period, on the other hand, fX[2] is a proxy variable to control for unobserved productivity
    L1l::Array{Float64}                     # lagged flexible inputs, e.g. labor, materials 
    k′::Array{Float64}                      # proxy input
    a::Array{Float64}                       # firm age
    L1a::Array{Float64}                       # firm age
    z::Array{Float64}                       # others, e.g. aggregate variables: whale population, total number of voyages, total number of agents
    L1z::Array{Float64}
    fid::Array{Int64}                     # firm id
    tid::Array{Int64}                     # year
    Φ̂::Array{Float64}
    L1Φ̂::Array{Float64}
    ω̂::Array{Float64}
    L1ω̂::Array{Float64}
    χ::Array{Int64}
    P::Array{Float64}
    IV::Array{Float64} 
end

@with_kw mutable struct DG_data <: Data
    t0::Int64                      # after the war, 18 June 1812 – 17 February 1815. That is, pre-1816 is the first stationary equilibrium
    t1::Int64                      # long enought time span. e.g. post-1925 is the second stationary equilibrium
    tP=1859                        # year of petroleum discovery in Pennsylvania 

    game_start::Int64
    game_end::Int64

    sMat::Array{Float64}
    PVec::Vector{Float64}
    WVec::Vector{Float64}
    KVec::Vector{Float64}
    QVec::Vector{Float64}
    NVec::Vector{Float64}
    S::Array{Float64}
    Q::Vector{Float64}
    πMat::Array{Float64}
    VMat::Array{Float64}
    ιMat0::Array{Float64}
    χMat0::Array{Float64}
    λVec0::Array{Float64}
    ιMat1::Array{Float64}
    χMat1::Array{Float64}
    λVec1::Array{Float64}
    Π::Array{Float64}
end


abstract type Model end                     # General type for any model
@with_kw mutable struct PF_model <: Model
    β::Array{Float64}                       # parameters of inputs
    σ::Float64
    ν::Float64
    F="CD"
end

@with_kw mutable struct DG_model <: Model
    # Production parameters
    β::Array{Float64}                       # parameters of production inputs
    F="CD"

    # Demand parameters
    α::Array{Float64}

    # Constants
    ρ=.9                          # discount factor
    xᵉ::Int64                     # entry state: Ω = Ωᵉ, A = 1
    δ::Float64                    # capacity depreciation rate

    # Parameters to search
    ψ
    Γ
    ϕ
    κ
    ϕᶠ

    scale_est = "Yes"

    # States
    state_space="KAΩ"
    K_size::Int64
    Ω_size::Int64
    A_size::Int64
    x_size::Int64

    ## Forward simulation
    #n_sim_paths=250
    #n_sim_paths_perturb = 500

    cost_struct=4

    post=1 # post-petroleum or not

    n_incumbents_pre=0
    n_incumbents_post=0

    obj_incumbents_pre::Float64
    obj_incumbents_post::Float64

    scale_Q::Int64
    scale_P::Int64

    N̄::Int64

    nboot::Int64

    options = DG_options()		# Additional options 
end     

# Options structure gathering additional model/estimation options
abstract type Option end                    # General type for any option
@with_kw mutable struct PF_options <: Option
    K::Int                                  # number of bootstrap repititions
    S::Int                                  # size of each bootstrap draws
    G::Int                                  # degree of polynomial approximation
    R=false                                 # bootstrap resampling replacement? true or false
    seed::Int = 1234
    algo=NelderMead()                       # optimization algorithm (default NelderMead() from Optim.jl)
    tol::Float64 = 1e-20                    # tolerance level for GMM objective function value  
    iter::Int = 1000                        # number of iterations for minimization 
    VA=true
end

## Parameters for program controls;precisions and tolerances
@with_kw mutable struct DG_options
    ζ₁=1e-3             # exponential parameter on ι and χ
    Ζ₁=1e-3             # addition parameter on ι and χ
    ζ₂=1e-6             # exponential parameter on λ
    Ζ₂=1e-6             # addition parameter on λ
    ϵ₀=1e-5             # if λ is smaller than this we treat it as zero, i.e. zero entry
    ϵ₁=5                # tolerance parameter for zero profit condition
    ϵ₂=1                # tolerance parameter for χ and ι
end  


