
using BenchmarkTools, Revise, CSV, DataFrames, DataFramesMeta, Parameters

# Set the working directory
if Sys.iswindows()
    path = raw"C:\Users\yangk\Dropbox\Research\whaling_commons\replication\navigating-the-commons"
elseif Sys.islinux()
    #path = raw"/u/project/rafey/ykyun/2022-Whaling-Regulation/"
    # path = raw"/mnt/phd/yyun/navigating-the-commons"
    path = dirname(@__DIR__)
end

# benchmark using 4 threads
println("The number of threads is $(Threads.nthreads())")

# Load function call
string(@__DIR__) in LOAD_PATH || push!(LOAD_PATH, @__DIR__);
using WhalesModule;
const M = WhalesModule;

m = M.DG_model(;  β=[],                                       # production function paratmeters
                F="CD",                                     # functional form of production function. e.g. `CD` is Cobb-Douglas
                α=[],                                       # demand curve parameters
                ρ=.9,                                       # time discounting factor
                xᵉ=1, δ=0.0,                                # ignore for now..
                ψ=1.0, Γ=[0.0;0.0], ϕ=0.0, κ=0.0, ϕᶠ=0.0,   # initial parameters. But need to be fixed later
                K_size=300,                                 # `K_size` refers to the incremental increase of capacity (fleet size). e.g. 300 means that the capacity space is [100, 400, 700, 1000, .....]. Maybe `K_width` is more appropriate naming. Let's fix later.
                A_size=15,                                  # `A_size` refers to the size of firm age space. It might be used later, but ignore for now.
                Ω_size=15,                                  # `Ω_size` refers to the size of productivity space. e.g. The total number of productivity grids is 15.
                x_size=0,                                   # `x_size` refers to the size of state space, as a combination of all state variables. The actual number will be allocated in the later part of the code.
                obj_incumbents_pre=0,obj_incumbents_post=0, # ignore for now..
                scale_Q=100, scale_P=1000,                  # scaling to avoid explosion in exponentials...
                N̄=400,                                      # maxmimum number of incumbents that the industry can accomodate
                nboot=100,                                  # ignore for now...
                state_space="KΩ")                           # state space is consisted of capacity and productivity

d = M.DG_data(; t0=1804, t1=1920,                             # initial year and last year
                game_start=1858, game_end=1920,             # need to be fixed... for now, they are anyway defined in the later part of the code
                sMat=[], PVec=[], WVec=[], KVec=[], QVec=[], NVec=[], S=[], Q=[], πMat=[], VMat=[],
                ιMat0=[],χMat0=[],λVec0=[],ιMat1=[],χMat1=[],λVec1=[],Π=[])

Tbar = d.t1 - d.t0 + 1 # last time period of non-stationary transition

d.Π = []

#========================================================================================#
#---------------------------------------- Setup -----------------------------------------#
#========================================================================================#

df = CSV.read(path*raw"/data/data_with_omega.csv", DataFrame)
rename!(df, :omega=>:ω)
dropmissing!(df, :ω)

# Get static parameters
PF_estimates = CSV.read(path*raw"/data/PF_CD.csv", DataFrame)
production_parameters = M.get_production_parameters(m.state_space, PF_estimates)

dm_estimates_pre = CSV.read(path*raw"/data/demand_curve_pre.csv", DataFrame)
dm_estimates_post = CSV.read(path*raw"/data/demand_curve_post.csv", DataFrame)
demand_parameters_pre = M.get_demand_parameters(dm_estimates_pre)
demand_parameters_post = M.get_demand_parameters(dm_estimates_post)
demand_parameters_post.αᵖ = demand_parameters_pre.αᵖ

# Construct state space
x, Ω_df, K_df, A_df, ω_trans = M.get_state_space!(m, df)

if m.state_space == "KAΩ"                
    println("The dimension of state space is $(m.K_size*m.A_size*m.Ω_size) with $(m.K_size) grid for K, $(m.A_size) grid for A, and $(m.Ω_size) grid for Ω")
elseif m.state_space == "KΩ"
    println("The dimension of state space is $(m.K_size*m.Ω_size) with $(m.K_size) grid for K and $(m.Ω_size) grid for Ω")
end

# Get aggregate variables
agg_output = CSV.read(path*raw"/data/agg_output.csv", DataFrame) 
agg_capacity = CSV.read(path*raw"/data/agg_capacity.csv", DataFrame) 
agg_incumbents = CSV.read(path*raw"/data/agg_incumbents.csv", DataFrame)
pop_df = CSV.read(path*raw"/data/whalepop_estimated.csv", DataFrame) 
usgdp_df = CSV.read(path*raw"/data/aggregate_gdp_pop_index.csv", DataFrame) 
whale_prices = CSV.read(path*raw"/data/demand_data.csv", DataFrame)[:,[:year,:price_whale_catchval_cpi]] 
dropmissing!(whale_prices)
rename!(whale_prices, :price_whale_catchval_cpi => :p)
petroleum = CSV.read(path*raw"/data/petroleum_price.csv", DataFrame)[:,[:year,:petroleum_price_ma]] 
rename!(petroleum, :petroleum_price_ma => :petroleum_price)
petroleum = [DataFrame(year=collect(d.t0:d.tP-1), petroleum_price=ones(d.tP-d.t0)); petroleum]

agg_variables = M.get_aggregate_variables(agg_output, agg_capacity, agg_incumbents, pop_df, usgdp_df, whale_prices, petroleum, d, Tbar)


# Get initial guess of expected industry state in each period directly from data
sMat_data = M.get_initial_expected_industry_state!(df, m.state_space, x, Ω_df, K_df, A_df, d)

# Get entrants' intial state
xᵉ, Π_ent, Kᵉˣ, Ωᵉˣ = M.get_entrants_state(df, x, Ω_df, K_df, m.state_space)



# For productivity transition
Π_base = reshape(ω_trans',(m.Ω_size*m.Ω_size))
Is = repeat(collect(1:m.x_size),inner=m.Ω_size)
if m.state_space == "KAΩ"                
    Xs = repeat(Π_base,outer=m.K_size*m.A_size)
elseif m.state_space == "KΩ"
    Xs = repeat(Π_base,outer=m.K_size)
end

#Is_Π, Js_Π = prepare_state_transition(m)
Is_Π, Js_Π, Xs_Π = M.prepare_state_transition(m, ω_trans)

m.options.ζ₁ = 2/3; m.options.Ζ₁ = 1.0
m.options.ζ₂ = 2/3; m.options.Ζ₂ = 1.0
m.options.ϵ₁ = 1e-10; m.options.ϵ₂ = 1e-10

Nᵖᵉ = 60
n_max = 200

m.cost_struct = 16

# ψ, κ, ϕᶠ, Γ for cost_struct == 16:
ψ = 1.0;
κ = 0.0; ϕ = 0.0;   
ϕᶠ = 0.0;  
Γ = [0.0;0.0;0.0;0.0;0.0];

m.κ = κ; m.ϕ = ϕ; m.Γ = Γ; m.ϕᶠ = ϕᶠ
m.ψ = 1.0 # set the logit scale parameter
println("Cost structure is $(m.cost_struct)")

m.scale_est = "Yes"
println("Estimate logit scale parameter? $(m.scale_est)")

#========================================================================================#
#-------------------------------------- Estimation --------------------------------------#
#========================================================================================#


# Estimate: golden-age
m.post = 0
d.game_start = 1831 
d.game_end = 1858
timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit = M.build_data(df, d, m, K_df, A_df, Ω_df, sMat_data)

K_space = K_df.K

#========================================================================================#

# @time M.solve_NOE(d, m)

@unpack QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z = agg_variables
K = x.K
Ω = x.Ω
if "A" ∈ names(x)
    A = x.A
end

if m.scale_est == "No"
    m.κ = θ[1]; m.ϕᶠ = θ[2]; m.Γ = θ[3:end]
else
    m.ψ = θ[1]; m.κ = θ[2]; m.ϕᶠ = θ[3]; m.Γ = θ[4:end]
end


θ₀ = reduce(vcat,[κ,ϕᶠ,Γ])
@time M.obj_func(θ₀)

pll = @time M.calc_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, K, production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space, Is, Xs, Π_ent, Xs_Π, Is_Π, Js_Π)


# @time M.calc_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, K, production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space, Is, Xs, Π_ent, Xs_Π, Is_Π, Js_Π)





@time nfxp_golden_age = M.optimize_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, x, production_parameters, demand_parameters_pre, demand_parameters_post, df, Nᵖᵉ, n_max, K_space, Is, Xs, Π_ent, Xs_Π, Is_Π, Js_Π)





# #========================================================================================#
# # within while loop
# #========================================================================================#
# λVec0 = ones(eltype(m.κ),Tbar + 1,m.K_size); λVec0[end] = 0            
# #d.ιMat0 = ones(Tbar+1,m.x_size,5)
# ιMat0 = ones(eltype(m.κ),Tbar+1,m.x_size,m.K_size)
# χMat0 = ones(eltype(m.κ),Tbar+1,m.x_size)
# λVec1 = ones(eltype(m.κ),Tbar+1,m.K_size); λVec1[end] = 0
# #d.ιMat1 = ones(Tbar+1,m.x_size,5)
# ιMat1 = ones(eltype(m.κ),Tbar+1,m.x_size,m.K_size)
# χMat1 = ones(eltype(m.κ),Tbar+1,m.x_size)

# # Initial matrices (period-by-state) and vectors (by period)
# sMat = zeros(eltype(m.κ),Tbar+1,m.x_size)   # industry state, i.e. the number of firms in each possible state
# PVec = zeros(eltype(m.κ),Tbar+1)            # price vector
# WVec = zeros(eltype(m.κ),Tbar+1)            # whale population in each period
# KVec = zeros(eltype(m.κ),Tbar+1)            # aggregate capacity in each period
# QVec = zeros(eltype(m.κ),Tbar+1)            # aggregate whale catch in each period
# QᶠVec = zeros(eltype(m.κ),Tbar+1)           # aggregate whale catch outside of US whaling
# NVec = zeros(eltype(m.κ),Tbar+1)            # aggregate number of incumbents in each period
# VMat = zeros(eltype(m.κ),Tbar+1,m.x_size)
# πMat = zeros(eltype(m.κ),Tbar+1,m.x_size)
# Q = zeros(eltype(m.κ),m.x_size)
# mc = zeros(eltype(m.κ),m.x_size)

# WVec[[1]] = W₀
# sMat[1,:] = sMat_data[1,:]
# NVec[[1]] .= max(sum(sMat[1,:]),1)
# KVec[[1]] .= sum(K .* sMat[1,:])
# calc_production!(Q, KVec, WVec, d, 1)
# QVec[[1]] .= sum(Q.*sMat[1,:])
# sMat[1:end-1,:] = sMat_data
# compute_aggregate_state!(NVec, KVec, QVec, QᶠVec, WVec, PVec, sMat, Q, d, m)

# n = 0
# Δ = 10; Δ₁ = 10; Δ₂ = 10; Δ₃ = 10; 


# sMat[1:end-1,:] = sMat_data
# @time compute_aggregate_state!(NVec, KVec, QVec, QᶠVec, WVec, PVec, sMat, Q, d, m)

# @time compute_backward!(ιMat1, χMat1, λVec1, πMat, VMat, Q, mc, KVec, WVec, PVec, d, m)
