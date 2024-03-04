
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
                state_space="KΩ"                           # state space is consisted of capacity and productivity
                )

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



# This is the main function to estimate the model. This takes a long time to run.
# @time nfxp_golden_age = M.optimize_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, x, production_parameters, demand_parameters_pre, demand_parameters_post, df, Nᵖᵉ, n_max, K_space, Is, Xs, Π_ent, Xs_Π, Is_Π, Js_Π)


#========================================================================================#
# Inside optimize_log_likelihood:

@unpack QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z = agg_variables
K = x.K
Ω = x.Ω
if "A" ∈ names(x)
    A = x.A
end


if m.scale_est == "No"
    θ₀ = reduce(vcat,[κ,ϕᶠ,Γ]) # Initial guess of θ=θ₀
    optimum = optimize(obj_func, 
                        θ₀, 
                        NelderMead(),
                        Optim.Options(outer_iterations = 1500,iterations=10000)) #Fminbox(NelderMead()) NelderMead() SimulatedAnnealing()
else # estimate logit scale parameter
    θ₀ = reduce(vcat,[m.ψ,m.κ,m.ϕᶠ,m.Γ]) # Initial guess of θ=θ₀
    if m.cost_struct == 1
        lower = [0.2,-Inf, -Inf, -Inf, 0.0,  -Inf, 0.0]
        upper = [5.0, Inf,  Inf, Inf, 10.0,  Inf, 10.0]
    elseif m.cost_struct == 4
        lower = [0.2,-Inf, -Inf, -Inf, -Inf]
        upper = [5.0, Inf,  Inf,  Inf,  Inf]
    elseif m.cost_struct == 11
        lower = [0.2,-Inf, -Inf, -Inf, -Inf,  0.0, -Inf, -Inf,  0.0]
        upper = [5.0, Inf,  Inf,  Inf,  Inf, 10.0,  Inf,  Inf, 10.0]
    elseif m.cost_struct == 14 
        lower = [0.2,-Inf, -Inf, 0,  -Inf, -Inf]
        upper = [5.0, Inf,  Inf, 10.0,  Inf,  Inf]
    elseif m.cost_struct == 15
        lower = [0.2,-Inf, -Inf, -Inf, -Inf, -Inf]
        upper = [5.0, Inf,  Inf,  Inf,  Inf,  Inf]
    elseif m.cost_struct == 16
        lower = [0.2,-Inf, -Inf,  0.0, -Inf, -Inf, -Inf, -Inf]
        upper = [5.0, Inf,  Inf, 10.0,  Inf,  Inf,  Inf,  Inf]
    end
end

#========================================================================================#
# Inside obj_func:

θ = θ₀ # First, set the initial guess of θ=θ₀

if m.scale_est == "No"
    m.κ = θ[1]; m.ϕᶠ = θ[2]; m.Γ = θ[3:end]
else
    m.ψ = θ[1]; m.κ = θ[2]; m.ϕᶠ = θ[3]; m.Γ = θ[4:end]
end


# Main function to benchmark. always benchmark with 4 threads.
M.calc_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, K, production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space, Is, Xs, Π_ent, Xs_Π, Is_Π, Js_Π) #Precompilation

@time M.calc_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, K, production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space, Is, Xs, Π_ent, Xs_Π, Is_Π, Js_Π)


#========================================================================================#
# Inside calc_log_likelihood:

@unpack QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z = agg_variables

@time exit_prob, invest_prob, entry_prob = M.solve_NOE(d, m, Tbar, agg_variables, sMat_data, K,production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space, Is, Xs, Π_ent, Xs_Π, Is_Π, Js_Π)[1:3]

j = 2
using SparseArrays
Π_Ω = sparse(Is, m.Js_vec[j], Xs, m.x_size, m.x_size)
@time Π_cond = M.calc_transition_productivity(m, Is, m.Js_vec[j], Xs,j, Π_Ω)


#=
feb 25:  15.008064 seconds (67.83 M allocations: 35.576 GiB, 28.84% gc time) -> 10.388061 seconds (54.77 M allocations: 32.615 GiB, 16.48% gc time)

=#
