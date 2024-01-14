# Set the working directory
if Sys.iswindows()
    path = raw"C:\Users\yangk\Dropbox\Research\whaling_commons\replication\navigating-the-commons"
elseif Sys.islinux()
    #path = raw"/u/project/rafey/ykyun/2022-Whaling-Regulation/"
    path = raw"/mnt/phd/yyun/navigating-the-commons"
end

println("The number of threads is $(Threads.nthreads())")

# Load function call
println("Load functions and packages")
include(path*raw"/code/0_mutable_struct.jl")
include(path*raw"/code/1_build_up.jl")
include(path*raw"/code/2_model_estimation.jl")
include(path*raw"/code/3_model_simulation.jl")

println("Set data and model parameters")

m = DG_model(;  β=[],                                       # production function paratmeters
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

d = DG_data(; t0=1804, t1=1920,                             # initial year and last year
                game_start=1858, game_end=1920,             # need to be fixed... for now, they are anyway defined in the later part of the code
                sMat=[], PVec=[], WVec=[], KVec=[], QVec=[], NVec=[], S=[], Q=[], πMat=[], VMat=[],
                ιMat0=[],χMat0=[],λVec0=[],ιMat1=[],χMat1=[],λVec1=[],Π=[])

Tbar = d.t1 - d.t0 + 1 # last time period of non-stationary transition

d.Π = []

#========================================================================================#
#---------------------------------------- Setup -----------------------------------------#
#========================================================================================#

df = CSV.read(path*raw"/data/data_with_omega.csv", DataFrame) 

# Get static parameters
include(path*raw"/code/get_production_parameters.jl")
include(path*raw"/code/get_demand_parameters.jl")

# Construct state space
include(path*raw"/code/get_state_space.jl")

if m.state_space == "KAΩ"                
    println("The dimension of state space is $(m.K_size*m.A_size*m.Ω_size) with $(m.K_size) grid for K, $(m.A_size) grid for A, and $(m.Ω_size) grid for Ω")
elseif m.state_space == "KΩ"
    println("The dimension of state space is $(m.K_size*m.Ω_size) with $(m.K_size) grid for K and $(m.Ω_size) grid for Ω")
end

# Get aggregate variables
include(path*raw"/code/get_aggregate_variables.jl")

# Get initial guess of expected industry state in each period directly from data
include(path*raw"/code/get_initial_expected_industry_state.jl")

# Get entrants' intial state
include(path*raw"/code/get_entrants_state.jl")

# Get initial beliefs
#include(path*raw"/code/get_initial_beliefs.jl")

# For productivity transition
Π_base = reshape(ω_trans',(m.Ω_size*m.Ω_size))
Is = repeat(collect(1:m.x_size),inner=m.Ω_size)
if m.state_space == "KAΩ"                
    Xs = repeat(Π_base,outer=m.K_size*m.A_size)
elseif m.state_space == "KΩ"
    Xs = repeat(Π_base,outer=m.K_size)
end

#Is_Π, Js_Π = prepare_state_transition(m)
Is_Π, Js_Π, Xs_Π = prepare_state_transition(m)

m.options.ζ₁ = 2/3; m.options.Ζ₁ = 1.0
m.options.ζ₂ = 2/3; m.options.Ζ₂ = 1.0
m.options.ϵ₁ = 1e-10; m.options.ϵ₂ = 1e-10

Nᵖᵉ = 60
n_max = 200

m.cost_struct = 16

if m.cost_struct == 1 || m.cost_struct == 9 || m.cost_struct == 10 || m.cost_struct == 12
    ψ = 1.0;
    Γ = [0.0,0.0,0.0,0.0];
    κ = 0.0; ϕ = 0.0;
    ϕᶠ = 0.0;
elseif m.cost_struct == 2
    ψ = 1.0;
    Γ = [0.0,0.0,0.0,0.0,0.0];
    κ = 0.0; ϕ = 0.0;
    ϕᶠ = 0.0;
elseif m.cost_struct == 3
    ψ = 1.0;
    Γ = [0.0;0.0];
    κ = 0.0; ϕ = 0.0;
    ϕᶠ = 0.0;
elseif m.cost_struct == 4 || m.cost_struct == 8
    ψ = 1.0;
    Γ = [0.0;0.0];#;0.0
    κ = 0.0; ϕ = 0.0;
    ϕᶠ = 0.0;
elseif m.cost_struct == 5
    ψ = 1.0;
    Γ = [0.0,0.0,0.0,0.0];
    κ = 0.0; ϕ = 0.0;
    ϕᶠ = 0.0;
elseif m.cost_struct == 11
    ψ = 1.0;
    Γ = [0.0;0.0;0.0;0.0;0.0;0.0];
    κ = 0.0; ϕ = 0.0;      
    ϕᶠ = 0.0;
elseif m.cost_struct == 13
    ψ = 1.0;
    κ = 0.0; ϕ = 0.0;   
    ϕᶠ = 0.0;  
    Γ = [0.0;0.0;0.0;0.0];
elseif m.cost_struct == 14 || m.cost_struct == 15
    ψ = 1.0;
    κ = 0.0; ϕ = 0.0;   
    ϕᶠ = 0.0;  
    Γ = [0.0;0.0;0.0];
elseif m.cost_struct == 16
    ψ = 1.0;
    κ = 0.0; ϕ = 0.0;   
    ϕᶠ = 0.0;  
    Γ = [0.0;0.0;0.0;0.0;0.0];
elseif m.cost_struct == 17
    ψ = 1.0;
    κ = 0.0; ϕ = 0.0;   
    ϕᶠ = 0.0;  
    Γ = [0.0;0.0;0.0;0.0];
end  
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
timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit = build_data(df, d, m)

@time nfxp_golden_age = optimize_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d,m)
θ̂_golden_age = nfxp_golden_age.minimizer
ll_golden_age = nfxp_golden_age.minimum
writedlm(path*raw"/output/param_golden_age_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv",  θ̂_golden_age, ',')
SE_golden_age = calc_std_errors(θ̂_golden_age, d, m)
writedlm(path*raw"/output/se_golden_age_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", SE_golden_age, ',')

# Estimate: early
m.post = 0
d.game_start = 1804
d.game_end = 1830
timevar, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)

@time nfxp_early = optimize_log_likelihood(timevar, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_early = nfxp_early.minimizer
ll_early = nfxp_early.minimum
writedlm(path*raw"/output/param_early_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv",  θ̂_early, ',')
SE_early = calc_std_errors(θ̂_early, d, m)
writedlm(path*raw"/output/se_early_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", SE_early, ',')

# Estimate: post-petroleum
m.post = 1
d.game_start = 1859
d.game_end = 1910
timevar, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)

@time nfxp_post = optimize_log_likelihood(timevar, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_post = nfxp_post.minimizer
ll_post = nfxp_post.minimum
writedlm(path*raw"/output/param_post_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv",  θ̂_post, ',')
SE_post = calc_std_errors(θ̂_post, d, m)
writedlm(path*raw"/output/se_post_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", SE_post, ',')

writedlm(path*raw"/output/est_dynamic_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", 
                                                                [θ̂_early θ̂_golden_age θ̂_post; 
                                                                SE_early SE_golden_age SE_post;
                                                                ll_early ll_golden_age ll_post], ',')

#================================ End of estimation ================================#

