# Set the working directory
if Sys.iswindows()
    path = raw"C:\Users\yangk\Dropbox\Research\2022-Whaling-Regulation"
elseif Sys.islinux()
    #path = raw"/u/project/rafey/ykyun/2022-Whaling-Regulation/"
    path = raw"/mnt/phd/yyun/2022-Whaling-Regulation"
end

println("The number of threads is $(Threads.nthreads())")

# Load function call
println("Load functions and packages")
include(path*raw"/Julia/0_mutable_struct.jl")
include(path*raw"/Julia/1_build_up.jl")
include(path*raw"/Julia/2_model_estimation.jl")
include(path*raw"/Julia/3_model_simulation.jl")

println("Set data and model parameters")

m = DG_model(;  β=[], F="CD", 
                α=[],
                ρ=.9, xᵉ=1, δ=0.0, 
                ψ=1.0, Γ=[0.0;0.0], ϕ=0.0, κ=0.0, ϕᶠ=0.0,
                K_size=300, A_size=15, Ω_size=15, x_size=0,
                obj_incumbents_pre=0,obj_incumbents_post=0,
                scale_Q=100, scale_P=1000,
                N̄=400,
                nboot=100,
                state_space="KΩ")

d = DG_data(; t0=1804, t1=1920, 
                game_start=1858, game_end=1920,
                sMat=[], PVec=[], WVec=[], KVec=[], QVec=[], NVec=[], S=[], Q=[], πMat=[], VMat=[],
                ιMat0=[],χMat0=[],λVec0=[],ιMat1=[],χMat1=[],λVec1=[],Π=[])

Tbar = d.t1 - d.t0 + 1 # last time period of non-stationary transition

d.Π = []

#========================================================================================#
#---------------------------------------- Setup -----------------------------------------#
#========================================================================================#

df = CSV.read(path*raw"/Data/Intermediate/data_with_omega.csv", DataFrame) 

# Get static parameters
include(path*raw"/Julia/get_production_parameters.jl")
include(path*raw"/Julia/get_demand_parameters.jl")

# Construct state space
include(path*raw"/Julia/get_state_space.jl")

if m.state_space == "KAΩ"                
    println("The dimension of state space is $(m.K_size*m.A_size*m.Ω_size) with $(m.K_size) grid for K, $(m.A_size) grid for A, and $(m.Ω_size) grid for Ω")
elseif m.state_space == "KΩ"
    println("The dimension of state space is $(m.K_size*m.Ω_size) with $(m.K_size) grid for K and $(m.Ω_size) grid for Ω")
end

# Get aggregate variables
include(path*raw"/Julia/get_aggregate_variables.jl")

# Get initial guess of expected industry state in each period directly from data
include(path*raw"/Julia/get_initial_expected_industry_state.jl")

# Get entrants' intial state
include(path*raw"/Julia/get_entrants_state.jl")

# Get initial beliefs
#include(path*raw"/Julia/get_initial_beliefs.jl")

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

#m.ψ = 0.15008418457487882;m.κ = 106479.40144738233; m.ϕᶠ = 41866.05402976554; m.Γ = [217.44931690391547, 0.06377888157296054,87.59773939904876,0.0001473452416852835];
#θ = reduce(vcat,[m.ψ,m.κ,m.ϕᶠ,m.Γ])
#========================================================================================#
#-------------------------------------- Estimation --------------------------------------#
#========================================================================================#

#------------------------------------- With ψ = 1.0 -------------------------------------#

# Estimate: golden-age
m.post = 0
d.game_start = 1831 
d.game_end = 1858
timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit = build_data(df, d, m)

@time nfxp_golden_age = optimize_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d,m)
θ̂_golden_age = nfxp_golden_age.minimizer
ll_golden_age = nfxp_golden_age.minimum
writedlm(path*raw"/Output/param_golden_age_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv",  θ̂_golden_age, ',')
SE_golden_age = calc_std_errors(θ̂_golden_age, d, m)
writedlm(path*raw"/Output/se_golden_age_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", SE_golden_age, ',')

# Estimate: early
m.post = 0
d.game_start = 1804
d.game_end = 1830
timevar, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)

@time nfxp_early = optimize_log_likelihood(timevar, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_early = nfxp_early.minimizer
ll_early = nfxp_early.minimum
writedlm(path*raw"/Output/param_early_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv",  θ̂_early, ',')
SE_early = calc_std_errors(θ̂_early, d, m)
writedlm(path*raw"/Output/se_early_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", SE_early, ',')

# Estimate: post-petroleum
m.post = 1
d.game_start = 1859
d.game_end = 1910
timevar, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)

@time nfxp_post = optimize_log_likelihood(timevar, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_post = nfxp_post.minimizer
ll_post = nfxp_post.minimum
writedlm(path*raw"/Output/param_post_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv",  θ̂_post, ',')
SE_post = calc_std_errors(θ̂_post, d, m)
writedlm(path*raw"/Output/se_post_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", SE_post, ',')

writedlm(path*raw"/Output/est_dynamic_cost16_cap300increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100_240108_4.csv", 
                                                                [θ̂_early θ̂_golden_age θ̂_post; 
                                                                SE_early SE_golden_age SE_post;
                                                                ll_early ll_golden_age ll_post], ',')

#================================ End of estimation ================================#

#====================================== Simulation ======================================#

# Early period
m.post = 0
d.game_start = 1804
d.game_end = 1830
θ̂_early = CSV.read(path*raw"/Output/param_early_cost11_cap350increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100.csv", DataFrame, header=false)[1]
m.ψ = θ̂_early[1]
m.κ = θ̂_early[2]#; ϕᵉⁿ_early = θ̂_early[1]
m.ϕᶠ = θ̂_early[3]#; ϕᵉˣ_early = θ̂_early[2]
m.Γ = θ̂_early[4:end]; 
exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE(d, m)

#m.Γ[1] = .01
#m.Γ[2] = 170
#m.Γ[3] = 30

# Golden-age
m.post = 0
d.game_start = 1831
d.game_end = 1858
θ̂_golden_age = CSV.read(path*raw"/Output/param_golden_age_cost1_cap400increm_age1_prod7_NM_w_bound_w_ec_ver10_scale100_240102_5.csv", DataFrame, header=false)[1]
m.ψ = θ̂_golden_age[1]
m.κ = θ̂_golden_age[2]#; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
m.ϕᶠ = θ̂_golden_age[3]#; ϕᵉˣ_golden_age = θ̂_golden_age[2]
m.Γ = θ̂_golden_age[4:end];
exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE(d, m)

# Post-petroleum period
m.post = 1
d.game_start = 1859
d.game_end = 1910
θ̂_post = CSV.read(path*raw"/Output/param_post_cost11_cap350increm_age1_prod15_NM_w_bound_w_ec_ver10_scale100.csv", DataFrame, header=false)[1]
m.κ = θ̂_post[1]#; ϕᵉⁿ_post = θ̂_post[1]
m.ϕᶠ = θ̂_post[2]#; ϕᵉˣ_post = θ̂_post[2]
m.Γ = θ̂_post[3:end];
exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE(d, m)

Random.seed!(1234)
nsim = 10_000
Nsim = zeros(Tbar+1,nsim); Ksim = zeros(Tbar+1,nsim); Qsim = zeros(Tbar+1,nsim); Wsim = zeros(Tbar+1,nsim); Psim = zeros(Tbar+1,nsim); PSsim = zeros(Tbar+1,nsim); CSsim = zeros(Tbar+1,nsim); SSsim = zeros(Tbar+1,nsim)
Threads.@threads for i = 1:nsim
    Nsim[:,i], Ksim[:,i], Qsim[:,i], Wsim[:,i], Psim[:,i], PSsim[:,i], CSsim[:,i], SSsim[:,i] = simulate_data(ϕᵉⁿ_early, ϕᵉˣ_early, γⁱⁿᵛ_early, γᵈⁱᵛ_early, 
                                                                                                            ϕᵉⁿ_golden_age, ϕᵉˣ_golden_age, γⁱⁿᵛ_golden_age, γᵈⁱᵛ_golden_age,
                                                                                                            ϕᵉⁿ_post, ϕᵉˣ_post, γⁱⁿᵛ_post, γᵈⁱᵛ_post,
                                                                                                            exit_prob_early, invest_prob_early, entry_prob_early, exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, exit_prob_post, invest_prob_post, entry_prob_post, 
                                                                                                            df, d, m)
end

writedlm(path*raw"/Output/Nsim.csv", Nsim, ',')
writedlm(path*raw"/Output/Ksim.csv", Ksim, ',')
writedlm(path*raw"/Output/Qsim.csv", Qsim, ',')
writedlm(path*raw"/Output/Wsim.csv", Wsim, ',')
writedlm(path*raw"/Output/Psim.csv", Psim, ',')
writedlm(path*raw"/Output/PSsim.csv", PSsim, ',')
writedlm(path*raw"/Output/CSsim.csv", CSsim, ',')
writedlm(path*raw"/Output/SSsim.csv", SSsim, ',')

Nsim_mean = mean(Nsim,dims=2); Nsim_sd = std(Nsim,dims=2); Nsim_lo = Nsim_mean - 1.96*Nsim_sd; Nsim_hi = Nsim_mean + 1.96*Nsim_sd
Ksim_mean = mean(Ksim,dims=2); Ksim_sd = std(Ksim,dims=2); Ksim_lo = Ksim_mean - 1.96*Ksim_sd; Ksim_hi = Ksim_mean + 1.96*Ksim_sd
Qsim_mean = mean(Qsim,dims=2); Qsim_sd = std(Qsim,dims=2); Qsim_lo = Qsim_mean - 1.96*Qsim_sd; Qsim_hi = Qsim_mean + 1.96*Qsim_sd
Wsim_mean = mean(Wsim,dims=2); Wsim_sd = std(Wsim,dims=2); Wsim_lo = Wsim_mean - 1.96*Wsim_sd; Wsim_hi = Wsim_mean + 1.96*Wsim_sd
Psim_mean = mean(Psim,dims=2); Psim_sd = std(Psim,dims=2); Psim_lo = Psim_mean - 1.96*Psim_sd; Psim_hi = Psim_mean + 1.96*Psim_sd

# Model fit
Qfig1, Qfig2 = model_fit_visualize(QData./1000, Qsim_mean./1000, Qsim_lo./1000, Qsim_hi./1000, d)
Qfig1 = plot(Qfig1, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_data.pdf")
Qfig2 = plot(Qfig2, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_model.pdf")

Pfig1, Pfig2 = model_fit_visualize(PData./1000, Psim_mean./1000, Psim_lo./1000, Psim_hi./1000, d)
Pfig1 = plot(Pfig1, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_data.pdf")
Pfig2 = plot(Pfig2, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_model.pdf")

Nfig1, Nfig2 = model_fit_visualize(NData, Nsim_mean, Nsim_lo, Nsim_hi, d)
Nfig1 = plot(Nfig1, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_data.pdf")
Nfig2 = plot(Nfig2, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_model.pdf")

Kfig1, Kfig2 = model_fit_visualize(KData./1000, Ksim_mean./1000, Ksim_lo./1000, Ksim_hi./1000, d)
Kfig1 = plot(Kfig1, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_data.pdf")
Kfig2 = plot(Kfig2, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_model.pdf")

Wfig1, Wfig2 = model_fit_visualize(WData./1_000_000, Wsim_mean./1_000_000, Wsim_lo./1_000_000, Wsim_hi./1_000_000, d)
Wfig1 = plot(Wfig1, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_data.pdf")
Wfig2 = plot(Wfig2, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_model.pdf")




#----------------------------------- No strategic interaction ----------------------------------#

# Early period
m.post = 0
d.game_start = 1804
d.game_end = 1830
θ̂_early = CSV.read(path*raw"/Output/est_params_dynamic_early_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_early[1]; ϕᵉⁿ_early = θ̂_early[1]
m.ϕ = θ̂_early[2]; ϕᵉˣ_early = θ̂_early[2]
m.Γ = θ̂_early[3:end]; γⁱⁿᵛ_early = θ̂_early[3]; γᵈⁱᵛ_early = θ̂_early[4]
exit_prob_early = zeros(Tbar+1,m.x_size)
invest_prob_early = zeros(Tbar+1,m.x_size,m.K_size)
entry_prob_early = zeros(Tbar+1)
for year = d.game_start-d.t0+1:d.game_end-d.t0+1
    exit_prob, invest_prob, entry_prob, NVec, KVec, QVec, WVec, PVec = solve_NOE_NSI(Γ, ϕ, κ, ψ, d, m, year)
    exit_prob_early[year,:] = exit_prob[year,:]
    invest_prob_early[year,:,:] = invest_prob[year,:,:]
    entry_prob_early[year] = entry_prob[year]
end
    
# Golden age period
m.post = 0
d.game_start = 1831
d.game_end = 1858
θ̂_golden_age = CSV.read(path*raw"/Output/est_params_dynamic_golden_age_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_golden_age[1]; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
m.ϕ = θ̂_golden_age[2]; ϕᵉˣ_golden_age = θ̂_golden_age[2]
m.Γ = θ̂_golden_age[3:end]; γⁱⁿᵛ_golden_age = θ̂_golden_age[3]; γᵈⁱᵛ_golden_age = θ̂_golden_age[4]
exit_prob_golden_age = zeros(Tbar+1,m.x_size)
invest_prob_golden_age = zeros(Tbar+1,m.x_size,m.K_size)
entry_prob_golden_age = zeros(Tbar+1)
for year = d.game_start-d.t0+1:d.game_end-d.t0+1
    exit_prob, invest_prob, entry_prob, NVec, KVec, QVec, WVec, PVec = solve_NOE_NSI(Γ, ϕ, κ, ψ, d, m, year)
    exit_prob_golden_age[year,:] = exit_prob[year,:]
    invest_prob_golden_age[year,:,:] = invest_prob[year,:,:]
    entry_prob_golden_age[year] = entry_prob[year]
end

# Post-petroleum period
m.post = 1
d.game_start = 1859
d.game_end = 1910
θ̂_post = CSV.read(path*raw"/Output/est_params_dynamic_post_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_post[1]; ϕᵉⁿ_post = θ̂_post[1]
m.ϕ = θ̂_post[2]; ϕᵉˣ_post = θ̂_post[2]
m.Γ = θ̂_post[3:end]; γⁱⁿᵛ_post = θ̂_post[3]; ; γᵈⁱᵛ_post = θ̂_post[4]
exit_prob_post = zeros(Tbar+1,m.x_size)
invest_prob_post = zeros(Tbar+1,m.x_size,m.K_size)
entry_prob_post = zeros(Tbar+1)
for year = d.game_start-d.t0+1:d.game_end-d.t0+1
    exit_prob, invest_prob, entry_prob, NVec, KVec, QVec, WVec, PVec = solve_NOE_NSI(Γ, ϕ, κ, ψ, d, m, year)
    exit_prob_post[year,:] = exit_prob[year,:]
    invest_prob_post[year,:,:] = invest_prob[year,:,:]
    entry_prob_post[year] = entry_prob[year]
end

Random.seed!(1234)
nsim = 1_00
Nsim = zeros(Tbar+1,nsim); Ksim = zeros(Tbar+1,nsim); Qsim = zeros(Tbar+1,nsim); Wsim = zeros(Tbar+1,nsim); Psim = zeros(Tbar+1,nsim); PSsim = zeros(Tbar+1,nsim); CSsim = zeros(Tbar+1,nsim); SSsim = zeros(Tbar+1,nsim)
Threads.@threads for i = 1:nsim
    Nsim[:,i], Ksim[:,i], Qsim[:,i], Wsim[:,i], Psim[:,i], PSsim[:,i], CSsim[:,i], SSsim[:,i] = simulate_data(df, d, m)
end

Nsim_mean = mean(Nsim,dims=2); Nsim_sd = std(Nsim,dims=2); Nsim_lo = Nsim_mean - 1.96*Nsim_sd; Nsim_hi = Nsim_mean + 1.96*Nsim_sd
Ksim_mean = mean(Ksim,dims=2); Ksim_sd = std(Ksim,dims=2); Ksim_lo = Ksim_mean - 1.96*Ksim_sd; Ksim_hi = Ksim_mean + 1.96*Ksim_sd
Qsim_mean = mean(Qsim,dims=2); Qsim_sd = std(Qsim,dims=2); Qsim_lo = Qsim_mean - 1.96*Qsim_sd; Qsim_hi = Qsim_mean + 1.96*Qsim_sd
Wsim_mean = mean(Wsim,dims=2); Wsim_sd = std(Wsim,dims=2); Wsim_lo = Wsim_mean - 1.96*Wsim_sd; Wsim_hi = Wsim_mean + 1.96*Wsim_sd
Psim_mean = mean(Psim,dims=2); Psim_sd = std(Psim,dims=2); Psim_lo = Psim_mean - 1.96*Psim_sd; Psim_hi = Psim_mean + 1.96*Psim_sd

# Model fit
Qfig1, Qfig2 = model_fit_visualize(QData./1000, Qsim_mean./1000, Qsim_lo./1000, Qsim_hi./1000, d)
Qfig1 = plot(Qfig1, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_data_NSI.pdf")
Qfig2 = plot(Qfig2, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_model_NSI.pdf")

Pfig1, Pfig2 = model_fit_visualize(PData./1000, Psim_mean./1000, Psim_lo./1000, Psim_hi./1000, d)
Pfig1 = plot(Pfig1, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_data_NSI.pdf")
Pfig2 = plot(Pfig2, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_model_NSI.pdf")

Nfig1, Nfig2 = model_fit_visualize(NData, Nsim_mean, Nsim_lo, Nsim_hi, d)
Nfig1 = plot(Nfig1, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_data_NSI.pdf")
Nfig2 = plot(Nfig2, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_model_NSI.pdf")

Kfig1, Kfig2 = model_fit_visualize(KData./1000, Ksim_mean./1000, Ksim_lo./1000, Ksim_hi./1000, d)
Kfig1 = plot(Kfig1, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_data_NSI.pdf")
Kfig2 = plot(Kfig2, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_model_NSI.pdf")

Wfig1, Wfig2 = model_fit_visualize(WData./1_000_000, Wsim_mean./1_000_000, Wsim_lo./1_000_000, Wsim_hi./1_000_000, d)
Wfig1 = plot(Wfig1, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_data_NSI.pdf")
Wfig2 = plot(Wfig2, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_model_NSI.pdf")

annual_PS_NSI = mean(PSsim,dims=2)
annual_CS_NSI = mean(CSsim,dims=2)
annual_SS_NSI = mean(SSsim,dims=2)
total_PS_NSI = mean(sum(PSsim,dims=1))
total_CS_NSI = mean(sum(CSsim,dims=1))
total_SS_NSI = mean(sum(SSsim,dims=1))


#----------------------------------- No endogenous price ----------------------------------#

# Early period
m.post = 0
d.game_start = 1804
d.game_end = 1830
θ̂_early = CSV.read(path*raw"/Output/est_params_dynamic_early_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_early[1]; ϕᵉⁿ_early = θ̂_early[1]
m.ϕ = θ̂_early[2]; ϕᵉˣ_early = θ̂_early[2]
m.Γ = θ̂_early[3:end]; γⁱⁿᵛ_early = θ̂_early[3]; γᵈⁱᵛ_early = θ̂_early[4]
exit_prob_early = zeros(Tbar+1,m.x_size)
invest_prob_early = zeros(Tbar+1,m.x_size,m.K_size)
entry_prob_early = zeros(Tbar+1)
#for year = d.game_start-d.t0+1:d.game_end-d.t0+1
exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE_NEP(Γ, ϕ, κ, ψ, d, m, year)
#    exit_prob_early[year,:] = exit_prob[year,:]
#    invest_prob_early[year,:,:] = invest_prob[year,:,:]
#    entry_prob_early[year] = entry_prob[year]
#end
#exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE_NSI_NEP(Γ, ϕ, κ, ψ, d, m, year)
    
# Golden age period
m.post = 0
d.game_start = 1831
d.game_end = 1858
θ̂_golden_age = CSV.read(path*raw"/Output/est_params_dynamic_golden_age_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_golden_age[1]; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
m.ϕ = θ̂_golden_age[2]; ϕᵉˣ_golden_age = θ̂_golden_age[2]
m.Γ = θ̂_golden_age[3:end]; γⁱⁿᵛ_golden_age = θ̂_golden_age[3]; γᵈⁱᵛ_golden_age = θ̂_golden_age[4]
exit_prob_golden_age = zeros(Tbar+1,m.x_size)
invest_prob_golden_age = zeros(Tbar+1,m.x_size,m.K_size)
entry_prob_golden_age = zeros(Tbar+1)
#for year = d.game_start-d.t0+1:d.game_end-d.t0+1
exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE_NEP(Γ, ϕ, κ, ψ, d, m, year)
#    exit_prob_golden_age[year,:] = exit_prob[year,:]
#    invest_prob_golden_age[year,:,:] = invest_prob[year,:,:]
#    entry_prob_golden_age[year] = entry_prob[year]
#end
#exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE_NSI_NEP(Γ, ϕ, κ, ψ, d, m, year)

# Post-petroleum period
m.post = 1
d.game_start = 1859
d.game_end = 1910
θ̂_post = CSV.read(path*raw"/Output/est_params_dynamic_post_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_post[1]; ϕᵉⁿ_post = θ̂_post[1]
m.ϕ = θ̂_post[2]; ϕᵉˣ_post = θ̂_post[2]
m.Γ = θ̂_post[3:end]; γⁱⁿᵛ_post = θ̂_post[3]; ; γᵈⁱᵛ_post = θ̂_post[4]
exit_prob_post = zeros(Tbar+1,m.x_size)
invest_prob_post = zeros(Tbar+1,m.x_size,m.K_size)
entry_prob_post = zeros(Tbar+1)
#for year = d.game_start-d.t0+1:d.game_end-d.t0+1
exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE_NEP(Γ, ϕ, κ, ψ, d, m, year)
#    exit_prob_post[year,:] = exit_prob[year,:]
#    invest_prob_post[year,:,:] = invest_prob[year,:,:]
#    entry_prob_post[year] = entry_prob[year]
#end
#exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE_NSI_NEP(Γ, ϕ, κ, ψ, d, m, year)

Random.seed!(1234)
nsim = 1_00
Nsim = zeros(Tbar+1,nsim); Ksim = zeros(Tbar+1,nsim); Qsim = zeros(Tbar+1,nsim); Wsim = zeros(Tbar+1,nsim); Psim = zeros(Tbar+1,nsim); PSsim = zeros(Tbar+1,nsim); CSsim = zeros(Tbar+1,nsim); SSsim = zeros(Tbar+1,nsim)
Threads.@threads for i = 1:nsim
    Nsim[:,i], Ksim[:,i], Qsim[:,i], Wsim[:,i], Psim[:,i], PSsim[:,i], CSsim[:,i], SSsim[:,i] = simulate_data(ϕᵉⁿ_early, ϕᵉˣ_early, γⁱⁿᵛ_early, γᵈⁱᵛ_early, 
                                                                                                            ϕᵉⁿ_golden_age, ϕᵉˣ_golden_age, γⁱⁿᵛ_golden_age, γᵈⁱᵛ_golden_age,
                                                                                                            ϕᵉⁿ_post, ϕᵉˣ_post, γⁱⁿᵛ_post, γᵈⁱᵛ_post,
                                                                                                            exit_prob_early, invest_prob_early, entry_prob_early, 
                                                                                                            exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, 
                                                                                                            exit_prob_post, invest_prob_post, entry_prob_post, 
                                                                                                            df, d, m)
end

Nsim_mean = mean(Nsim,dims=2); Nsim_sd = std(Nsim,dims=2); Nsim_lo = Nsim_mean - 1.96*Nsim_sd; Nsim_hi = Nsim_mean + 1.96*Nsim_sd
Ksim_mean = mean(Ksim,dims=2); Ksim_sd = std(Ksim,dims=2); Ksim_lo = Ksim_mean - 1.96*Ksim_sd; Ksim_hi = Ksim_mean + 1.96*Ksim_sd
Qsim_mean = mean(Qsim,dims=2); Qsim_sd = std(Qsim,dims=2); Qsim_lo = Qsim_mean - 1.96*Qsim_sd; Qsim_hi = Qsim_mean + 1.96*Qsim_sd
Wsim_mean = mean(Wsim,dims=2); Wsim_sd = std(Wsim,dims=2); Wsim_lo = Wsim_mean - 1.96*Wsim_sd; Wsim_hi = Wsim_mean + 1.96*Wsim_sd
Psim_mean = mean(Psim,dims=2); Psim_sd = std(Psim,dims=2); Psim_lo = Psim_mean - 1.96*Psim_sd; Psim_hi = Psim_mean + 1.96*Psim_sd

# Model fit
Qfig1, Qfig2 = model_fit_visualize(QData./1000, Qsim_mean./1000, Qsim_lo./1000, Qsim_hi./1000, d)
Qfig1 = plot(Qfig1, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_data_NSI.pdf")
Qfig2 = plot(Qfig2, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_model_NSI.pdf")

Pfig1, Pfig2 = model_fit_visualize(PData./1000, Psim_mean./1000, Psim_lo./1000, Psim_hi./1000, d)
Pfig1 = plot(Pfig1, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_data_NSI.pdf")
Pfig2 = plot(Pfig2, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_model_NSI.pdf")

Nfig1, Nfig2 = model_fit_visualize(NData, Nsim_mean, Nsim_lo, Nsim_hi, d)
Nfig1 = plot(Nfig1, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_data_NSI.pdf")
Nfig2 = plot(Nfig2, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_model_NSI.pdf")

Kfig1, Kfig2 = model_fit_visualize(KData./1000, Ksim_mean./1000, Ksim_lo./1000, Ksim_hi./1000, d)
Kfig1 = plot(Kfig1, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_data_NSI.pdf")
Kfig2 = plot(Kfig2, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_model_NSI.pdf")

Wfig1, Wfig2 = model_fit_visualize(WData./1_000_000, Wsim_mean./1_000_000, Wsim_lo./1_000_000, Wsim_hi./1_000_000, d)
Wfig1 = plot(Wfig1, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_data_NSI.pdf")
Wfig2 = plot(Wfig2, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_model_NSI.pdf")

annual_PS_NSI = mean(PSsim,dims=2)
annual_CS_NSI = mean(CSsim,dims=2)
annual_SS_NSI = mean(SSsim,dims=2)
total_PS_NSI = mean(sum(PSsim,dims=1))
total_CS_NSI = mean(sum(CSsim,dims=1))
total_SS_NSI = mean(sum(SSsim,dims=1))


#------------------------------------- Policy experiments: Tax -------------------------------------#
# Social discount factor
discount = (.97).^(collect(1:Tbar+1) .- 1)

# Read the baseline simulated paths without regulation
PSsim_base = CSV.File(path*raw"/Output/PSsim.csv"; header=false) |> Tables.matrix
CSsim_base = CSV.File(path*raw"/Output/CSsim.csv"; header=false) |> Tables.matrix
SSsim_base = CSV.File(path*raw"/Output/SSsim.csv"; header=false) |> Tables.matrix

# Net present value of surpluses
PSsim_base = PSsim_base.*discount
CSsim_base = CSsim_base.*discount
SSsim_base = SSsim_base.*discount

total_PS_base = mean(sum(PSsim_base,dims=1))
total_CS_base = mean(sum(CSsim_base,dims=1))
total_SS_base = mean(sum(SSsim_base,dims=1))

# Setting for tax
tax_lo = 1.00; tax_hi = 1.80
tax_width = 0.05
tax_num = length(collect(tax_lo:tax_width:tax_hi))

total_PS_entry_tax = zeros(tax_num)
total_CS_entry_tax = zeros(tax_num)
total_SS_entry_tax = zeros(tax_num)

total_PS_inv_tax = zeros(tax_num)
total_CS_inv_tax = zeros(tax_num)
total_SS_inv_tax = zeros(tax_num)

total_PS_pigouvian_tax = zeros(tax_num)
total_CS_pigouvian_tax = zeros(tax_num)
total_SS_pigouvian_tax = zeros(tax_num)

Nsim_mean = zeros(Tbar+1,tax_num); Nsim_sd = zeros(Tbar+1,tax_num); Nsim_lo = zeros(Tbar+1,tax_num); Nsim_hi = zeros(Tbar+1,tax_num)
Ksim_mean = zeros(Tbar+1,tax_num); Ksim_sd = zeros(Tbar+1,tax_num); Ksim_lo = zeros(Tbar+1,tax_num); Ksim_hi = zeros(Tbar+1,tax_num)
Qsim_mean = zeros(Tbar+1,tax_num); Qsim_sd = zeros(Tbar+1,tax_num); Qsim_lo = zeros(Tbar+1,tax_num); Qsim_hi = zeros(Tbar+1,tax_num)
Wsim_mean = zeros(Tbar+1,tax_num); Wsim_sd = zeros(Tbar+1,tax_num); Wsim_lo = zeros(Tbar+1,tax_num); Wsim_hi = zeros(Tbar+1,tax_num)
Psim_mean = zeros(Tbar+1,tax_num); Psim_sd = zeros(Tbar+1,tax_num); Psim_lo = zeros(Tbar+1,tax_num); Psim_hi = zeros(Tbar+1,tax_num)

# Entry tax
i = 0
for tax_rate = tax_lo:tax_width:tax_hi
    i += 1
    # Early period
    m.post = 0
    d.game_start = 1804
    d.game_end = 1830
    θ̂_early = CSV.read(path*raw"/Output/est_params_dynamic_early_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_early[1]; ϕᵉⁿ_early = θ̂_early[1]
    m.ϕ = θ̂_early[2]; ϕᵉˣ_early = θ̂_early[2]
    m.Γ = θ̂_early[3:end]; γⁱⁿᵛ_early = θ̂_early[3]; γᵈⁱᵛ_early = θ̂_early[4]
    m.κ *= tax_rate; ϕᵉⁿ_early *= tax_rate
    exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

    # Golden-age
    m.post = 0
    d.game_start = 1831
    d.game_end = 1858
    θ̂_golden_age = CSV.read(path*raw"/Output/est_params_dynamic_golden_age_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_golden_age[1]; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
    m.ϕ = θ̂_golden_age[2]; ϕᵉˣ_golden_age = θ̂_golden_age[2]
    m.Γ = θ̂_golden_age[3:end]; γⁱⁿᵛ_golden_age = θ̂_golden_age[3]; ; γᵈⁱᵛ_golden_age = θ̂_golden_age[4]
    m.κ *= tax_rate; ϕᵉⁿ_golden_age *= tax_rate
    exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

    # Post-petroleum period
    m.post = 1
    d.game_start = 1859
    d.game_end = 1910
    θ̂_post = CSV.read(path*raw"/Output/est_params_dynamic_post_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_post[1]; ϕᵉⁿ_post = θ̂_post[1]
    m.ϕ = θ̂_post[2]; ϕᵉˣ_post = θ̂_post[2]
    m.Γ = θ̂_post[3:end]; γⁱⁿᵛ_post = θ̂_post[3]; ; γᵈⁱᵛ_post = θ̂_post[4]
    m.κ *= tax_rate; ϕᵉⁿ_post *= tax_rate
    exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

    Random.seed!(1234)
    nsim = 1_000
    Nsim = zeros(Tbar+1,nsim); Ksim = zeros(Tbar+1,nsim); Qsim = zeros(Tbar+1,nsim); Wsim = zeros(Tbar+1,nsim); Psim = zeros(Tbar+1,nsim); PSsim = zeros(Tbar+1,nsim); CSsim = zeros(Tbar+1,nsim); SSsim = zeros(Tbar+1,nsim)
    Threads.@threads for i = 1:nsim
        #println("sim: $i")
        Nsim[:,i], Ksim[:,i], Qsim[:,i], Wsim[:,i], Psim[:,i], PSsim[:,i], CSsim[:,i], SSsim[:,i] = simulate_data(ϕᵉⁿ_early, ϕᵉˣ_early, γⁱⁿᵛ_early, γᵈⁱᵛ_early, 
                                                                                                                ϕᵉⁿ_golden_age, ϕᵉˣ_golden_age, γⁱⁿᵛ_golden_age, γᵈⁱᵛ_golden_age,
                                                                                                                ϕᵉⁿ_post, ϕᵉˣ_post, γⁱⁿᵛ_post, γᵈⁱᵛ_post,
                                                                                                                exit_prob_early, invest_prob_early, entry_prob_early, exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, exit_prob_post, invest_prob_post, entry_prob_post, 
                                                                                                                df, d, m)
    end

    Nsim_mean[:,i] = mean(Nsim,dims=2); Nsim_sd[:,i] = std(Nsim,dims=2); Nsim_lo[:,i] = Nsim_mean[:,i] - 1.96*Nsim_sd[:,i]; Nsim_hi[:,i] = Nsim_mean[:,i] + 1.96*Nsim_sd[:,i]
    Ksim_mean[:,i] = mean(Ksim,dims=2); Ksim_sd[:,i] = std(Ksim,dims=2); Ksim_lo[:,i] = Ksim_mean[:,i] - 1.96*Ksim_sd[:,i]; Ksim_hi[:,i] = Ksim_mean[:,i] + 1.96*Ksim_sd[:,i]
    Qsim_mean[:,i] = mean(Qsim,dims=2); Qsim_sd[:,i] = std(Qsim,dims=2); Qsim_lo[:,i] = Qsim_mean[:,i] - 1.96*Qsim_sd[:,i]; Qsim_hi[:,i] = Qsim_mean[:,i] + 1.96*Qsim_sd[:,i]
    Wsim_mean[:,i] = mean(Wsim,dims=2); Wsim_sd[:,i] = std(Wsim,dims=2); Wsim_lo[:,i] = Wsim_mean[:,i] - 1.96*Wsim_sd[:,i]; Wsim_hi[:,i] = Wsim_mean[:,i] + 1.96*Wsim_sd[:,i]
    Psim_mean[:,i] = mean(Psim,dims=2); Psim_sd[:,i] = std(Psim,dims=2); Psim_lo[:,i] = Psim_mean[:,i] - 1.96*Psim_sd[:,i]; Psim_hi[:,i] = Psim_mean[:,i] + 1.96*Psim_sd[:,i]

    PSsim = PSsim.*discount
    CSsim = CSsim.*discount
    SSsim = SSsim.*discount

    total_PS_entry_tax[i] = mean(sum(PSsim,dims=1))
    total_CS_entry_tax[i] = mean(sum(CSsim,dims=1))
    total_SS_entry_tax[i] = mean(sum(SSsim,dims=1))
end

# Investment tax
i = 0
for tax_rate = tax_lo:tax_width:tax_hi
    i += 1
    # Early period
    m.post = 0
    d.game_start = 1804
    d.game_end = 1830
    θ̂_early = CSV.read(path*raw"/Output/est_params_dynamic_early_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_early[1]; ϕᵉⁿ_early = θ̂_early[1]
    m.ϕ = θ̂_early[2]; ϕᵉˣ_early = θ̂_early[2]
    m.Γ = θ̂_early[3:end]; γⁱⁿᵛ_early = θ̂_early[3]; γᵈⁱᵛ_early = θ̂_early[4]
    m.Γ[1] *= tax_rate; γⁱⁿᵛ_early *= tax_rate
    exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

    # Golden-age
    m.post = 0
    d.game_start = 1831
    d.game_end = 1858
    θ̂_golden_age = CSV.read(path*raw"/Output/est_params_dynamic_golden_age_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_golden_age[1]; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
    m.ϕ = θ̂_golden_age[2]; ϕᵉˣ_golden_age = θ̂_golden_age[2]
    m.Γ = θ̂_golden_age[3:end]; γⁱⁿᵛ_golden_age = θ̂_golden_age[3]; γᵈⁱᵛ_golden_age = θ̂_golden_age[4]
    m.Γ[1] *= tax_rate; γⁱⁿᵛ_golden_age *= tax_rate
    exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

    # Post-petroleum period
    m.post = 1
    d.game_start = 1859
    d.game_end = 1910
    θ̂_post = CSV.read(path*raw"/Output/est_params_dynamic_post_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_post[1]; ϕᵉⁿ_post = θ̂_post[1]
    m.ϕ = θ̂_post[2]; ϕᵉˣ_post = θ̂_post[2]
    m.Γ = θ̂_post[3:end]; γⁱⁿᵛ_post = θ̂_post[3]; γᵈⁱᵛ_post = θ̂_post[4]
    m.Γ[1] *= tax_rate; γⁱⁿᵛ_post *= tax_rate
    exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

    Random.seed!(1234)
    nsim = 1_000
    Nsim = zeros(Tbar+1,nsim); Ksim = zeros(Tbar+1,nsim); Qsim = zeros(Tbar+1,nsim); Wsim = zeros(Tbar+1,nsim); Psim = zeros(Tbar+1,nsim); PSsim = zeros(Tbar+1,nsim); CSsim = zeros(Tbar+1,nsim); SSsim = zeros(Tbar+1,nsim)
    Threads.@threads for i = 1:nsim
        #println("sim: $i")
        Nsim[:,i], Ksim[:,i], Qsim[:,i], Wsim[:,i], Psim[:,i], PSsim[:,i], CSsim[:,i], SSsim[:,i] = simulate_data(ϕᵉⁿ_early, ϕᵉˣ_early, γⁱⁿᵛ_early, γᵈⁱᵛ_early, 
                                                                                                            ϕᵉⁿ_golden_age, ϕᵉˣ_golden_age, γⁱⁿᵛ_golden_age, γᵈⁱᵛ_golden_age,
                                                                                                            ϕᵉⁿ_post, ϕᵉˣ_post, γⁱⁿᵛ_post, γᵈⁱᵛ_post,
                                                                                                            exit_prob_early, invest_prob_early, entry_prob_early, exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, exit_prob_post, invest_prob_post, entry_prob_post, 
                                                                                                            df, d, m)
    end

    Nsim_mean[:,i] = mean(Nsim,dims=2); Nsim_sd[:,i] = std(Nsim,dims=2); Nsim_lo[:,i] = Nsim_mean[:,i] - 1.96*Nsim_sd[:,i]; Nsim_hi[:,i] = Nsim_mean[:,i] + 1.96*Nsim_sd[:,i]
    Ksim_mean[:,i] = mean(Ksim,dims=2); Ksim_sd[:,i] = std(Ksim,dims=2); Ksim_lo[:,i] = Ksim_mean[:,i] - 1.96*Ksim_sd[:,i]; Ksim_hi[:,i] = Ksim_mean[:,i] + 1.96*Ksim_sd[:,i]
    Qsim_mean[:,i] = mean(Qsim,dims=2); Qsim_sd[:,i] = std(Qsim,dims=2); Qsim_lo[:,i] = Qsim_mean[:,i] - 1.96*Qsim_sd[:,i]; Qsim_hi[:,i] = Qsim_mean[:,i] + 1.96*Qsim_sd[:,i]
    Wsim_mean[:,i] = mean(Wsim,dims=2); Wsim_sd[:,i] = std(Wsim,dims=2); Wsim_lo[:,i] = Wsim_mean[:,i] - 1.96*Wsim_sd[:,i]; Wsim_hi[:,i] = Wsim_mean[:,i] + 1.96*Wsim_sd[:,i]
    Psim_mean[:,i] = mean(Psim,dims=2); Psim_sd[:,i] = std(Psim,dims=2); Psim_lo[:,i] = Psim_mean[:,i] - 1.96*Psim_sd[:,i]; Psim_hi[:,i] = Psim_mean[:,i] + 1.96*Psim_sd[:,i]

    PSsim = PSsim.*discount
    CSsim = CSsim.*discount
    SSsim = SSsim.*discount

    total_PS_inv_tax[i] = mean(sum(PSsim,dims=1))
    total_CS_inv_tax[i] = mean(sum(CSsim,dims=1))
    total_SS_inv_tax[i] = mean(sum(SSsim,dims=1))
end

# Pigouvian tax per whale
i = 0
for tax_rate = tax_lo:tax_width:tax_hi
    i += 1
    # Early period
    m.post = 0
    d.game_start = 1804
    d.game_end = 1830
    θ̂_early = CSV.read(path*raw"/Output/est_params_dynamic_early_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_early[1]; ϕᵉⁿ_early = θ̂_early[1]
    m.ϕ = θ̂_early[2]; ϕᵉˣ_early = θ̂_early[2]
    m.Γ = θ̂_early[3:end]; γⁱⁿᵛ_early = θ̂_early[3]; γᵈⁱᵛ_early = θ̂_early[4]
    m.Γ[1] *= tax_rate; γⁱⁿᵛ_early *= tax_rate
    exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE_pigouvian_tax(Γ, ϕ, κ, ψ, d, m, tax_rate)

    # Golden-age
    m.post = 0
    d.game_start = 1831
    d.game_end = 1858
    θ̂_golden_age = CSV.read(path*raw"/Output/est_params_dynamic_golden_age_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_golden_age[1]; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
    m.ϕ = θ̂_golden_age[2]; ϕᵉˣ_golden_age = θ̂_golden_age[2]
    m.Γ = θ̂_golden_age[3:end]; γⁱⁿᵛ_golden_age = θ̂_golden_age[3]; γᵈⁱᵛ_golden_age = θ̂_golden_age[4]
    m.Γ[1] *= tax_rate; γⁱⁿᵛ_golden_age *= tax_rate
    exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE_pigouvian_tax(Γ, ϕ, κ, ψ, d, m, tax_rate)

    # Post-petroleum period
    m.post = 1
    d.game_start = 1859
    d.game_end = 1910
    θ̂_post = CSV.read(path*raw"/Output/est_params_dynamic_post_cost11_cap25.csv", DataFrame, header=false)[1]
    m.κ = θ̂_post[1]; ϕᵉⁿ_post = θ̂_post[1]
    m.ϕ = θ̂_post[2]; ϕᵉˣ_post = θ̂_post[2]
    m.Γ = θ̂_post[3:end]; γⁱⁿᵛ_post = θ̂_post[3]; γᵈⁱᵛ_post = θ̂_post[4]
    m.Γ[1] *= tax_rate; γⁱⁿᵛ_post *= tax_rate
    exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE_pigouvian_tax(Γ, ϕ, κ, ψ, d, m, tax_rate)

    Random.seed!(1234)
    nsim = 1_000
    Nsim = zeros(Tbar+1,nsim); Ksim = zeros(Tbar+1,nsim); Qsim = zeros(Tbar+1,nsim); Wsim = zeros(Tbar+1,nsim); Psim = zeros(Tbar+1,nsim); PSsim = zeros(Tbar+1,nsim); CSsim = zeros(Tbar+1,nsim); SSsim = zeros(Tbar+1,nsim)
    Threads.@threads for i = 1:nsim
        #println("sim: $i")
        Nsim[:,i], Ksim[:,i], Qsim[:,i], Wsim[:,i], Psim[:,i], PSsim[:,i], CSsim[:,i], SSsim[:,i] = simulate_data(ϕᵉⁿ_early, ϕᵉˣ_early, γⁱⁿᵛ_early, γᵈⁱᵛ_early, 
                                                                                                            ϕᵉⁿ_golden_age, ϕᵉˣ_golden_age, γⁱⁿᵛ_golden_age, γᵈⁱᵛ_golden_age,
                                                                                                            ϕᵉⁿ_post, ϕᵉˣ_post, γⁱⁿᵛ_post, γᵈⁱᵛ_post,
                                                                                                            exit_prob_early, invest_prob_early, entry_prob_early, exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, exit_prob_post, invest_prob_post, entry_prob_post, 
                                                                                                            df, d, m)    end

    Nsim_mean[:,i] = mean(Nsim,dims=2); Nsim_sd[:,i] = std(Nsim,dims=2); Nsim_lo[:,i] = Nsim_mean[:,i] - 1.96*Nsim_sd[:,i]; Nsim_hi[:,i] = Nsim_mean[:,i] + 1.96*Nsim_sd[:,i]
    Ksim_mean[:,i] = mean(Ksim,dims=2); Ksim_sd[:,i] = std(Ksim,dims=2); Ksim_lo[:,i] = Ksim_mean[:,i] - 1.96*Ksim_sd[:,i]; Ksim_hi[:,i] = Ksim_mean[:,i] + 1.96*Ksim_sd[:,i]
    Qsim_mean[:,i] = mean(Qsim,dims=2); Qsim_sd[:,i] = std(Qsim,dims=2); Qsim_lo[:,i] = Qsim_mean[:,i] - 1.96*Qsim_sd[:,i]; Qsim_hi[:,i] = Qsim_mean[:,i] + 1.96*Qsim_sd[:,i]
    Wsim_mean[:,i] = mean(Wsim,dims=2); Wsim_sd[:,i] = std(Wsim,dims=2); Wsim_lo[:,i] = Wsim_mean[:,i] - 1.96*Wsim_sd[:,i]; Wsim_hi[:,i] = Wsim_mean[:,i] + 1.96*Wsim_sd[:,i]
    Psim_mean[:,i] = mean(Psim,dims=2); Psim_sd[:,i] = std(Psim,dims=2); Psim_lo[:,i] = Psim_mean[:,i] - 1.96*Psim_sd[:,i]; Psim_hi[:,i] = Psim_mean[:,i] + 1.96*Psim_sd[:,i]

    PSsim = PSsim.*discount
    CSsim = CSsim.*discount
    SSsim = SSsim.*discount

    total_PS_pigouvian_tax[i] = mean(sum(PSsim,dims=1))
    total_CS_pigouvian_tax[i] = mean(sum(CSsim,dims=1))
    total_SS_pigouvian_tax[i] = mean(sum(SSsim,dims=1))
end

ΔPS_entry_tax = (total_PS_entry_tax .- total_PS_base)./1000#./abs(total_PS_base)# # rescaling to billion dollars
ΔCS_entry_tax = (total_CS_entry_tax .- total_CS_base)./1000#./abs(total_CS_base)#
ΔSS_entry_tax = (total_SS_entry_tax .- total_SS_base)./1000#./abs(total_SS_base)#

ΔPS_inv_tax = (total_PS_inv_tax .- total_PS_base)./1000#./abs(total_PS_base)# # rescaling to billion dollars
ΔCS_inv_tax = (total_CS_inv_tax .- total_CS_base)./1000#./abs(total_CS_base)#
ΔSS_inv_tax = (total_SS_inv_tax .- total_SS_base)./1000#./abs(total_SS_base)#

ΔPS_pigouvian_tax = (total_PS_pigouvian_tax .- total_PS_base)./1000#./abs(total_PS_base) # rescaling to billion dollars
ΔCS_pigouvian_tax = (total_CS_pigouvian_tax .- total_CS_base)./1000#./abs(total_CS_base)
ΔSS_pigouvian_tax = (total_SS_pigouvian_tax .- total_SS_base)./1000#./abs(total_SS_base)

scheme = ColorSchemes.Paired_10
palette = scheme.colors[1:end]

plot(collect(round(tax_lo-1,digits=2)*100:tax_width*100:round(tax_hi-1,digits=2)*100), ΔSS_entry_tax, label="",
                              xlabel="Tax rate (%)",
                              ylabel="Change in social surplus",
                              linecolor=palette[10], 
                              linewidth=2, 
                              frame=:box)
savefig(path*raw"/Output/Figures/welfare_entry_tax.pdf")

plot(collect(round(tax_lo-1,digits=2)*100:tax_width*100:round(tax_hi-1,digits=2)*100), ΔSS_inv_tax, label="",
                              xlabel="Tax rate (%)",
                              ylabel="Change in social surplus",
                              linecolor=palette[10],
                              linewidth=2, 
                              frame=:box)                              
savefig(path*raw"/Output/Figures/welfare_investment_tax.pdf")

plot(collect(round(tax_lo-1,digits=2)*100:tax_width*100:round(tax_hi-1,digits=2)*100), ΔSS_pigouvian_tax, label="",
                              xlabel="Tax rate (%)",
                              ylabel="Change in social surplus",
                              linecolor=palette[10],
                              linewidth=2, 
                              frame=:box)                              
savefig(path*raw"/Output/Figures/welfare_pigouvian_tax.pdf")


plot(collect(round(tax_lo-1,digits=2)*100:tax_width*100:round(tax_hi-1,digits=2)*100), ΔSS_pigouvian_tax, 
                              label="Pigouvian tax per whale catch",
                              xlabel="Tax rate (%)",
                              ylabel="Change in social surplus",
                              linecolor=palette[10],
                              linewidth=2, 
                              frame=:box,
                              legend=:topleft,
                              linestyle=:dash)    
plot!(collect(round(tax_lo-1,digits=2)*100:tax_width*100:round(tax_hi-1,digits=2)*100), ΔSS_entry_tax, 
                              label="Entry tax",
                              xlabel="Tax rate (%)",
                              ylabel="Change in social surplus",
                              linecolor=palette[8], 
                              linewidth=2, 
                              frame=:box)
plot!(collect(round(tax_lo-1,digits=2)*100:tax_width*100:round(tax_hi-1,digits=2)*100), ΔSS_inv_tax, 
                              label="Investment tax",
                              xlabel="Tax rate (%)",
                              ylabel="Change in welfare",
                              linecolor=palette[2],
                              linewidth=2, 
                              frame=:box,
                              linestyle=:dot)  
savefig(path*raw"/Output/Figures/welfare_taxes.pdf")




#------------------------------------------ Market outcome of entry/investment taxes/subsidies ------------------------------------------------#
# Early period
m.post = 0
d.game_start = 1804
d.game_end = 1830
θ̂_early = CSV.read(path*raw"/Output/est_params_dynamic_early_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_early[1]; ϕᵉⁿ_early = θ̂_early[1]
m.ϕ = θ̂_early[2]; ϕᵉˣ_early = θ̂_early[2]
m.Γ = θ̂_early[3:end]; γⁱⁿᵛ_early = θ̂_early[3]; γᵈⁱᵛ_early = θ̂_early[4]
m.κ *= 1.03; ϕᵉⁿ_early *= 1.03
m.ϕ *= 1.03; ϕᵉˣ_early *= 1.03
m.Γ[1] *= 1.03; γⁱⁿᵛ_early *= 1.03
m.Γ[2] *= 0.97; γᵈⁱᵛ_early *= 0.97
exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

# Golden-age
m.post = 0
d.game_start = 1831
d.game_end = 1858
θ̂_golden_age = CSV.read(path*raw"/Output/est_params_dynamic_golden_age_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_golden_age[1]; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
m.ϕ = θ̂_golden_age[2]; ϕᵉˣ_golden_age = θ̂_golden_age[2]
m.Γ = θ̂_golden_age[3:end]; γⁱⁿᵛ_golden_age = θ̂_golden_age[3]; γᵈⁱᵛ_golden_age = θ̂_golden_age[4]
m.κ *= 1.03; ϕᵉⁿ_golden_age *= 1.03
m.ϕ *= 1.03; ϕᵉˣ_golden_age *= 1.03
m.Γ[1] *= 1.03; γⁱⁿᵛ_golden_age *= 1.03
m.Γ[2] *= 0.97; γᵈⁱᵛ_golden_age *= 0.97
exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

# Post-petroleum period
m.post = 1
d.game_start = 1859
d.game_end = 1910
θ̂_post = CSV.read(path*raw"/Output/est_params_dynamic_post_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_post[1]; ϕᵉⁿ_post = θ̂_post[1]
m.ϕ = θ̂_post[2]; ϕᵉˣ_post = θ̂_post[2]
m.Γ = θ̂_post[3:end]; γⁱⁿᵛ_post = θ̂_post[3]; γᵈⁱᵛ_post = θ̂_post[4]
m.κ *= 1.03; ϕᵉⁿ_post *= 1.03
m.ϕ *= 1.03; ϕᵉˣ_post *= 1.03
m.Γ[1] *= 1.03; γⁱⁿᵛ_post *= 1.03
m.Γ[2] *= 0.97; γᵈⁱᵛ_post *= 0.97
exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE(Γ, ϕ, κ, ψ, d, m)

Random.seed!(1234)
nsim = 5_000
Nsim = zeros(Tbar+1,nsim); Ksim = zeros(Tbar+1,nsim); Qsim = zeros(Tbar+1,nsim); Wsim = zeros(Tbar+1,nsim); Psim = zeros(Tbar+1,nsim); PSsim = zeros(Tbar+1,nsim); CSsim = zeros(Tbar+1,nsim); SSsim = zeros(Tbar+1,nsim)
Threads.@threads for i = 1:nsim
    #println("sim: $i")
    Nsim[:,i], Ksim[:,i], Qsim[:,i], Wsim[:,i], Psim[:,i], PSsim[:,i], CSsim[:,i], SSsim[:,i] = simulate_data(ϕᵉⁿ_early, ϕᵉˣ_early, γⁱⁿᵛ_early, γᵈⁱᵛ_early, 
                                                                                                        ϕᵉⁿ_golden_age, ϕᵉˣ_golden_age, γⁱⁿᵛ_golden_age, γᵈⁱᵛ_golden_age,
                                                                                                        ϕᵉⁿ_post, ϕᵉˣ_post, γⁱⁿᵛ_post, γᵈⁱᵛ_post,
                                                                                                        exit_prob_early, invest_prob_early, entry_prob_early, exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, exit_prob_post, invest_prob_post, entry_prob_post, 
                                                                                                        df, d, m)
end

Nsim_mean = mean(Nsim,dims=2); Nsim_sd = std(Nsim,dims=2); Nsim_lo = Nsim_mean - 1.96*Nsim_sd; Nsim_hi = Nsim_mean + 1.96*Nsim_sd
Ksim_mean = mean(Ksim,dims=2); Ksim_sd = std(Ksim,dims=2); Ksim_lo = Ksim_mean - 1.96*Ksim_sd; Ksim_hi = Ksim_mean + 1.96*Ksim_sd
Qsim_mean = mean(Qsim,dims=2); Qsim_sd = std(Qsim,dims=2); Qsim_lo = Qsim_mean - 1.96*Qsim_sd; Qsim_hi = Qsim_mean + 1.96*Qsim_sd
Wsim_mean = mean(Wsim,dims=2); Wsim_sd = std(Wsim,dims=2); Wsim_lo = Wsim_mean - 1.96*Wsim_sd; Wsim_hi = Wsim_mean + 1.96*Wsim_sd
Psim_mean = mean(Psim,dims=2); Psim_sd = std(Psim,dims=2); Psim_lo = Psim_mean - 1.96*Psim_sd; Psim_hi = Psim_mean + 1.96*Psim_sd

# Market outcomes with entry/investment taxes

# Read the baseline simulated paths without regulation
Qsim_base = CSV.File(path*raw"/Output/Qsim.csv"; header=false) |> Tables.matrix
Psim_base = CSV.File(path*raw"/Output/Psim.csv"; header=false) |> Tables.matrix
Nsim_base = CSV.File(path*raw"/Output/Nsim.csv"; header=false) |> Tables.matrix
Ksim_base = CSV.File(path*raw"/Output/Ksim.csv"; header=false) |> Tables.matrix
Wsim_base = CSV.File(path*raw"/Output/Wsim.csv"; header=false) |> Tables.matrix
Nsim_base_mean = mean(Nsim_base,dims=2); Nsim_base_sd = std(Nsim_base,dims=2); Nsim_base_lo = Nsim_base_mean - 1.96*Nsim_base_sd; Nsim_base_hi = Nsim_base_mean + 1.96*Nsim_base_sd
Ksim_base_mean = mean(Ksim_base,dims=2); Ksim_base_sd = std(Ksim_base,dims=2); Ksim_base_lo = Ksim_base_mean - 1.96*Ksim_base_sd; Ksim_base_hi = Ksim_base_mean + 1.96*Ksim_base_sd
Qsim_base_mean = mean(Qsim_base,dims=2); Qsim_base_sd = std(Qsim_base,dims=2); Qsim_base_lo = Qsim_base_mean - 1.96*Qsim_base_sd; Qsim_base_hi = Qsim_base_mean + 1.96*Qsim_base_sd
Wsim_base_mean = mean(Wsim_base,dims=2); Wsim_base_sd = std(Wsim_base,dims=2); Wsim_base_lo = Wsim_base_mean - 1.96*Wsim_base_sd; Wsim_base_hi = Wsim_base_mean + 1.96*Wsim_base_sd
Psim_base_mean = mean(Psim_base,dims=2); Psim_base_sd = std(Psim_base,dims=2); Psim_base_lo = Psim_base_mean - 1.96*Psim_base_sd; Psim_base_hi = Psim_base_mean + 1.96*Psim_base_sd

Qfig1, Qfig2 = model_fit_visualize_counterfactuals(QData./1000, Qsim_base_mean./1000, Qsim_mean./1000, Qsim_lo./1000, Qsim_hi./1000, d)
Qfig1 = plot(Qfig1, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_data_entry_investment_taxes.pdf")
Qfig2 = plot(Qfig2, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_model_entry_investment_taxes.pdf")

Pfig1, Pfig2 = model_fit_visualize_counterfactuals(PData./1000, Psim_base_mean./1000, Psim_mean./1000, Psim_lo./1000, Psim_hi./1000, d)
Pfig1 = plot(Pfig1, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_data_entry_investment_taxes.pdf")
Pfig2 = plot(Pfig2, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_model_entry_investment_taxes.pdf")

Nfig1, Nfig2 = model_fit_visualize_counterfactuals(NData, Nsim_base_mean, Nsim_mean, Nsim_lo, Nsim_hi, d)
Nfig1 = plot(Nfig1, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_data_entry_investment_taxes.pdf")
Nfig2 = plot(Nfig2, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_model_entry_investment_taxes.pdf")

Kfig1, Kfig2 = model_fit_visualize_counterfactuals(KData./1000, Ksim_base_mean./1000, Ksim_mean./1000, Ksim_lo./1000, Ksim_hi./1000, d)
Kfig1 = plot(Kfig1, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_data_entry_investment_taxes.pdf")
Kfig2 = plot(Kfig2, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_model_entry_investment_taxes.pdf")

Wfig1, Wfig2 = model_fit_visualize_counterfactuals(WData./1_000_000, Wsim_base_mean./1_000_000, Wsim_mean./1_000_000, Wsim_lo./1_000_000, Wsim_hi./1_000_000, d)
Wfig1 = plot(Wfig1, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_data_entry_investment_taxes.pdf")
Wfig2 = plot(Wfig2, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_model_entry_investment_taxes.pdf")






























#------------------------------------- Consumer surplus -------------------------------------#
Psize = size(PData,1)


[QData[1:Psize] PData]

CS = zeros(Tbar+1)
for t = 1:1909 - d.t0 + 1
    if t < 1859 - d.t0 + 1
        cs(q,t) = exp.(- (αᵖᵒᵖ_pre/αᵖ_pre)*log.(Pop[t]) - (αᵍᵈᵖ_pre/αᵖ_pre)*log.(GDP[t]) - (αᵗ_pre/αᵖ_pre)*(d.t0 - 1800 - 1 + t) .- (α₀_pre/αᵖ_pre))*(1/(1+(1/αᵖ_pre)))*q^(1+(1/αᵖ_pre)) - PData[t]*q
    else 
        cs(q,t) = exp.(- (αᵖᵒᵖ_post/αᵖ_post)*log.(Pop[t]) - (αᵍᵈᵖ_post/αᵖ_post)*log.(GDP[t]) - (αᵗ_post/αᵖ_post)*(d.t0 - 1800 - 1 + t).- (PVec_base/αᵖ_post))*(1/(1+(1/αᵖ_post)))*q^(1+(1/αᵖ_post)) - PData[t]*q
    end
    CS[t] = cs(QData[t],t)
end

CS = CS./m.scale


#------------------------------------- Producer surplus -------------------------------------#
G = copy(x)
rename!(G, :K => :K′, :K_index => :K′_index)
df_tmp.Kprime = exp.(df_tmp.kprime)

price_df = copy(whale_prices)
rename!(price_df, :year => :tid, :p => :P)
tmp_df = leftjoin(df_tmp, price_df; on = :tid)
sort!(tmp_df, [:fid, :tid])
tmp_df.Q = exp.(tmp_df.q)

u0 = copy(tmp_df[1805 .≤ df_tmp.tid .≤ 1910, ["fid","tid","Q","K_index","A","Ω_index","exit","Kprime","P"]])
u0.K′_index = findnearest(K_df.K, u0.Kprime) # for next period fleet size
u0 = leftjoin(u0, unique(G[:,["K′_index","K′"]]); on = :K′_index)
u0 = leftjoin(u0, unique(x[:,["K_index","K"]]); on = :K_index)

u0.inv = (u0.K′ - (1 - m.δ)*u0.K)
u0 = coalesce.(u0, 0.0)
u0.inv_cost = zeros(Float64,size(u0,1))
u0.exit_cost = zeros(Float64,size(u0,1))
u0.entry_cost = zeros(Float64,size(u0,1))

u0.inv_cost[u0.inv .≥ 0 .&& u0.tid .< 1831] = (γⁱⁿᵛ_early*(u0.inv[u0.inv .≥ 0 .&& u0.tid .< 1831]).^2)/m.scale 
u0.inv_cost[u0.inv .< 0 .&& u0.tid .< 1831] = (γᵈⁱᵛ_early*(u0.inv[u0.inv .< 0 .&& u0.tid .< 1831]).^2)/m.scale 
u0.exit_cost[u0.exit .== 1 .&& u0.tid .< 1831] .= ϕᵉˣ_early
u0.entry_cost[u0.A .== 1 .&& u0.tid .< 1831] .= ϕᵉⁿ_early

u0.inv_cost[u0.inv .≥ 0 .&& 1831 .≤ u0.tid .< 1859] = (γⁱⁿᵛ_golden_age*(u0.inv[u0.inv .≥ 0 .&& 1831 .≤ u0.tid .< 1859]).^2)/m.scale 
u0.inv_cost[u0.inv .< 0 .&& 1831 .≤ u0.tid .< 1859] = (γᵈⁱᵛ_golden_age*(u0.inv[u0.inv .< 0 .&& 1831 .≤ u0.tid .< 1859]).^2)/m.scale 
u0.exit_cost[u0.exit .== 1 .&& 1831 .≤ u0.tid .< 1859] .= ϕᵉˣ_golden_age
u0.entry_cost[u0.A .== 1 .&& 1831 .≤ u0.tid .< 1859] .= ϕᵉⁿ_golden_age

u0.inv_cost[u0.inv .≥ 0 .&& 1859 .≤ u0.tid] = (γⁱⁿᵛ_post*(u0.inv[u0.inv .≥ 0 .&& 1859 .≤ u0.tid]).^2)/m.scale 
u0.inv_cost[u0.inv .< 0 .&& 1859 .≤ u0.tid] = (γᵈⁱᵛ_post*(u0.inv[u0.inv .< 0 .&& 1859 .≤ u0.tid]).^2)/m.scale 
u0.exit_cost[u0.exit .== 1 .&& 1859 .≤ u0.tid] .= ϕᵉˣ_post
u0.entry_cost[u0.A .== 1 .&& 1859 .≤ u0.tid] .= ϕᵉⁿ_post

u0.profit = u0.P.*u0.Q*.3/m.scale
u0.surplus = u0.profit - u0.inv_cost + u0.exit_cost - u0.entry_cost
PS = sum(u0.profit) - sum(u0.inv_cost) + sum(u0.exit_cost) - sum(u0.entry_cost)

PS = combine(groupby(u0[:,[:tid, :surplus]], :tid), :surplus .=> sum).surplus_sum

plot(collect(d.t0:d.t0+104-1),CS[1:104])
plot!(collect(d.t0:d.t0+104-1),PS)

sum(CS)
sum(PS)

#-------------------------------------- Social planner --------------------------------------#

annual_PS = mean(PSsim,dims=2)
annual_CS = mean(CSsim,dims=2)
annual_SS = mean(SSsim,dims=2)
total_PS = mean(sum(PSsim,dims=1))
total_CS = mean(sum(CSsim,dims=1))
total_SS = mean(sum(SSsim,dims=1))

SP_PS, SP_PS_play = findmax(sum(PSsim,dims=1))
SP_PS_play = SP_PS_play[2]
SP_CS, SP_CS_play = findmax(sum(CSsim,dims=1))
SP_CS_play = SP_CS_play[2]
SP_SS, SP_SS_play = findmax(sum(SSsim,dims=1))
SP_SS_play = SP_SS_play[2]

(total_SS-SP_SS)/SP_SS*100

Qfig1,Qfig2 = model_fit_visualize(QData./1000, Qsim[:,SP_CS_play]./1000, Qsim[:,SP_CS_play]./1000, Qsim[:,SP_CS_play]./1000, d)
Qfig = plot(Qfig, ylabel="Number (Thousands)")

Kfig1,Kfig2 = model_fit_visualize(KData./1000, Ksim[:,SP_CS_play]./1000, Ksim[:,SP_CS_play]./1000, Ksim[:,SP_CS_play]./1000, d)
Kfig = plot(Kfig, ylabel="Number (Thousands)")

Wfig1,Wfig2 = model_fit_visualize(WData./1000, Wsim[:,SP_CS_play]./1000, Wsim[:,SP_CS_play]./1000, Wsim[:,SP_CS_play]./1000, d)
Wfig = plot(Wfig, ylabel="Number (Thousands)")

Pfig,Pfig2 = model_fit_visualize(PData./1000, Psim[:,SP_CS_play]./1000, Psim[:,SP_CS_play]./1000, Psim[:,SP_SS_play]./1000, d)
Pfig = plot(Pfig, ylabel="1880 dollars (Thousands)")





#============================================================================================#
#------------------------------------ Policy experiments ------------------------------------#
#============================================================================================#

#------------------------------------ Differentiated Pigouvian tax ------------------------------------#
# Early period
m.post = 0
d.game_start = 1804
d.game_end = 1830
θ̂_early = CSV.read(path*raw"/Output/est_params_dynamic_early_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_early[1]; ϕᵉⁿ_early = θ̂_early[1]
m.ϕ = θ̂_early[2]; ϕᵉˣ_early = θ̂_early[2]
m.Γ = θ̂_early[3:end]; γⁱⁿᵛ_early = θ̂_early[3]; ; γᵈⁱᵛ_early = θ̂_early[4]
exit_prob_early, invest_prob_early, entry_prob_early, NVec, KVec, QVec, WVec, PVec = solve_NOE_diff_tax_exit(Γ, ϕ, κ, ψ, d, m)

# Golden-age
m.post = 0
d.game_start = 1831
d.game_end = 1858
θ̂_golden_age = CSV.read(path*raw"/Output/est_params_dynamic_golden_age_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_golden_age[1]; ϕᵉⁿ_golden_age = θ̂_golden_age[1]
m.ϕ = θ̂_golden_age[2]; ϕᵉˣ_golden_age = θ̂_golden_age[2]
m.Γ = θ̂_golden_age[3:end]; γⁱⁿᵛ_golden_age = θ̂_golden_age[3]; ; γᵈⁱᵛ_golden_age = θ̂_golden_age[4]
exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, NVec, KVec, QVec, WVec, PVec = solve_NOE_diff_tax_exit(Γ, ϕ, κ, ψ, d, m)

# Post-petroleum period
m.post = 1
d.game_start = 1859
d.game_end = 1910
θ̂_post = CSV.read(path*raw"/Output/est_params_dynamic_post_cost11_cap25.csv", DataFrame, header=false)[1]
m.κ = θ̂_post[1]; ϕᵉⁿ_post = θ̂_post[1]
m.ϕ = θ̂_post[2]; ϕᵉˣ_post = θ̂_post[2]
m.Γ = θ̂_post[3:end]; γⁱⁿᵛ_post = θ̂_post[3]; ; γᵈⁱᵛ_post = θ̂_post[4]
exit_prob_post, invest_prob_post, entry_prob_post, NVec, KVec, QVec, WVec, PVec = solve_NOE_diff_tax_exit(Γ, ϕ, κ, ψ, d, m)

Random.seed!(1234)
nsim = 50#_000
Nsim_diff_tax = zeros(Tbar+1,nsim); Ksim_diff_tax = zeros(Tbar+1,nsim); Qsim_diff_tax = zeros(Tbar+1,nsim); 
Wsim_diff_tax = zeros(Tbar+1,nsim); Psim_diff_tax = zeros(Tbar+1,nsim); 
PSsim_diff_tax = zeros(Tbar+1,nsim); CSsim_diff_tax = zeros(Tbar+1,nsim); SSsim_diff_tax = zeros(Tbar+1,nsim)
Threads.@threads for i = 1:nsim
    Nsim_diff_tax[:,i], Ksim_diff_tax[:,i], Qsim_diff_tax[:,i], Wsim_diff_tax[:,i], Psim_diff_tax[:,i], PSsim_diff_tax[:,i], CSsim_diff_tax[:,i], SSsim_diff_tax[:,i] = simulate_data(NVec, KVec, QVec, WVec, PVec, df, d, m)
end

Nsim_mean_diff_tax = mean(Nsim_diff_tax,dims=2); Nsim_sd_diff_tax = std(Nsim_diff_tax,dims=2); Nsim_lo_diff_tax = Nsim_mean_diff_tax - 1.96*Nsim_sd_diff_tax; Nsim_hi_diff_tax = Nsim_mean_diff_tax + 1.96*Nsim_sd_diff_tax
Ksim_mean_diff_tax = mean(Ksim_diff_tax,dims=2); Ksim_sd_diff_tax = std(Ksim_diff_tax,dims=2); Ksim_lo_diff_tax = Ksim_mean_diff_tax - 1.96*Ksim_sd_diff_tax; Ksim_hi_diff_tax = Ksim_mean_diff_tax + 1.96*Ksim_sd_diff_tax
Qsim_mean_diff_tax = mean(Qsim_diff_tax,dims=2); Qsim_sd_diff_tax = std(Qsim_diff_tax,dims=2); Qsim_lo_diff_tax = Qsim_mean_diff_tax - 1.96*Qsim_sd_diff_tax; Qsim_hi_diff_tax = Qsim_mean_diff_tax + 1.96*Qsim_sd_diff_tax
Wsim_mean_diff_tax = mean(Wsim_diff_tax,dims=2); Wsim_sd_diff_tax = std(Wsim_diff_tax,dims=2); Wsim_lo_diff_tax = Wsim_mean_diff_tax - 1.96*Wsim_sd_diff_tax; Wsim_hi_diff_tax = Wsim_mean_diff_tax + 1.96*Wsim_sd_diff_tax
Psim_mean_diff_tax = mean(Psim_diff_tax,dims=2); Psim_sd_diff_tax = std(Psim_diff_tax,dims=2); Psim_lo_diff_tax = Psim_mean_diff_tax - 1.96*Psim_sd_diff_tax; Psim_hi_diff_tax = Psim_mean_diff_tax + 1.96*Psim_sd_diff_tax

writedlm(path*raw"/Output/Nsim_diff_tax.csv", Nsim_diff_tax, ',')
writedlm(path*raw"/Output/Ksim_diff_tax.csv", Ksim_diff_tax, ',')
writedlm(path*raw"/Output/Qsim_diff_tax.csv", Qsim_diff_tax, ',')
writedlm(path*raw"/Output/Wsim_diff_tax.csv", Wsim_diff_tax, ',')
writedlm(path*raw"/Output/Psim_diff_tax.csv", Psim_diff_tax, ',')
writedlm(path*raw"/Output/PSsim_diff_tax.csv", PSsim_diff_tax, ',')
writedlm(path*raw"/Output/CSsim_diff_tax.csv", CSsim_diff_tax, ',')
writedlm(path*raw"/Output/SSsim_diff_tax.csv", SSsim_diff_tax, ',')

# Read the simulated paths from estimated model
Nsim = CSV.File(path*raw"/Output/Nsim.csv") |> Tables.matrix
Ksim = CSV.File(path*raw"/Output/Nsim.csv") |> Tables.matrix
Qsim = CSV.File(path*raw"/Output/Nsim.csv") |> Tables.matrix
Wsim = CSV.File(path*raw"/Output/Nsim.csv") |> Tables.matrix
Psim = CSV.File(path*raw"/Output/Nsim.csv") |> Tables.matrix

Nsim_mean = mean(Nsim,dims=2); Nsim_sd = std(Nsim,dims=2); Nsim_lo = Nsim_mean - 1.96*Nsim_sd; Nsim_hi = Nsim_mean + 1.96*Nsim_sd
Ksim_mean = mean(Ksim,dims=2); Ksim_sd = std(Ksim,dims=2); Ksim_lo = Ksim_mean - 1.96*Ksim_sd; Ksim_hi = Ksim_mean + 1.96*Ksim_sd
Qsim_mean = mean(Qsim,dims=2); Qsim_sd = std(Qsim,dims=2); Qsim_lo = Qsim_mean - 1.96*Qsim_sd; Qsim_hi = Qsim_mean + 1.96*Qsim_sd
Wsim_mean = mean(Wsim,dims=2); Wsim_sd = std(Wsim,dims=2); Wsim_lo = Wsim_mean - 1.96*Wsim_sd; Wsim_hi = Wsim_mean + 1.96*Wsim_sd
Psim_mean = mean(Psim,dims=2); Psim_sd = std(Psim,dims=2); Psim_lo = Psim_mean - 1.96*Psim_sd; Psim_hi = Psim_mean + 1.96*Psim_sd

# Model fit
Qfig1, Qfig2 = model_fit_visualize(QData./1000, Qsim_mean_diff_tax./1000, Qsim_lo_diff_tax./1000, Qsim_hi_diff_tax./1000, d)
Qfig1 = plot(Qfig1, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_data.pdf")
Qfig2 = plot(Qfig2, ylabel="Number (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_Q_model.pdf")

Pfig1, Pfig2 = model_fit_visualize(PData./1000, Psim_mean_diff_tax./1000, Psim_lo_diff_tax./1000, Psim_hi_diff_tax./1000, d)
Pfig1 = plot(Pfig1, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_data.pdf")
Pfig2 = plot(Pfig2, ylabel="1880 dollars (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_P_model.pdf")

Nfig1, Nfig2 = model_fit_visualize(NData, Nsim_mean_diff_tax, Nsim_lo_diff_tax, Nsim_hi_diff_tax, d)
Nfig1 = plot(Nfig1, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_data.pdf")
Nfig2 = plot(Nfig2, ylabel="Number"); savefig(path*raw"/Output/Figures/model_fit_N_model.pdf")

Kfig1, Kfig2 = model_fit_visualize(KData./1000, Ksim_mean_diff_tax./1000, Ksim_lo_diff_tax./1000, Ksim_hi_diff_tax./1000, d)
Kfig1 = plot(Kfig1, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_data.pdf")
Kfig2 = plot(Kfig2, ylabel="Tonnage (Thousands)"); savefig(path*raw"/Output/Figures/model_fit_K_model.pdf")

Wfig1, Wfig2 = model_fit_visualize(WData./1_000_000, Wsim_mean_diff_tax./1_000_000, Wsim_lo_diff_tax./1_000_000, Wsim_hi_diff_tax./1_000_000, d)
Wfig1 = plot(Wfig1, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_data.pdf")
Wfig2 = plot(Wfig2, ylabel="Number (Millions)"); savefig(path*raw"/Output/Figures/model_fit_W_model.pdf")






#============================================================================================#
#------------------------------------- Robustness check -------------------------------------#
#============================================================================================#

#------------------------------------- With ψ = sqrt(10) -------------------------------------#
m.ψ = sqrt(10) # set the logit scale parameter

# Estimate: early
m.post = 0
d.game_start = 1816
d.game_end = 1830
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_early = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_early = nfxp_early.minimizer
ll_early = nfxp_early.minimum
writedlm(path*raw"/Output/est_params_dynamic_early_10.csv",  θ̂_early, ',')
SE_early = calc_std_errors(θ̂_early, d, m)
writedlm(path*raw"/Output/est_se_early_10.csv", SE_early, ',')

# Estimate: golden-age
m.post = 0
d.game_start = 1831
d.game_end = 1858
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_golden_age = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_golden_age = nfxp_golden_age.minimizer
ll_golden_age = nfxp_golden_age.minimum
writedlm(path*raw"/Output/est_params_dynamic_golden_age_10.csv",  θ̂_golden_age, ',')
SE_golden_age = calc_std_errors(θ̂_golden_age, d, m)
writedlm(path*raw"/Output/est_se_golden_age_10.csv", SE_golden_age, ',')


# Estimate: post-petroleum
m.post = 1
d.game_start = 1859
d.game_end = 1910
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_post = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_post = nfxp_post.minimizer
ll_post = nfxp_post.minimum
writedlm(path*raw"/Output/est_params_dynamic_post_10.csv",  θ̂_post, ',')
SE_post = calc_std_errors(θ̂_post, d, m)
writedlm(path*raw"/Output/est_se_post_10.csv", SE_post, ',')

writedlm(path*raw"/Output/est_dynamic_costs_10.csv", [θ̂_early θ̂_golden_age θ̂_post; 
                                                   SE_early SE_golden_age SE_post;
                                                   ll_early ll_golden_age ll_post], ',')

#------------------------------------- With ψ = 1/sqrt(10) -------------------------------------#
m.ψ = 1/sqrt(10) # set the logit scale parameter

# Estimate: early
m.post = 0
d.game_start = 1816
d.game_end = 1830
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_early = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_early = nfxp_early.minimizer
ll_early = nfxp_early.minimum
writedlm(path*raw"/Output/est_params_dynamic_early_1_10.csv",  θ̂_early, ',')
SE_early = calc_std_errors(θ̂_early, d, m)
writedlm(path*raw"/Output/est_se_early_1_10.csv", SE_early, ',')

# Estimate: golden-age
m.post = 0
d.game_start = 1831
d.game_end = 1858
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_golden_age = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_golden_age = nfxp_golden_age.minimizer
ll_golden_age = nfxp_golden_age.minimum
writedlm(path*raw"/Output/est_params_dynamic_golden_age_1_10.csv",  θ̂_golden_age, ',')
SE_golden_age = calc_std_errors(θ̂_golden_age, d, m)
writedlm(path*raw"/Output/est_se_golden_age_1_10.csv", SE_golden_age, ',')


# Estimate: post-petroleum
m.post = 1
d.game_start = 1859
d.game_end = 1910
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_post = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_post = nfxp_post.minimizer
ll_post = nfxp_post.minimum
writedlm(path*raw"/Output/est_params_dynamic_post_1_10.csv",  θ̂_post, ',')
SE_post = calc_std_errors(θ̂_post, d, m)
writedlm(path*raw"/Output/est_se_post_1_10.csv", SE_post, ',')

writedlm(path*raw"/Output/est_dynamic_costs_1_10.csv", [θ̂_early θ̂_golden_age θ̂_post; 
                                                   SE_early SE_golden_age SE_post;
                                                   ll_early ll_golden_age ll_post], ',')


#------------------------------------- Estimate ψ together -------------------------------------#

# Estimate: early
m.post = 0
d.game_start = 1816
d.game_end = 1830
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_early = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_early = nfxp_early.minimizer
ll_early = nfxp_early.minimum
writedlm(path*raw"/Output/est_params_dynamic_early_with_scale.csv",  θ̂_early, ',')
SE_early = calc_std_errors(θ̂_early, d, m)
writedlm(path*raw"/Output/est_se_early_with_scale.csv", SE_early, ',')

# Estimate: golden-age
m.post = 0
d.game_start = 1831
d.game_end = 1858
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_golden_age = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_golden_age = nfxp_golden_age.minimizer
ll_golden_age = nfxp_golden_age.minimum
writedlm(path*raw"/Output/est_params_dynamic_golden_age_with_scale.csv",  θ̂_golden_age, ',')
SE_golden_age = calc_std_errors(θ̂_golden_age, d, m)
writedlm(path*raw"/Output/est_se_golden_age_with_scale.csv", SE_golden_age, ',')


# Estimate: post-petroleum
m.post = 1
d.game_start = 1859
d.game_end = 1910
time, state, decision, decision_entry, decision_pe_quit = build_data(df, d, m)
m.cost_struct = 4

@time nfxp_post = optimize_log_likelihood(time, state, decision, decision_entry, decision_pe_quit, d,m)
θ̂_post = nfxp_post.minimizer
ll_post = nfxp_post.minimum
writedlm(path*raw"/Output/est_params_dynamic_post_with_scale.csv",  θ̂_post, ',')
SE_post = calc_std_errors(θ̂_post, d, m)
writedlm(path*raw"/Output/est_se_post_with_scale.csv", SE_post, ',')

writedlm(path*raw"/Output/est_dynamic_costs_with_scale.csv", [θ̂_early θ̂_golden_age θ̂_post; 
                                                                SE_early SE_golden_age SE_post;
                                                                ll_early ll_golden_age ll_post], ',')
#============================================================================================#
#--------------------------------- End of robustness check ----------------------------------#
#============================================================================================#
