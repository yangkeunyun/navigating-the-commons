function simulate_data(ϕᵉⁿ_early, ϕᵉˣ_early, γⁱⁿᵛ_early, γᵈⁱᵛ_early, 
    ϕᵉⁿ_golden_age, ϕᵉˣ_golden_age, γⁱⁿᵛ_golden_age, γᵈⁱᵛ_golden_age,
    ϕᵉⁿ_post, ϕᵉˣ_post, γⁱⁿᵛ_post, γᵈⁱᵛ_post,
    exit_prob_early, invest_prob_early, entry_prob_early, exit_prob_golden_age, invest_prob_golden_age, entry_prob_golden_age, exit_prob_post, invest_prob_post, entry_prob_post, 
    df, d, m)
PVec = zeros(eltype(m.κ),Tbar+1)            # price vector
WVec = zeros(eltype(m.κ),Tbar+1)            # whale population in each period
KVec = zeros(eltype(m.κ),Tbar+1)            # aggregate capacity in each period
QVec = zeros(eltype(m.κ),Tbar+1)            # aggregate whale catch in each period
QᶠVec = zeros(eltype(m.κ),Tbar+1)            # aggregate whale catch outside of US whaling
NVec = zeros(eltype(m.κ),Tbar+1)            # aggregate number of incumbents in each period
WVec[[1]] = W₀
# Simulation of the game starts from here
u1 = df_tmp[df_tmp.tid .== 1804, :]
#df[df.tid .== 1804, [:K,:Q]]
u1.time = u1.tid .- d.t0 .+ 1 

state_capacity = Int.(findnearest(K_df.K, u1.K))
state_age = Int.(findnearest(A_df.A, u1.A))
state_productivity = Int.(findnearest(Ω_df.Ω, u1.Ω))
if m.state_space == "KAΩ"
state = [state_capacity state_age state_productivity]
state = ((state[:,1] .- 1)*m.A_size*m.Ω_size .+ (state[:,2].-1)*m.Ω_size .+ state[:,3])
elseif m.state_space == "KΩ"
state = [state_capacity state_productivity]
state = ((state[:,1] .- 1)*m.Ω_size .+ state[:,2])    
end
decision_capacity_prime = Int.(findnearest(K_df.K, exp.(u1.kprime)))
decision_investment = decision_capacity_prime
decision_stay = 1 .- u1.exit
decision = decision_investment.*decision_stay

time = u1.time

u1.time = time; u1.state = state
u1 = u1[:,[:fid, :time, :state]]
u1 = leftjoin(u1, x[:,[:K,:K_index,:Ω,:Ω_index,:state]]; on = :state)
colnames = ["fid", "time", "K_index", "Ω_index"]
sim_data = DataFrame([name => Int[] for name in colnames])
Q = zeros(eltype(m.κ),m.x_size)
PVec = zeros(eltype(m.κ),Tbar+1)

G = copy(x)
rename!(G, :K_index => :decision)

PS = zeros(eltype(m.κ),Tbar+1)
PVec_base = log(QData[1859 - d.t0 + 1]) .- αᵖ_post*log(PData[1859 - d.t0 + 1]) - αᵖᵒᵖ_post*log(Pop[1859 - d.t0 + 1]) - αᵍᵈᵖ_post*log.(GDP[1859 - d.t0 + 1]) - αᵖᵉᵗ_post*log(Pet[1859 - d.t0 + 1]) - αᵗ_post*(d.t0 - 1800 - 1 + 1859 - d.t0 + 1) 
CS = zeros(Tbar+1)

for t = 1804 - d.t0 + 1:1910 - d.t0 + 1
#println("time: $t")
# Data for the initial period of the game
u0 = copy(u1)

NVec[[t]] .= size(u0,1)
KVec[[t]] .= sum(u0.K)
K_missing = sum(df[df.tid .== t + d.t0 - 1 .&& ismissing.(df.fid), [:K]].K)
KVec[[t]] .+= K_missing
if KVec[t] < 562 # minimum aggregate fleet size in Tower (1907)
KVec[t] = 562
end

calc_production!(Q, KVec, WVec, d, t)
x.Q = Q
u0 = leftjoin(u0,unique(x[:,[:K_index,:Ω_index,:Q]]); on = [:K_index,:Ω_index])
QVec[[t]] .= sum(u0.Q)
Q_missing = sum(df[df.tid .== t + d.t0 - 1 .&& ismissing.(df.fid), [:Q]].Q)
QVec[[t]] .+= Q_missing
if isequal(QVec[t], NaN)
QVec[t] = QVec[t-1] 
end    
if QVec[t] < 60 # minimum aggregate whale catch calculated from data of Tower (1907)
QVec[t] = 60
end

WVec[[t+1]] .= WVec[t] .+ r_max*WVec[t]*(1 .- (WVec[t]./W₀).^z) .- QVec[t]*1.71
if WVec[t+1] < 0
QVec[[t]] .= max((WVec[t] + r_max*WVec[t]*(1 - (WVec[t]/W₀[1])^z))/1.71 - WMin, .1)
WVec[[t+1]] .= WVec[t]
end    

if t == 1859 - d.t0 + 1 
PVec_base = log(QData[t]) .- αᵖ_post*log(PData[t]) - αᵖᵒᵖ_post*log(Pop[t]) - αᵍᵈᵖ_post*log.(GDP[t]) - αᵖᵉᵗ_post*log(Pet[t]) - αᵗ_post*(d.t0 - 1800 - 1 + t) 
end

if t < 1859 - d.t0 + 1
PVec[[t]] .= exp.((1/αᵖ_pre).*log.(QVec[t]) - (αᵖᵒᵖ_pre/αᵖ_pre)*log.(Pop[t]) 
            - (αᵍᵈᵖ_pre/αᵖ_pre)*log.(GDP[t]) - (αᵗ_pre/αᵖ_pre)*(d.t0 - 1800 - 1 + t)
            .- (α₀_pre/αᵖ_pre))
else 
PVec[[t]] .= exp.((1/αᵖ_post).*log.(QVec[t]) - (αᵖᵒᵖ_post/αᵖ_post)*log.(Pop[t]) 
            - (αᵍᵈᵖ_post/αᵖ_post)*log.(GDP[t]) - (αᵗ_post/αᵖ_post)*(d.t0 - 1800 - 1 + t) 
            - (αᵖᵉᵗ_post/αᵖ_post)*Post1859[t].*log.(Pet[t]) .- (PVec_base/αᵖ_post))
end

time = u0.time
if m.state_space == "KAΩ"
state = [K_index A_index Ω_index]
state = ((state[:,1] .- 1)*m.A_size*m.Ω_size .+ (state[:,2].-1)*m.Ω_size .+ state[:,3])
elseif m.state_space == "KΩ"
state = [u0.K_index u0.Ω_index]
state = ((state[:,1] .- 1)*m.Ω_size .+ state[:,2])    
end

# Get incumbent decisions based on conditional choice probabilities
if t < 1831 - d.t0 + 1
Pr_exit = [exit_prob_early[t,i] for (t,i) in zip(time, state)]
Pr_inv = [invest_prob_early[t,i,:] for (t,i) in zip(time, state)]
elseif 1831 - d.t0 + 1 ≤ t < 1859 - d.t0 + 1
Pr_exit = [exit_prob_golden_age[t,i] for (t,i) in zip(time, state)]
Pr_inv = [invest_prob_golden_age[t,i,:] for (t,i) in zip(time, state)]
else
Pr_exit = [exit_prob_post[t,i] for (t,i) in zip(time, state)]
Pr_inv = [invest_prob_post[t,i,:] for (t,i) in zip(time, state)]
end
σ = [[Pr_exit[i]; Pr_inv[i]] for i in 1:length(Pr_exit)]
u0.decision = rand.(Categorical.(σ)) .- 1

# Get productivity transition
Pr_prod = [row |> x -> Vector(x) for row in eachrow(ω_trans[u0.Ω_index[u0.decision .> 0], :])]

next_prod = rand.(Categorical.(Pr_prod))

# New entrants
Nᵖᵉ = max.(m.N̄ .- size(u0,1),0)
if t < 1831 - d.t0 + 1
Nᵉ = rand.(Bernoulli(entry_prob_early[t]), Nᵖᵉ)
elseif 1831 - d.t0 + 1 ≤ t < 1859 - d.t0 + 1
Nᵉ = rand.(Bernoulli(entry_prob_golden_age[t]), Nᵖᵉ)
else
Nᵉ = rand.(Bernoulli(entry_prob_post[t]), Nᵖᵉ)
end

tmp = ones(Int64, sum(Nᵉ), 3) .* [t; Kᵉ; Ωᵉ]'

# Get costs
u0 = leftjoin(u0, unique(G[:,["K","decision"]]); on = :decision, makeunique = true)
u0.inv = (u0.K_1 - (1 - m.δ)*u0.K)
u0 = coalesce.(u0, 0.0)
u0.inv_cost = zeros(Float64,size(u0,1))
u0.exit_cost = zeros(Float64,size(u0,1))
if t < 1831 - d.t0 + 1
u0.inv_cost[u0.inv .≥ 0] = γⁱⁿᵛ_early*((u0.inv[u0.inv .≥ 0]).^2)/m.scale 
u0.inv_cost[u0.inv .< 0] = γᵈⁱᵛ_early*((u0.inv[u0.inv .< 0]).^2)/m.scale 
u0.exit_cost[u0.decision .== 0] .= ϕᵉˣ_early
# Get profits
u0.profit = PVec[t]*u0.Q*.3/m.scale
if size(u0,1) == 0
PS[t] = 0
else    
PS[t] = sum(u0.profit) - sum(u0.inv_cost) + sum(u0.exit_cost) - sum(Nᵉ)*ϕᵉⁿ_early
end
elseif 1831 - d.t0 + 1 ≤ t < 1859 - d.t0 + 1
u0.inv_cost[u0.inv .≥ 0] = γⁱⁿᵛ_golden_age*((u0.inv[u0.inv .≥ 0]).^2)/m.scale 
u0.inv_cost[u0.inv .< 0] = γᵈⁱᵛ_golden_age*((u0.inv[u0.inv .< 0]).^2)/m.scale 
u0.exit_cost[u0.decision .== 0] .= ϕᵉˣ_golden_age
# Get profits
u0.profit = PVec[t]*u0.Q*.3/m.scale
if size(u0,1) == 0
PS[t] = 0
else
PS[t] = sum(u0.profit) - sum(u0.inv_cost) + sum(u0.exit_cost) - sum(Nᵉ)*ϕᵉⁿ_golden_age
end
else
u0.inv_cost[u0.inv .≥ 0] = γⁱⁿᵛ_post*((u0.inv[u0.inv .≥ 0]).^2)/m.scale 
u0.inv_cost[u0.inv .< 0] = γᵈⁱᵛ_post*((u0.inv[u0.inv .< 0]).^2)/m.scale 
u0.exit_cost[u0.decision .== 0] .= ϕᵉˣ_post
# Get profits
u0.profit = PVec[t]*u0.Q*.3/m.scale
if size(u0,1) == 0
PS[t] = 0
else
PS[t] = sum(u0.profit) - sum(u0.inv_cost) + sum(u0.exit_cost) - sum(Nᵉ)*ϕᵉⁿ_post
end
end

# Assuming unique(u0.fid) contains the existing integers
if t == 1
existing_fid = unique(u0.fid)    
else
existing_fid = unique(sim_data.fid)
end

# Calculate the total number of new fid you want to generate
total_new_fid = sum(Nᵉ)

# Initialize an empty array to store the new unique fid
new_unique_fid = Int[]

# Generate new unique fid
while length(new_unique_fid) < total_new_fid    
new_fid = rand(1:maximum(existing_fid) + total_new_fid)
if !(new_fid in existing_fid || new_fid in new_unique_fid)
push!(new_unique_fid, new_fid)
end
end

# Ensure you have exactly the desired number of new unique fid
new_unique_fid = new_unique_fid[1:total_new_fid]
#unique([existing_fid;new_unique_fid])
tmp = [new_unique_fid tmp]

# Next period
u1 = u0[u0.decision .> 0, [:fid, :time, :decision]]
u1.Ω_index = next_prod
rename!(u1, :decision => :K_index)
# Convert the matrix 'tmp' into a DataFrame
tmp = DataFrame(tmp, ["fid","time","K_index","Ω_index"])
# Append 'tmp' to 'u0' vertically
u1 = vcat(u1, tmp)
u1.time .+= 1

u1 = leftjoin(u1,unique(x[:,[:K,:K_index]]); on = :K_index)
u1 = leftjoin(u1,unique(x[:,[:Ω,:Ω_index]]); on = :Ω_index)

if t == 1804 - d.t0 + 1
u0 = u0[:,[:fid,:time,:K_index,:Ω_index,:K,:Ω]]
sim_data = vcat(u0,u1)
else
sim_data = [sim_data; u1]
end
sort!(sim_data,[:fid,:time])

u1 = sim_data[sim_data.time .== t + 1,:] # Simulated data for the next period

#println(t)

#CS[t] = cs(PVec[t],QVec[t],t)
if t < 1859 - d.t0 + 1
CS[t] = exp.(- (αᵖᵒᵖ_pre/αᵖ_pre)*log.(Pop[t]) - (αᵍᵈᵖ_pre/αᵖ_pre)*log.(GDP[t]) - (αᵗ_pre/αᵖ_pre)*(d.t0 - 1800 - 1 + t) .- (α₀_pre/αᵖ_pre))*(1/(1+(1/αᵖ_pre)))*QVec[t]^(1+(1/αᵖ_pre)) - PVec[t]*QVec[t]
else 
CS[t] = exp.(- (αᵖᵒᵖ_post/αᵖ_post)*log.(Pop[t]) - (αᵍᵈᵖ_post/αᵖ_post)*log.(GDP[t]) - (αᵗ_post/αᵖ_post)*(d.t0 - 1800 - 1 + t).- (PVec_base/αᵖ_post))*(1/(1+(1/αᵖ_post)))*QVec[t]^(1+(1/αᵖ_post)) - PVec[t]*QVec[t]
end
end

CS = CS./m.scale
SS = CS + PS

return NVec, KVec, QVec, WVec, PVec, PS, CS, SS
end


function model_fit_visualize(data, model, model_lo, model_hi, d)
scheme = ColorSchemes.Purples_3
palette = scheme.colors[1:end]
#data = QData; model = Qsim_mean; model_lo = Qsim_lo; model_hi = Qsim_hi
#data = PData; model = Psim_mean; model_lo = Psim_lo; model_hi = Psim_hi
if data == PData./1000
Psize = size(data[1804-d.t0+1:end],1)
fig_data = plot(
collect(1804:1804+Psize-1),
data[1804-d.t0+1:end], 
frame=:box, 
xlabel="Year", 
#title="Data", 
linecolor="black", 
linewidth=2, 
xticks=1800:15:1920, 
label="Data",
legend=:topleft,
size=(700,350), dpi=1000, margin=3mm)    
else
fig_data = plot(collect(1804:1910-1),
data[1804-d.t0+1:1910-d.t0], 
frame=:box, 
xlabel="Year", 
#title="Data", 
titlefont=font(12),
linecolor="black", 
linewidth=2, 
xticks=1800:15:1920,
label="Data",
size=(700,350), dpi=1000, margin=3mm)
end
fig_model = plot(fig_data,
collect(1804:1910-1),
model[1804-d.t0+1:1910-d.t0], 
frame=:box, 
xlabel="Year", 
#title="Model", 
titlefont=font(12),
linecolor=palette[3], 
linewidth=2, 
linestyle=:dash, 
xticks=1800:15:1920,
label="Model",
size=(700,350), dpi=1000, margin=3mm)
plot!(fig_model, 
collect(1804:1910-1),
model_lo[1804-d.t0+1:1910-d.t0], fillrange = model_hi[1804-d.t0+1:1910-d.t0], fillalpha = 0.30, c=palette[2],
frame=:box, 
linestyle=:dash, 
linecolor=palette[2],
label="",
size=(700,350), dpi=1000, margin=3mm)
plot!(fig_model,
collect(1804:1910-1),
model_hi[1804-d.t0+1:1910-d.t0], 
frame=:box, 
linestyle=:dash, 
linecolor=palette[2],
label="",
size=(700,350), dpi=1000, margin=3mm)

# Get the y-limits of the data plot
if data == PData./1000
data_ylims_min, data_ylims_max = extrema(data[1804-d.t0+1:end])
else
data_ylims_min, data_ylims_max = extrema(data[1804-d.t0+1:1910-d.t0])
end
model_lo_ylims_min, model_lo_ylims_max = extrema(model_lo[1804-d.t0+1:1910-d.t0])
model_hi_ylims_min, model_hi_ylims_max = extrema(model_hi[1804-d.t0+1:1910-d.t0])
ylims_min = min(data_ylims_min,model_lo_ylims_min,model_hi_ylims_min)
ylims_max = max(data_ylims_max,model_lo_ylims_max,model_hi_ylims_max)
# Set the same y-axis limits for both plots
ylims!(fig_data, ylims_min, ylims_max)
ylims!(fig_model, ylims_min, ylims_max)
#fig = plot(fig_data, fig_model, 
#        legend = false, 
#        yrotation = 90, 
#        layout=grid(1,2), 
#        size=(800,300), 
#        dpi=1000, margin=3mm, topmargin=0.5mm)#
return fig_data, fig_model
end

function model_fit_visualize_counterfactuals(data, model_base, model, model_lo, model_hi, d)
scheme = ColorSchemes.Purples_3
palette = scheme.colors[1:end]
#data = QData; model = Qsim_mean; model_lo = Qsim_lo; model_hi = Qsim_hi
#data = PData; model = Psim_mean; model_lo = Psim_lo; model_hi = Psim_hi
if data == PData./1000
Psize = size(data[1804-d.t0+1:end],1)
fig_data = plot(
collect(1804:1804+Psize-1),
data[1804-d.t0+1:end], 
frame=:box, 
xlabel="Year", 
#title="Data", 
linecolor="black", 
linewidth=2, 
xticks=1800:15:1920, 
label="Data",
legend=:topleft,
size=(700,350), dpi=1000, margin=3mm)    
else
fig_data = plot(collect(1804:1910-1),
data[1804-d.t0+1:1910-d.t0], 
frame=:box, 
xlabel="Year", 
#title="Data", 
titlefont=font(12),
linecolor="black", 
linewidth=2, 
xticks=1800:15:1920,
label="Data",
size=(700,350), dpi=1000, margin=3mm)
end
fig_model = plot(fig_data,
collect(1804:1910-1),
model_base[1804-d.t0+1:1910-d.t0], 
frame=:box, 
xlabel="Year", 
#title="Model", 
titlefont=font(12),
linecolor=palette[3], 
linewidth=2, 
linestyle=:dash, 
xticks=1800:15:1920,
label="Model without tax/subsidy",
size=(700,350), dpi=1000, margin=3mm)
plot!(fig_model,
collect(1804:1910-1),
model[1804-d.t0+1:1910-d.t0], 
frame=:box, 
xlabel="Year", 
#title="Model", 
titlefont=font(12),
linecolor=:orange, 
linewidth=2, 
linestyle=:dot, 
xticks=1800:15:1920,
label="Model with taxes/subsidies",
size=(700,350), dpi=1000, margin=3mm)
#plot!(fig_model, 
#        collect(1804:1910-1),
#        model_lo[1804-d.t0+1:1910-d.t0], fillrange = model_hi[1804-d.t0+1:1910-d.t0], fillalpha = 0.30, c=palette[2],
#        frame=:box, 
#        linestyle=:dash, 
#        linecolor=palette[2],
#        label="",
#        size=(700,350), dpi=1000, margin=3mm)
#plot!(fig_model,
#        collect(1804:1910-1),
#        model_hi[1804-d.t0+1:1910-d.t0], 
#        frame=:box, 
#        linestyle=:dash, 
#        linecolor=palette[2],
#        label="",
#        size=(700,350), dpi=1000, margin=3mm)

# Get the y-limits of the data plot
if data == PData./1000
data_ylims_min, data_ylims_max = extrema(data[1804-d.t0+1:end])
else
data_ylims_min, data_ylims_max = extrema(data[1804-d.t0+1:1910-d.t0])
end
model_lo_ylims_min, model_lo_ylims_max = extrema(model_lo[1804-d.t0+1:1910-d.t0])
model_hi_ylims_min, model_hi_ylims_max = extrema(model_hi[1804-d.t0+1:1910-d.t0])
ylims_min = min(data_ylims_min,model_lo_ylims_min,model_hi_ylims_min)
ylims_max = max(data_ylims_max,model_lo_ylims_max,model_hi_ylims_max)
# Set the same y-axis limits for both plots
ylims!(fig_data, ylims_min, ylims_max)
ylims!(fig_model, ylims_min, ylims_max)
#fig = plot(fig_data, fig_model, 
#        legend = false, 
#        yrotation = 90, 
#        layout=grid(1,2), 
#        size=(800,300), 
#        dpi=1000, margin=3mm, topmargin=0.5mm)#
return fig_data, fig_model
end