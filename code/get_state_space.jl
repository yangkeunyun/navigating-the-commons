println("Construct state space")

df_tmp = copy(df)
# Productivity space
rename!(df_tmp, :omega=>:ω)
dropmissing!(df_tmp, :ω)
ω_min = percentile(df_tmp.ω,4); ω_max = percentile(df_tmp.ω,99)
#ω_min = percentile(exp.(df_tmp.ω),.5); ω_max = percentile(exp.(df_tmp.ω),99.5)
#ω_min = minimum(df_tmp.ω); ω_max = maximum(df_tmp.ω)
#density(df_tmp.ω)
ω_width = (ω_max - ω_min)/(m.Ω_size - 1)
ω_space = collect(LinRange(ω_min,ω_max,m.Ω_size))
Ω_df = DataFrame(Ω=ω_space, Ω_index=collect(1:m.Ω_size))

#df_tmp.Ω_index = 1
#for p = 1:m.Ω_size-2
#    df_tmp.Ω_index[ω_space[p] .≤ df_tmp.ω .< ω_space[p+1]] .= p+1
#end
#df_tmp.Ω_index[ω_space[14] .≤ df_tmp.ω] .= 15
#density(df_tmp.Ω_index)

Ω_df.Ω = exp.(Ω_df.Ω)

# Productivity transition matrix directly from data
df_tmp.Ω = exp.(df_tmp.ω)
df_tmp.Ω_index = findnearest(Ω_df.Ω, df_tmp.Ω)
ω_df = copy(df_tmp)#density(df.Ω_index) #freqtable(df.Ω_index)
#df.Ω_index = levelcode.(cut(df.Ω_index, m.Ω_size))
#density(dropmissing(df, [:Ω_index]).Ω_index)
ω_transition_df = combine(groupby(ω_df, :fid), :Ω_index => Base.Fix2(lead, 1) => :Ω′_index)
ω_transition_df.Ω_index = ω_df.Ω_index
ω_transition_df2 = combine(groupby(ω_df, :fid), :ω => Base.Fix2(lead, 1) => :ω′)
ω_transition_df.ω = ω_df.ω
ω_transition_df.ω′ = ω_transition_df2.ω′
#ω_transition_df = CSV.read(path*raw"/Data/Intermediate/omega_transition.csv", DataFrame) 
dropmissing!(ω_transition_df, [:Ω_index, :Ω′_index])

ω_transition_df = transform(groupby(ω_transition_df, [:Ω_index, :Ω′_index]), nrow => :cat_countmar) # get the number of occurrences for each transition
ω_transition_df = transform(groupby(ω_transition_df, :Ω_index), nrow => :cat_count) # get the number of occurrences for each state

ω_transition_df[:,:prob] = ω_transition_df[:,:cat_countmar]./ω_transition_df[:,:cat_count] # get the transition probability 
trans = unique(ω_transition_df[:,[:Ω_index,:Ω′_index,:prob]]) # remove unnecessary rows and columns  
sort!(trans, :Ω′_index) # Caution! It is needed in Linux
reshape_trans = unstack(trans,:Ω′_index,:prob)
sort!(reshape_trans, :Ω_index)
# Transform into a matrix format
reshape_trans = coalesce.(reshape_trans, 0.0)
#ω_trans = convert(Matrix{Float64},reshape_trans[:,2:end]) # Finally, convert the dataframe to a matrix
ω_trans = Tables.matrix(reshape_trans[:,2:end])

function generate_smoothed_transition_matrix(data_matrix, smoothing_factor)
    grid_size = size(data_matrix, 1)
    transition_matrix = zeros(Float64, grid_size, grid_size)

    for i in 1:grid_size
        adjusted_probabilities = zeros(Float64, grid_size)
        for j in 1:grid_size
            # Adjust the smoothing factor based on the direction of the transition
            direction_factor = 1.0 / max(1, abs(i - j))

            # Laplace smoothing with adjusted factor
            adjusted_probabilities[j] = (data_matrix[i, j] + smoothing_factor * direction_factor)
        end

        # Normalize the adjusted probabilities
        transition_matrix[i, :] = adjusted_probabilities / sum(adjusted_probabilities)
    end

    return transition_matrix
end

# Using Laplace smoothing with a smoothing factor of 0.1
smoothing_factor = 0.05
ω_trans = generate_smoothed_transition_matrix(ω_trans, smoothing_factor)

#VMat[d.game_start - d.t0 + 2,:]
#=
ω_lag = combine(groupby(df, :fid), :omega => Base.Fix2(lag, 1) => :L1omega)
df_tmp2 = hcat(df, ω_lag, makeunique=true)
dropmissing(df_tmp2, :fid)
reg_omega = reg(df_tmp2, @formula(omega ~ L1omega))
#reg_omega = reg(df_tmp2, @formula(omega ~ L1omega + fe(fid)))
θ_λ = coef(reg_omega)[2]

μ = mean(df_tmp.ω)
σ² = (std(df_tmp.ω))^2
znum = m.Ω_size
zdev = 3

ω_trans = tauchen(μ,σ²,θ_λ,znum,zdev)[2]
=#
#sum(ω_trans,dims=1)

using Statistics
# Capacity space
df_tmp.K = exp.(df_tmp.k)
Ton = df_tmp[!,:K]
K_min = minimum(Ton); K_max = maximum(Ton)+Statistics.std(Ton)
#K_min = percentile(Ton,1); K_max = percentile(Ton,99)
#K_min = 80; K_max = 8_000
K_width = (K_max - K_min)/(m.K_size - 1)
#K_space = collect(LinRange(K_min,K_max,m.K_size))
#K_space = collect(minimum(df_tmp.K)-10:150:maximum(df_tmp.K)+10)
K_space = collect(100:m.K_size:K_max)
m.K_size = length(K_space)
K_df = DataFrame(K=K_space, K_index=collect(1:m.K_size))

# Age space
#df.A = exp.(df.A)
Age = df.A
Age = skipmissing(Age)
A_lo = minimum(Age); A_hi = m.A_size
A_space = collect(A_lo:A_hi) # collect(A_lo:A_hi)
A_df = DataFrame(A=A_space, A_index=collect(1:m.A_size))

# Individual state space, x = (Ω, K, A)
if m.state_space == "KAΩ"
    x = crossjoin(K_df,A_df,Ω_df)
    x_index = unique(x[:,[:K_index, :A_index, :Ω_index]])
    x.state = collect(1:size(x_index,1))
    m.x_size = size(x,1)
    x_space = collect(1:m.x_size)
    K = x.K
    A = x.A
    Ω = x.Ω
elseif m.state_space == "KΩ"
    x = crossjoin(K_df,Ω_df)
    x_index = unique(x[:,[:K_index, :Ω_index]])
    x.state = collect(1:size(x_index,1))
    m.x_size = size(x,1)
    x_space = collect(1:m.x_size)
    K = x.K
    Ω = x.Ω
end

