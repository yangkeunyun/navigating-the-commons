

function generate_smoothed_transition_matrix(
    data_matrix::Matrix{Float64},
    smoothing_factor::Float64
    )::Matrix{Float64}

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


function get_state_space!(
    m::DG_model,
    df::DataFrame
    )::Tuple{DataFrame, DataFrame, DataFrame, DataFrame, Matrix{Float64}}

    println("Constructing state space")

    # Productivity space

    ω_min = percentile(df.ω,4); ω_max = percentile(df.ω,99)
    ω_space = collect(LinRange(ω_min,ω_max,m.Ω_size))
    Ω_df = DataFrame(Ω=ω_space, Ω_index=collect(1:m.Ω_size))

    Ω_df.Ω = exp.(Ω_df.Ω)

    # Productivity transition matrix directly from data
    df.Ω = exp.(df.ω)
    df.Ω_index = findnearest(Ω_df.Ω, df.Ω)
    ω_df = copy(df)#density(df.Ω_index) #freqtable(df.Ω_index)
    ω_transition_df = combine(groupby(ω_df, :fid), :Ω_index => Base.Fix2(lead, 1) => :Ω′_index)
    ω_transition_df.Ω_index = ω_df.Ω_index
    ω_transition_df2 = combine(groupby(ω_df, :fid), :ω => Base.Fix2(lead, 1) => :ω′)
    ω_transition_df.ω = ω_df.ω
    ω_transition_df.ω′ = ω_transition_df2.ω′
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
    ω_trans = Tables.matrix(reshape_trans[:,2:end])

    # Using Laplace smoothing with a smoothing factor of 0.1
    smoothing_factor = 0.05
    ω_trans = generate_smoothed_transition_matrix(ω_trans, smoothing_factor)

    # Capacity space
    df.K = exp.(df.k)
    Ton = df[!,:K]
    K_max = maximum(Ton)+Statistics.std(Ton)
    K_space = collect(100:m.K_size:K_max)
    m.K_size = length(K_space)
    K_df = DataFrame(K=K_space, K_index=collect(1:m.K_size))

    # Age space
    Age = df.A
    Age = skipmissing(Age)
    A_lo = minimum(Age); A_hi = m.A_size
    A_space = collect(A_lo:A_hi) # collect(A_lo:A_hi)
    A_df = DataFrame(A=A_space, A_index=collect(1:m.A_size))

    if m.state_space == "KAΩ"
        x = crossjoin(K_df,A_df,Ω_df)
        x_index = unique(x[:,[:K_index, :A_index, :Ω_index]])
        x.state = collect(1:size(x_index,1))
        m.x_size = size(x,1)
    elseif m.state_space == "KΩ"
        x = crossjoin(K_df,Ω_df)
        x_index = unique(x[:,[:K_index, :Ω_index]])
        x.state = collect(1:size(x_index,1))
        m.x_size = size(x,1)
    end

    m.social_surplus_stay = zeros(eltype(m.Γ),m.x_size,m.K_size)
    m.pr_capacity_prime = zeros(m.x_size,m.K_size)
    m.EVᵐᵃˣ = zeros(eltype(m.Γ),m.x_size)
    m.Ω_size_collect_tr = collect(1:m.Ω_size)'
    m.Js_vec = [zeros(Int, m.x_size*m.Ω_size) for _ in 1:m.K_size]
    m.Π_Ω_vec = [sparse(1:m.x_size, 1:m.x_size, 0.) for _ in 1:m.K_size]
    m.cost_vec = [zeros(m.x_size) for _ in 1:m.K_size]

    return x, Ω_df, K_df, A_df, ω_trans
    
end

