function findnearest(A::AbstractArray,values)
    index = zeros(size(values))
    for (i,i_val) in enumerate(values)
        index[i] = findmin(abs.(A .- i_val))[2]
    end
    return index
end

function prepare_state_transition(m)
    if m.state_space == "KAΩ"
        Is_Π = repeat(collect(1:m.x_size),inner=m.K_size*m.Ω_size)
        Π_base = repeat(ω_trans[1,:],outer=m.K_size)'
        for i = 2:m.Ω_size
            Π_base = [Π_base; repeat(ω_trans[i,:],outer=m.K_size)']
        end
        Π_base = reshape(Π_base',(m.Ω_size*m.Ω_size*m.K_size))
        Xs_Π = repeat(Π_base,outer=m.K_size*m.A_size)

        Js_Π = zeros(Int, m.x_size*m.Ω_size*m.K_size)
        for l = 1:m.x_size
            A_ind = unique(A[Is[1+(l-1)*m.Ω_size:l*m.Ω_size]])[1]
            if A_ind == m.A_size 
                A_ind -= 1 # A_size is absorbing state
            end
            for k = 1:m.K_size
                Js_Π[1+(l-1)*m.Ω_size*m.K_size+(k-1)*m.Ω_size:(l-1)*m.Ω_size*m.K_size+(k)*m.Ω_size] = 
                    collect(1:m.Ω_size)' .+ (A_ind .+ 1 .- 1)'*m.Ω_size .+ m.A_size*m.Ω_size*(k - 1)
            end
        end
    elseif m.state_space == "KΩ"
        Is_Π = repeat(collect(1:m.x_size),inner=m.K_size*m.Ω_size)
        Π_base = repeat(ω_trans[1,:],outer=m.K_size)'
        for i = 2:m.Ω_size
            Π_base = [Π_base; repeat(ω_trans[i,:],outer=m.K_size)']
        end
        Π_base = reshape(Π_base',(m.Ω_size*m.Ω_size*m.K_size))
        Xs_Π = repeat(Π_base,outer=m.K_size)

        Js_Π = zeros(Int, m.x_size*m.Ω_size*m.K_size)
        for l = 1:m.x_size
            for k = 1:m.K_size
                Js_Π[1+(l-1)*m.Ω_size*m.K_size+(k-1)*m.Ω_size:(l-1)*m.Ω_size*m.K_size+(k)*m.Ω_size] = 
                    collect(1:m.Ω_size)' .+ m.Ω_size*(k - 1)
            end
        end
    end
    return Is_Π, Js_Π, Xs_Π
end

function resample_data(df, d)
    df_resample = copy(df)
    df_resample = df_resample[d.game_start .≤ df_resample.tid .≤ d.game_end, :]
    mat= convert(Matrix, df_resample[[:tid,:K,:A,:Ω,:kprime,:exit]])
    ind = collect(1:size(df_resample,1))
    ind_boot = sample(ind, size(df_resample,1)) # Indexes of bootstrap resamples  
    mat_boot = mat[ind_boot,:] 
    nms_boot = ["tid", "K", "A", "Ω", "kprime", "exit"]
    data_boot = DataFrame(mat_boot, nms_boot)
    return data_boot
end

function build_data(df, d, m)
    df_tmp = copy(df)
    rename!(df_tmp, :omega=>:ω)
    dropmissing!(df_tmp, :ω)
    df_tmp.time = df_tmp.tid .- d.t0 .+ 1 
    df_tmp = df_tmp[d.game_start .≤ df_tmp.tid .≤ d.game_end, :]

    state_capacity = Int.(findnearest(K_df.K, df_tmp.K))
    state_age = Int.(findnearest(A_df.A, df_tmp.A))
    state_productivity = Int.(findnearest(Ω_df.Ω, exp.(df_tmp.ω)))
    if m.state_space == "KAΩ"
        state = [state_capacity state_age state_productivity]
        state = ((state[:,1] .- 1)*m.A_size*m.Ω_size .+ (state[:,2].-1)*m.Ω_size .+ state[:,3])
    elseif m.state_space == "KΩ"
        state = [state_capacity state_productivity]
        state = ((state[:,1] .- 1)*m.Ω_size .+ state[:,2])    
    end
    decision_capacity_prime = Int.(findnearest(K_df.K, df_tmp.Kprime))
    decision_investment = decision_capacity_prime
    
    decision_stay = 1 .- df_tmp.exit
    decision = decision_investment.*decision_stay

    timevar = df_tmp.time

    # Number of entry from data
    df_tmp = copy(df)
    rename!(df_tmp, :omega=>:ω)
    dropmissing!(df_tmp, :ω)
    df_tmp.time = df_tmp.tid .- d.t0 .+ 1 
    df_tmp_K = df_tmp[d.game_start+1 .≤ df_tmp.tid .≤ d.game_end+1 .&& df_tmp.A .== 1, [:time, :K]]
    sort!(df_tmp_K, [:time, :K])
    df_tmp_K.time .-= 1
    df_tmp_K.K_index = Int.(findnearest(K_df.K, df_tmp_K.K))

    Nᵖᵉ = max.(m.N̄ .- sum(sMat_data[d.game_start-d.t0+1:d.game_end-d.t0+1,:],dims=2),0)

    timevar_ent = df_tmp_K.time 
    decision_ent = df_tmp_K.K_index 

    n_entrants = combine(groupby(df_tmp_K,:time), nrow => :K)
    rename!(n_entrants, :K => :entrants)
    decision_pe_quit = Nᵖᵉ .- n_entrants.entrants

    return timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit
end

std_norm_cdf(x::T) where {T <: Real} = 0.5 * erfc(-x/sqrt(2))
std_norm_cdf(x::Array{T}) where {T <: Real} = 0.5 .* erfc(-x./sqrt(2))

function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::T3=3) where {T1 <: Real, T2 <: Real, T3 <: Real}

    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    #y = range(ω_min, stop=ω_max, length=N)
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # fill in the middle columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    # NOTE: I need to shift this vector after finding probabilities
    #       because when finding the probabilities I use a function
    #       std_norm_cdf that assumes its input argument is distributed
    #       N(0, 1). After adding the mean E[y] is no longer 0, so
    #       I would be passing elements with the wrong distribution.
    #
    #       It is ok to do after the fact because adding this constant to each
    #       term effectively shifts the entire distribution. Because the
    #       normal distribution is symmetric and we just care about relative
    #       distances between points, the probabilities will be the same.
    #
    #       I could have shifted it before, but then I would need to evaluate
    #       the cdf with a function that allows the distribution of input
    #       arguments to be [μ/(1 - ρ), 1] instead of [0, 1]

    yy = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable

    # renormalize. In some test cases the rows sum to something that is 2e-15
    # away from 1.0, which caused problems in the MarkovChain constructor
    Π = Π./sum(Π, dims = 2)

    n, m = size(Π)

    return yy, Π
    #MarkovChain(Π, yy)
end

function tauchen(μ,σ²,ρ,znum,zdev)
    
    # Input:  zₜ₊₁ = ρzₜ + ωₜ₊₁, where ωₜ ~ N(μ,σ²) : note that stationary distribution of z is N( μ/(1-ρ), σ²/sqrt(1-ρ²) )
    #         Create transition matrix for the discrete markov process approximation of AR(1) process, by Tauchens method
    #         zdev:  the upper and lower points of the AR (1) variable are calculed as μ +- zdev*σ (3 by Tauchen)
    #         znum:  number of states in discretization of z (must be an odd number)

    # Output: AR (1) process z
    #         transition matrix: tranmat Πz,z' (given z, transtion to z')
    
    sigma  = sqrt(σ²);                 # stddev of ω
    zmean  = μ/(1-ρ);                  # expected value of z (mean)
    zsigma = sigma/sqrt(1-ρ^2);        # stddev of z
      
    z = zmean*ones(znum,1) .+ LinRange(-zdev*zsigma, zdev*zsigma, znum);   
    d = z[2]-z[1];                     # Note: all the points are equidistant by construction

    zi = z*ones(1,znum);
    zj = ones(znum,1)*transpose(z);
    
    # conditional on zᵢ, probability it falls in zj,  Note: mⱼ = (zⱼ₊₁ + zⱼ)/2
    P_part_up   = cdf.(Normal(),( ( zj .+ d/2 .- ( ρ.*zi .+ μ )) ./ sigma ));        
    P_part_down = cdf.(Normal(),( ( zj .- d/2 .- ( ρ.*zi .+ μ )) ./ sigma ));
  
    # change end points 
    P         = P_part_up - P_part_down;
    P[:,1]    = P_part_up[:,1];
    P[:,znum] = 1 .-P_part_down[:,znum];
    
    return z, P
end 
