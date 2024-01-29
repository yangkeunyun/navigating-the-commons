
function calc_transition_productivity(capacity_prime, m)
    #capacity_prime = capacity_prime_cond
    #= Compute transition probabilities =#
    Js = zeros(Int, m.x_size*m.Ω_size)
    if m.state_space == "KAΩ"
        Threads.@threads for l = 1:m.x_size
            A_ind = unique(A[Is[1+(l-1)*m.Ω_size:l*m.Ω_size]])[1]
            if A_ind == m.A_size 
                A_ind -= 1 # A_size is absorbing state
            end
            Js[1+(l-1)*m.Ω_size:l*m.Ω_size] = collect(1:m.Ω_size)' .+ (A_ind .+ 1 .- 1)'*m.Ω_size .+ m.A_size*m.Ω_size*(capacity_prime[l] - 1)
        end
    elseif m.state_space == "KΩ"
        Threads.@threads for l = 1:m.x_size
            Js[1+(l-1)*m.Ω_size:l*m.Ω_size] = collect(1:m.Ω_size)' .+ m.Ω_size*(capacity_prime[l] - 1)
        end    
    end
    Π_Ω = sparse(Is, Js, Xs, m.x_size,m.x_size) # the model is solved using policy iteration with a sparse transition matrix since ∃ (x_size) states
    return Π_Ω
end


function calc_transition_matrix!(Π, ιMat1, m, t)
    Xs_Π_update = zeros(eltype(m.Γ),size(Xs_Π,1))
    if m.state_space == "KAΩ"
        Threads.@threads for l = 1:m.x_size
            A_ind = unique(A[Is[1+(l-1)*m.Ω_size:l*m.Ω_size]])[1]
            if A_ind == m.A_size 
                A_ind -= 1 # A_size is absorbing state
            end
            for k = 1:m.K_size
                Xs_Π_update[1+(l-1)*m.Ω_size*m.K_size+(k-1)*m.Ω_size:(l-1)*m.Ω_size*m.K_size+(k)*m.Ω_size] = 
                    Xs_Π[1+(l-1)*m.Ω_size*m.K_size+(k-1)*m.Ω_size:(l-1)*m.Ω_size*m.K_size+(k)*m.Ω_size]*ιMat1[t,l,k]
            end
        end
    elseif m.state_space == "KΩ"
        Threads.@threads for l = 1:m.x_size
            for k = 1:m.K_size
                Xs_Π_update[1+(l-1)*m.Ω_size*m.K_size+(k-1)*m.Ω_size:(l-1)*m.Ω_size*m.K_size+(k)*m.Ω_size] = 
                    Xs_Π[1+(l-1)*m.Ω_size*m.K_size+(k-1)*m.Ω_size:(l-1)*m.Ω_size*m.K_size+(k)*m.Ω_size]*ιMat1[t,l,k]
            end
        end
    end 
    Π .= sparse(Is_Π, Js_Π, Xs_Π_update, m.x_size,m.x_size) # the model is solved using policy iteration with a sparse transition matrix since ∃ (x_size) states
    return nothing
end

# 0.000152 seconds (27 allocations: 11.250 KiB)
function calc_production!(
    Q::Vector{Float64},
    KVec::Vector{Float64},
    WVec::Vector{Float64},
    d::DG_data,
    t::Int64,
    state_space::String,
    production_parameters::ProductionParameters,
    x::DataFrame
)::Nothing
    @unpack βₖ, βₐ, βᵏ, βʷ, βᵗ, β₀, λ = production_parameters
    K = x.K
    Ω = x.Ω
    if "A" ∈ names(x)
        A = x.A
    end


    if state_space == "KAΩ"
        Q .= exp.(βₖ.*log.(K) + βₐ.*log.(A) + log.(Ω) .+ βᵏ.*log.(KVec[t]) .+ βʷ.*log.(WVec[t]) .+ βᵗ.*(d.t0 - 1800 - 1 + t) .+ β₀)
    elseif state_space == "KΩ"
        Q .= exp.(βₖ.*log.(K) + log.(Ω) .+ βᵏ.*log.(KVec[t]) .+ βʷ.*log.(WVec[t]) .+ βᵗ.*(d.t0 - 1800 - 1 + t) .+ β₀)
    end
    return nothing
end    

function compute_aggregate_state!(
    NVec::Vector{Float64},
    KVec::Vector{Float64},
    QVec::Vector{Float64},
    QᶠVec::Vector{Float64},
    WVec::Vector{Float64},
    PVec::Vector{Float64},
    sMat::Matrix{Float64},
    Q::Vector{Float64},
    d::DG_data,
    m::DG_model,
    production_parameters::ProductionParameters,
    agg_variables::AggregateVariables,
    demand_parameters_pre::DemandParameters,
    demand_parameters_post::DemandParameters,
    x::DataFrame,
    df::DataFrame
)

    @unpack QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z = agg_variables
    
    @unpack αᵖ, αᵖᵒᵖ, αᵍᵈᵖ, αᵗ, α₀ = demand_parameters_pre
    αᵖ_pre, αᵖᵒᵖ_pre, αᵍᵈᵖ_pre, αᵗ_pre, α₀_pre = αᵖ, αᵖᵒᵖ, αᵍᵈᵖ, αᵗ, α₀

    @unpack αᵖ, αᵖᵒᵖ, αᵍᵈᵖ, αᵗ, αᵖᵉᵗ = demand_parameters_post
    αᵖ_post, αᵖᵒᵖ_post, αᵍᵈᵖ_post, αᵗ_post, αᵖᵉᵗ_post = αᵖ, αᵖᵒᵖ, αᵍᵈᵖ, αᵗ, αᵖᵉᵗ

    K = x.K
    Ω = x.Ω
    if "A" ∈ names(x)
        A = x.A
    end

    PVec_base = log(QData[1859 - d.t0 + 1]) .- αᵖ_post*log(PData[1859 - d.t0 + 1]) - αᵖᵒᵖ_post*log(Pop[1859 - d.t0 + 1]) - αᵍᵈᵖ_post*log.(GDP[1859 - d.t0 + 1]) - αᵖᵉᵗ_post*log(Pet[1859 - d.t0 + 1]) - αᵗ_post*(d.t0 - 1800 - 1 + 1859 - d.t0 + 1) 
    for t = 1:(d.game_end - d.t0 + 6) #d.game_start -d.t0 
        NVec[[t]] .= sum(sMat[t,:])
        KVec[[t]] .= sum(K .* sMat[t,:])
        K_missing = sum(df[df.tid .== t + d.t0 - 1 .&& ismissing.(df.fid), [:K]].K)
        KVec[[t]] .+= K_missing
        if KVec[t] < 562 # minimum aggregate fleet size in Tower (1907)
            KVec[t] = 562
        end

        calc_production!(Q, KVec, WVec, d, t, m.state_space, production_parameters, x)
        QVec[[t]] .= sum(Q .* sMat[t,:])
        #QVec[[t]] .= sum(Q .* sMat_data[t,:])
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
        
        if t ≤ (d.game_end - d.t0 + 1)
            if m.post == 0
                PVec[[t]] .= exp.((1/αᵖ_pre).*log.(QVec[t]) - (αᵖᵒᵖ_pre/αᵖ_pre)*log.(Pop[t]) - (αᵍᵈᵖ_pre/αᵖ_pre)*log.(GDP[t]) - (αᵗ_pre/αᵖ_pre)*(d.t0 - 1800 - 1 + t) .- (α₀_pre/αᵖ_pre))
            else 
                PVec[[t]] .= exp.((1/αᵖ_post).*log.(QVec[t]) - (αᵖᵒᵖ_post/αᵖ_post)*log.(Pop[t]) - (αᵍᵈᵖ_post/αᵖ_post)*log.(GDP[t]) - (αᵗ_post/αᵖ_post)*(d.t0 - 1800 - 1 + t) - (αᵖᵉᵗ_post/αᵖ_post)*Post1859[t].*log.(Pet[t]) .- (PVec_base/αᵖ_post))
            end            
        else
            PVec[[t]] .= (1+(((PVec[t-1] / PVec[t-2] - 1) + (PVec[t-2] / PVec[t-3] - 1) + (PVec[t-3] / PVec[t-4] - 1) + (PVec[t-4] / PVec[t-5] - 1) + (PVec[t-5] / PVec[t-6] - 1))/5))*PVec[t-1]
        end
    end
    
    return nothing
end


#0.000127 seconds (43 allocations: 53.250 KiB)
function calc_cost(ι, m)
    #ι = invest_cond
    cost = zeros(eltype(m.Γ),m.x_size)
    if m.cost_struct == 1
        cost[vec(ι) .> 0] = (m.Γ[1]*ι[ι .> 0] + m.Γ[2].*(ι[ι .> 0]).^2 )
        cost[vec(ι) .< 0] = (m.Γ[3]*ι[ι .< 0] + m.Γ[4].*(ι[ι .< 0]).^2 )
    elseif m.cost_struct == 2
        cost[vec(ι) .> 0] = (m.Γ[1]*ι[ι .> 0] + m.Γ[2].*(ι[ι .> 0]).^2 )
        cost[vec(ι) .< 0] = (m.Γ[3]*ι[ι .< 0] + m.Γ[4].*(ι[ι .< 0]).^2 )
    elseif m.cost_struct == 3
        cost[vec(ι) .≥ 0] = (m.Γ[1]*ι[ι .≥ 0])
        cost[vec(ι) .< 0] = (m.Γ[2]*ι[ι .< 0])
    elseif m.cost_struct == 4
        cost[vec(ι) .> 0] = (m.Γ[1]*(ι[ι .> 0]).^2)
        cost[vec(ι) .< 0] = (m.Γ[2]*(ι[ι .< 0]).^2)
    elseif m.cost_struct == 5
        cost[vec(ι) .> 0] = (m.Γ[1]*(ι[ι .> 0]).^2) + m.Γ[2]*(ι[ι .> 0])
        cost[vec(ι) .< 0] = (m.Γ[1]*(ι[ι .< 0]).^2) + m.Γ[3]*(ι[ι .< 0])
    elseif m.cost_struct == 6
        cost[vec(ι) .≥ 0] = m.Γ[1] .+ m.Γ[2]*ι[ι .≥ 0]./K[ι .≥ 0]
        cost[vec(ι) .< 0] = m.Γ[3] .+ m.Γ[4]*ι[ι .< 0]./K[ι .< 0]
    elseif m.cost_struct == 7
        cost[vec(ι) .≥ 0] = m.Γ[1] .+ m.Γ[2]*ι[ι .≥ 0]./K[ι .≥ 0] + m.Γ[3]*(ι[ι .≥ 0]./K[ι .≥ 0]).^2 
        cost[vec(ι) .< 0] = m.Γ[4] .+ m.Γ[5]*ι[ι .< 0]./K[ι .< 0] - m.Γ[6]*(ι[ι .< 0]./K[ι .< 0]).^2 
    elseif m.cost_struct == 8
        cost[vec(ι) .≥ 0] = m.Γ[1]*K[ι .≥ 0].*(ι[ι .≥ 0]./K[ι .≥ 0]).^2 
        cost[vec(ι) .< 0] = m.Γ[2]*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2 
    elseif m.cost_struct == 9
        cost[vec(ι) .> 0] = m.Γ[1]*m.scale_P*m.scale_Q + m.Γ[2]*K[ι .> 0].*(ι[ι .> 0]./K[ι .> 0]).^2 
        cost[vec(ι) .< 0] = m.Γ[3]*m.scale_P*m.scale_Q + m.Γ[4]*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2 
    elseif m.cost_struct == 10
        cost[vec(ι) .> 0] = m.Γ[1] .+ ι[ι .> 0] + m.Γ[2]*(ι[ι .> 0]).^2 
        cost[vec(ι) .< 0] = m.Γ[3] .+ ι[ι .< 0] + m.Γ[4]*(ι[ι .< 0]).^2     
    elseif m.cost_struct == 11
        cost[vec(ι) .> 0] = m.Γ[1] .+ m.Γ[2]*ι[ι .> 0] + m.Γ[3]*(ι[ι .> 0]).^2
        cost[vec(ι) .< 0] = m.Γ[4] .+ m.Γ[5]*ι[ι .< 0] + m.Γ[6]*(ι[ι .< 0]).^2
    elseif m.cost_struct == 12
        cost[vec(ι) .> 0] = m.Γ[1] .+ (m.Γ[2]*K[ι .> 0].*(ι[ι .> 0]./K[ι .> 0]).^2)*m.scale_Q
        cost[vec(ι) .< 0] = m.Γ[3] .+ (m.Γ[4]*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2)*m.scale_Q
    elseif m.cost_struct == 13
        cost[vec(ι) .> 0] = (m.Γ[1]/2)*K[ι .> 0].*(ι[ι .> 0]./K[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0]
        cost[vec(ι) .< 0] = (m.Γ[3]/2)*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2 + m.Γ[4]*ι[ι .< 0]
    elseif m.cost_struct == 14
        cost[vec(ι) .> 0] = m.Γ[1]*(ι[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0]
        cost[vec(ι) .< 0] = m.Γ[1]*(ι[ι .< 0]).^2 + m.Γ[3]*ι[ι .< 0]
    elseif m.cost_struct == 15
        cost[vec(ι) .> 0] = m.Γ[1]*ι[ι .> 0]
        cost[vec(ι) .< 0] = m.Γ[2]*K[ι .< 0] + m.Γ[3]*(ι[ι .< 0]).^2
    elseif m.cost_struct == 16
        cost[vec(ι) .> 0] = m.Γ[1]*(ι[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0] .+ m.Γ[4]
        cost[vec(ι) .< 0] = m.Γ[1]*(ι[ι .< 0]).^2 + m.Γ[3]*ι[ι .< 0] .+ m.Γ[5]
    elseif m.cost_struct == 17
        cost[vec(ι) .> 0] = m.Γ[1]*(ι[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0] + m.Γ[4]./K[ι .> 0]
        cost[vec(ι) .< 0] = m.Γ[1]*(ι[ι .< 0]).^2 + m.Γ[3]*ι[ι .< 0] + m.Γ[4]./K[ι .< 0]
    end
    return cost
end

# 0.000112 seconds (18 allocations: 32.797 KiB)
function calc_cost_ext(ι, m)
    #ι = invest_cond
    if m.cost_struct == 1
        cost = (m.Γ[3]*ι[ι .< 0] + m.Γ[4].*(ι[ι .< 0]).^2 )
    elseif m.cost_struct == 2
        cost = (m.Γ[3]*ι[ι .< 0] + m.Γ[4].*(ι[ι .< 0]).^2 )
    elseif m.cost_struct == 3
        cost = (m.Γ[2]*ι[ι .< 0])
    elseif m.cost_struct == 4
        cost = (m.Γ[2]*(ι[ι .< 0]).^2)
    elseif m.cost_struct == 5
        cost = (m.Γ[1]*(ι[ι .< 0]).^2) + m.Γ[3]*(ι[ι .< 0])
    elseif m.cost_struct == 6
        cost = m.Γ[3] .+ m.Γ[4]*ι[ι .< 0]./K[ι .< 0]
    elseif m.cost_struct == 7
        cost = m.Γ[4] .+ m.Γ[5]*ι[ι .< 0]./K[ι .< 0] - m.Γ[6]*(ι[ι .< 0]./K[ι .< 0]).^2 
    elseif m.cost_struct == 8
        cost = m.Γ[2]*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2 
    elseif m.cost_struct == 9
        cost = m.Γ[3]*m.scale_P*m.scale_Q + m.Γ[4]*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2 
    elseif m.cost_struct == 10
        cost = m.Γ[3] .+ ι[ι .< 0] + m.Γ[4]*(ι[ι .< 0]).^2     
    elseif m.cost_struct == 11
        cost = m.Γ[4] .+ m.Γ[5]*ι[ι .< 0] + m.Γ[6]*(ι[ι .< 0]).^2
    elseif m.cost_struct == 12
        cost = m.Γ[3] .+ (m.Γ[4]*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2)*m.scale_Q
    elseif m.cost_struct == 13
        cost = (m.Γ[3]/2)*K[ι .< 0].*(ι[ι .< 0]./K[ι .< 0]).^2 + m.Γ[4]*ι[ι .< 0]
    elseif m.cost_struct == 14
        cost = m.Γ[1]*(ι[ι .< 0]).^2 + m.Γ[3]*ι[ι .< 0]
    elseif m.cost_struct == 15
        cost = m.Γ[2]*K[ι .< 0] + m.Γ[3]*(ι[ι .< 0]).^2
    elseif m.cost_struct == 16
        cost = m.Γ[1]*(ι[ι .< 0]).^2 + m.Γ[3]*ι[ι .< 0] # remove adjustment fixed cost to identify exit fixed cost
    elseif m.cost_struct == 17
        cost = m.Γ[1]*(ι[ι .< 0]).^2 + m.Γ[3]*ι[ι .< 0] + m.Γ[4]./K[ι .< 0]
    end
    return cost
end


function calc_cost_ent(ι, m)
    #ι = invest_cond
    if m.cost_struct == 1
        cost = (m.Γ[1]*ι[ι .> 0] + m.Γ[2].*(ι[ι .> 0]).^2 )
    elseif m.cost_struct == 2
        cost = (m.Γ[1]*ι[ι .> 0] + m.Γ[2].*(ι[ι .> 0]).^2 )
    elseif m.cost_struct == 3
        cost = (m.Γ[1]*ι[ι .≥ 0])
    elseif m.cost_struct == 4
        cost = (m.Γ[1]*(ι[ι .> 0]).^2)
    elseif m.cost_struct == 5
        cost = (m.Γ[1]*(ι[ι .> 0]).^2) + m.Γ[2]*(ι[ι .> 0])
    elseif m.cost_struct == 6
        cost = m.Γ[1] .+ m.Γ[2]*ι[ι .≥ 0]./K[ι .≥ 0]
    elseif m.cost_struct == 7
        cost = m.Γ[1] .+ m.Γ[2]*ι[ι .≥ 0]./K[ι .≥ 0] + m.Γ[3]*(ι[ι .≥ 0]./K[ι .≥ 0]).^2 
    elseif m.cost_struct == 8
        cost = m.Γ[1]*K[ι .≥ 0].*(ι[ι .≥ 0]./K[ι .≥ 0]).^2 
    elseif m.cost_struct == 9
        cost = m.Γ[1]*m.scale_P*m.scale_Q + m.Γ[2]*K[ι .> 0].*(ι[ι .> 0]./K[ι .> 0]).^2 
    elseif m.cost_struct == 10
        cost = m.Γ[1] .+ ι[ι .> 0] + m.Γ[2]*(ι[ι .> 0]).^2 
    elseif m.cost_struct == 11
        cost = m.Γ[2]*ι[ι .> 0] + m.Γ[3]*(ι[ι .> 0]).^2 # remove investment fixed cost to identify entry fixed cost
    elseif m.cost_struct == 12
        cost = m.Γ[1] .+ (m.Γ[2]*K[ι .> 0].*(ι[ι .> 0]./K[ι .> 0]).^2)*m.scale_Q
    elseif m.cost_struct == 13
        cost = (m.Γ[1]/2)*K[ι .> 0].*(ι[ι .> 0]./K[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0]
    elseif m.cost_struct == 14
        cost = m.Γ[1]*(ι[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0]
    elseif m.cost_struct == 15
        cost = m.Γ[1]*ι[ι .> 0]
    elseif m.cost_struct == 16
        cost = m.Γ[1]*(ι[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0] # remove adjustment fixed cost to identify entry fixed cost
    elseif m.cost_struct == 17
        cost = m.Γ[1]*(ι[ι .> 0]).^2 + m.Γ[2]*ι[ι .> 0] + m.Γ[4]./K[ι .> 0]
    end
    return cost
end

        
function compute_distribution!(sMat, ιMat1, λVec1, d, m)
    Π = spzeros(eltype(m.Γ),m.x_size, m.x_size)
    for t = (d.game_start - d.t0 + 1):(d.game_end - d.t0 + 6)
        calc_transition_matrix!(Π, ιMat1,m,t)
        sMat[t+1,:] .= Π'sMat[t,:]
        Nᵖᵉ = max.(m.N̄ - sum(sMat[t,:]),0)
        tmp_ent = repeat(λVec1[t,:], inner=m.Ω_size)
        tmp_pro = repeat(Π_ent, outer=m.K_size)
        sMat[t+1, :] += Nᵖᵉ*tmp_ent.*tmp_pro # new entrants
    end  
    return nothing
end


# 0.002739 seconds (26.27 k allocations: 14.199 MiB)
function get_EV_CCP_inc(period_profit, V, m, K_space, K)
    #V = VMat[t+1,:] # for test
    #period_profit = πMat[t,:]
    soical_surplus_stay = zeros(eltype(m.Γ),m.x_size,m.K_size)
    pr_capacity_prime = zeros(m.x_size,m.K_size)
    EVᵐᵃˣ = zeros(eltype(m.Γ),m.x_size)
    Threads.@threads for j = 1:m.K_size 
        capacity_prime_cond = j*ones(Int,m.x_size)
        invest_cond = K_space[capacity_prime_cond] - (1-m.δ)*K
        cost = calc_cost(invest_cond, m)
        Π_cond = calc_transition_productivity(capacity_prime_cond,m) 
        soical_surplus_stay[:,j] = -(cost)/(m.scale_P*m.scale_Q) + m.ρ*Π_cond*V 
    end

    EVᵐᵃˣ = maximum((soical_surplus_stay)/m.ψ, dims=2)
    Υ = exp.(((m.ϕᶠ .- calc_cost_ext(- (1-m.δ)*K, m))/(m.scale_P*m.scale_Q))/m.ψ .- EVᵐᵃˣ) .+ sum(exp.(((soical_surplus_stay))/m.ψ .- EVᵐᵃˣ), dims=2)
    social_surplus = period_profit + m.ψ*(0.5772156649 .+ EVᵐᵃˣ + log.(Υ)) # Euler's constant ≃ 0.5772156649
    pr_capacity_prime = exp.(((soical_surplus_stay))/m.ψ .- EVᵐᵃˣ)./Υ 
    pr_exit = exp.(((m.ϕᶠ .- calc_cost_ext(- (1-m.δ)*K, m))/(m.scale_P*m.scale_Q))/m.ψ .- EVᵐᵃˣ)./Υ 

    return social_surplus, pr_capacity_prime, pr_exit
end

# 0.000508 seconds (666 allocations: 26.969 KiB)
function get_CCP_ent(V, m)    
    #V = VMat[t+1,:] # for test
    soical_surplus_ent = zeros(eltype(m.Γ),1,m.K_size)
    pr_capacity_prime = zeros(1,m.K_size)
    EVᵐᵃˣ = zeros(eltype(m.Γ),m.x_size)
    Threads.@threads for j = 1:m.K_size 
        invest_cond = K_space[j]
        cost = calc_cost_ent(invest_cond, m)
        soical_surplus_ent[:,j] .= -(cost)/(m.scale_P*m.scale_Q) + m.ρ*Π_ent'*V[1+(j-1)*m.Ω_size:j*m.Ω_size] - (m.κ/(m.scale_P*m.scale_Q))
    end

    EVᵐᵃˣ = maximum((soical_surplus_ent)/m.ψ, dims=2)
    Υ = exp.((0.0/(m.scale_P*m.scale_Q))/m.ψ .- EVᵐᵃˣ) .+ sum(exp.(((soical_surplus_ent))/m.ψ  .- EVᵐᵃˣ), dims=2)
    pr_entry = exp.(((soical_surplus_ent))/m.ψ .- EVᵐᵃˣ)./Υ 
    
    return pr_entry
end


function compute_backward!(
    ιMat1::Array{Float64, 3},
    χMat1::Matrix{Float64},
    λVec1::Matrix{Float64},
    πMat::Matrix{Float64},
    VMat::Matrix{Float64},
    Q::Vector{Float64},
    mc::Vector{Float64},
    KVec::Vector{Float64},
    WVec::Vector{Float64},
    PVec::Vector{Float64},
    d::DG_data,
    m::DG_model,
    production_parameters::ProductionParameters,
    agg_variables::AggregateVariables,
    x::DataFrame,
    K_space::Vector{Float64},
    K::Vector{Float64}
)

    for t = (d.game_end - d.t0 + 6):-1:(d.game_start - d.t0 + 1) #Tbar:-1:(d.game_start - d.t0 + 1)
        # Calculate the number of whales harvested and profit
        calc_production!(Q, KVec, WVec, d, t, m.state_space, production_parameters, x)
        #calc_mc!(mc, KVec, WVec, d, t)
        πMat[t,:] = ((PVec[t]/m.scale_P)*(Q./m.scale_Q))*.3 #- (PVec[t]./m.scale_P)*mc

        # Calculate stationary value functions
        if t == (d.game_end - d.t0 + 6) #Tbar
        #if t == Tbar
            Ṽ = ((1)/(1-m.ρ))*πMat[t,:] #(m.ρ^(d.game_end - d.t0 + 1 + 1)/(1-m.ρ))*π̃(1-m.ρ^10)
            VMat[t+1,:] = Ṽ
        end

        # Incumbents' problem
        VMat[t,:], ιMat1[t,:,:], χMat1[t,:] = get_EV_CCP_inc(πMat[t,:], VMat[t+1,:], m, K_space, K)

        # Potential entrants' problem
        λVec1[t,:] = get_CCP_ent(VMat[t+1,:], m)    
    end
    return nothing
end

function solve_NOE(
    d::DG_data,
    m::DG_model,
    Tbar::Int64,
    agg_variables::AggregateVariables,
    sMat_data::Matrix{Float64},
    K::Vector{Float64},
    production_parameters::ProductionParameters,
    x::DataFrame,
    demand_parameters_pre::DemandParameters,
    demand_parameters_post::DemandParameters,
    df::DataFrame,
    n_max::Int64,
    K_space::Vector{Float64}
)
    @unpack QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z = agg_variables
    if m.cost_struct == 1 || m.cost_struct == 10
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ⁺₁ = $(m.Γ[1]), γ⁺₂ = $(m.Γ[2]), γ⁻₁ = $(m.Γ[3]), γ⁻₂ = $(m.Γ[4])")
    elseif m.cost_struct == 2 || m.cost_struct == 9 || m.cost_struct == 12
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ⁺ = $(m.Γ[1]), γ⁺₂ = $(m.Γ[2]), γ⁻ = $(m.Γ[3]), γ⁻₂ = $(m.Γ[4]), γᶠ = $(m.Γ[5])")
    elseif m.cost_struct == 3
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ⁺ = $(m.Γ[1]), γ⁻ = $(m.Γ[2])")
    elseif m.cost_struct == 4 || m.cost_struct == 8
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ⁺ = $(m.Γ[1]), γ⁻ = $(m.Γ[2])")#, γᶠ = $(m.Γ[3])
    elseif m.cost_struct == 5
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ² = $(m.Γ[1]), γ⁺ = $(m.Γ[2]), γ⁻ = $(m.Γ[3])")   
    elseif m.cost_struct == 11
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ⁺ = $(m.Γ[1]), γ⁺₁ = $(m.Γ[2]), γ⁺₂ = $(m.Γ[3]), γ⁻ = $(m.Γ[4]), γ⁻₁ = $(m.Γ[5]), γ⁻₂ = $(m.Γ[6])")
    elseif m.cost_struct == 13
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ²₊ = $(m.Γ[1]), γ¹₊ = $(m.Γ[2]), γ²₋ = $(m.Γ[3]), γ¹₋ = $(m.Γ[4])")  
    elseif m.cost_struct == 14
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ² = $(m.Γ[1]), γ⁺ = $(m.Γ[2]), γ⁻ = $(m.Γ[3])")   
    elseif m.cost_struct == 15
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ¹₊ = $(m.Γ[1]), γ¹₋ = $(m.Γ[2]), γ²₋ = $(m.Γ[3])")  
    elseif m.cost_struct == 16
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ² = $(m.Γ[1]), γ⁺ = $(m.Γ[2]), γ⁻ = $(m.Γ[3]), γᶠ₊ = $(m.Γ[4]), γᶠ₋ = $(m.Γ[5])")    
    elseif m.cost_struct == 17
        println("\nParameters: ψ = $(m.ψ), κ = $(m.κ), ϕ = $(m.ϕ), ϕᶠ = $(m.ϕᶠ), γ² = $(m.Γ[1]), γ⁺ = $(m.Γ[2]), γ⁻ = $(m.Γ[3]), γᶠ = $(m.Γ[4])")    
    end
  
    λVec0 = ones(eltype(m.κ),Tbar + 1,m.K_size); λVec0[end] = 0            
    #d.ιMat0 = ones(Tbar+1,m.x_size,5)
    ιMat0 = ones(eltype(m.κ),Tbar+1,m.x_size,m.K_size)
    χMat0 = ones(eltype(m.κ),Tbar+1,m.x_size)
    λVec1 = ones(eltype(m.κ),Tbar+1,m.K_size); λVec1[end] = 0
    #d.ιMat1 = ones(Tbar+1,m.x_size,5)
    ιMat1 = ones(eltype(m.κ),Tbar+1,m.x_size,m.K_size)
    χMat1 = ones(eltype(m.κ),Tbar+1,m.x_size)

    # Initial matrices (period-by-state) and vectors (by period)
    sMat = zeros(eltype(m.κ),Tbar+1,m.x_size)   # industry state, i.e. the number of firms in each possible state
    PVec = zeros(eltype(m.κ),Tbar+1)            # price vector
    WVec = zeros(eltype(m.κ),Tbar+1)            # whale population in each period
    KVec = zeros(eltype(m.κ),Tbar+1)            # aggregate capacity in each period
    QVec = zeros(eltype(m.κ),Tbar+1)            # aggregate whale catch in each period
    QᶠVec = zeros(eltype(m.κ),Tbar+1)           # aggregate whale catch outside of US whaling
    NVec = zeros(eltype(m.κ),Tbar+1)            # aggregate number of incumbents in each period
    VMat = zeros(eltype(m.κ),Tbar+1,m.x_size)
    πMat = zeros(eltype(m.κ),Tbar+1,m.x_size)
    Q = zeros(eltype(m.κ),m.x_size)
    mc = zeros(eltype(m.κ),m.x_size)
    
    WVec[[1]] = W₀
    sMat[1,:] = sMat_data[1,:]
    NVec[[1]] .= max(sum(sMat[1,:]),1)
    KVec[[1]] .= sum(K .* sMat[1,:])
    calc_production!(Q, KVec, WVec, d, 1, m.state_space, production_parameters, x)
    QVec[[1]] .= sum(Q.*sMat[1,:])
    sMat[1:end-1,:] = sMat_data
    compute_aggregate_state!(NVec, KVec, QVec, QᶠVec, WVec, PVec, sMat, Q, d, m, production_parameters, agg_variables, demand_parameters_pre, demand_parameters_post, x, df)

    n = 0
    Δ = 10; Δ₁ = 10; Δ₂ = 10; Δ₃ = 10; 
    
    while ((Δ₁ > m.options.ϵ₁ || Δ₂ > m.options.ϵ₁ || Δ₃ > m.options.ϵ₁) && n ≤ n_max)
        n += 1
        if n == 1
            sMat[1:end-1,:] = sMat_data
            compute_aggregate_state!(NVec, KVec, QVec, QᶠVec, WVec, PVec, sMat, Q, d, m, production_parameters, agg_variables, demand_parameters_pre, demand_parameters_post, x, df)
        else 
            compute_distribution!(sMat, ιMat1, λVec1, d, m)
            compute_aggregate_state!(NVec, KVec, QVec, QᶠVec, WVec, PVec, sMat, Q, d, m, production_parameters, agg_variables, demand_parameters_pre, demand_parameters_post, x, df)
        end

        compute_backward!(ιMat1, χMat1, λVec1, πMat, VMat, Q, mc, KVec, WVec, PVec, d, m, production_parameters, agg_variables, x, K_space, K)

        ιMat1 = (ιMat0 + ιMat1)/2
        χMat1 = (χMat0 + χMat1)/2
        λVec1 = (λVec0 + λVec1)/2
    
        Δ₁ = maximum(abs.(ιMat0 - ιMat1))
        Δ₂ = maximum(abs.(χMat0 - χMat1))
        Δ₃ = maximum(abs.(λVec0 - λVec1))

        Δ = sum(abs.(χMat0 - χMat1)) + sum(abs.(ιMat0 - ιMat1)) + sum(abs.(λVec0 - λVec1))

        # Update policies
        ιMat0 .= ιMat1
        χMat0 .= χMat1
        λVec0 .= λVec1

        #println("n = $n; Δ: $Δ, Δ₁: $(Δ₁), Δ₂: $Δ₂, Δ₃: $Δ₃")
    end

    if n > n_max
        println("Stragtegies do not converge; Δ: $Δ, Δ₁: $Δ₁, Δ₂: $Δ₂, Δ₃: $Δ₃")
    else
        println("Stragtegies converge at n = $n; Δ: $Δ, Δ₁: $Δ₁, Δ₂: $Δ₂, Δ₃: $Δ₃")
    end

    return χMat1, ιMat1, λVec1, NVec, KVec, QVec, WVec, PVec
end

function calc_log_likelihood(
    timevar::Vector{Int64},
    state::Vector{Int64},
    decision::Vector{Int64},
    timevar_ent::Vector{Int64},
    decision_ent::Vector{Int64},
    decision_pe_quit::Matrix{Float64},
    d::DG_data,
    m::DG_model,
    Tbar::Int64,
    agg_variables::AggregateVariables,
    sMat_data::Matrix{Float64},
    K::Vector{Float64},
    production_parameters::ProductionParameters,
    x::DataFrame,
    demand_parameters_pre::DemandParameters,
    demand_parameters_post::DemandParameters,
    df::DataFrame,
    n_max::Int64,
    K_space::Vector{Float64}
)::Float64

    @unpack QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z = agg_variables

    exit_prob, invest_prob, entry_prob = solve_NOE(d, m, Tbar, agg_variables, sMat_data, K,production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space)[1:3]
    
    exit_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:][exit_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:] .== 0.0] = exit_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:][exit_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:] .== 0.0] .+ 1e-25
    invest_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:,:][invest_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:,:] .== 0.0] = invest_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:,:][invest_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:,:] .== 0.0] .+ 1e-25
    pll_inc = [if σ == 0 log(exit_prob[convert(Int,t),convert(Int,i)]) 
                else log(invest_prob[convert(Int,t),convert(Int,i),convert(Int,σ)]) end 
                for (t,i,σ) in zip(timevar, state, decision)]

    entry_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:][entry_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:] .== 0.0] = entry_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:][entry_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:] .== 0.0] .+ 1e-25
    pll_pe_ent = [log(entry_prob[convert(Int,t),convert(Int,σ)]) for (t,σ) in zip(timevar_ent, decision_ent)]
    
    quit_prob = 1 .- sum(entry_prob[d.game_start-d.t0+1:d.game_end - d.t0+1,:],dims=2)
    quit_prob[quit_prob .< 0.0] .= 1e-25; quit_prob[quit_prob .> 1.0] .= 1

    pll_pe_quit = decision_pe_quit[1:end].*log.(quit_prob)

    ll = sum(pll_inc) + sum(pll_pe_ent) + sum(pll_pe_quit)

    return ll
end    
#sortslices(tmp, dims = 1)
#tmp[tmp[:,1].==-Inf,:]

## Outer loop: search for parameters θ
function optimize_log_likelihood(
    timevar::Vector{Int64},
    state::Vector{Int64},
    decision::Vector{Int64},
    timevar_ent::Vector{Int64},
    decision_ent::Vector{Int64},
    decision_pe_quit::Matrix{Float64},
    d::DG_data,
    m::DG_model,
    Tbar::Int64,
    agg_variables::AggregateVariables,
    sMat_data::Matrix{Float64},
    x::DataFrame,
    production_parameters::ProductionParameters,
    demand_parameters_pre::DemandParameters,
    demand_parameters_post::DemandParameters,
    df::DataFrame,
    Nᵖᵉ::Int64,
    n_max::Int64,
    K_space::Vector{Float64}
)

    @unpack QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z = agg_variables
    K = x.K
    Ω = x.Ω
    if "A" ∈ names(x)
        A = x.A
    end
    
    function obj_func(θ)
        if m.scale_est == "No"
            m.κ = θ[1]; m.ϕᶠ = θ[2]; m.Γ = θ[3:end]
        else
            m.ψ = θ[1]; m.κ = θ[2]; m.ϕᶠ = θ[3]; m.Γ = θ[4:end]
        end

        pll = @time calc_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, K, production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space)
        if m.state_space == "KAΩ"                
            println("The dimension of state space is $(m.K_size*m.A_size*m.Ω_size)")
        elseif m.state_space == "KΩ"
            println("The dimension of state space is $(m.K_size*m.Ω_size) with $(m.K_size) grid for K and $(m.Ω_size) grid for Ω")
        end

        if d.game_start == 1804
            println("Cost structure: $(m.cost_struct); Scale: $(m.scale_Q); Ealry period partial log-likelihood: $pll")
        elseif 1804 < d.game_start < 1859
            println("Cost structure: $(m.cost_struct); Scale: $(m.scale_Q); Golden-age period partial log-likelihood: $pll")
        else
            println("Cost structure: $(m.cost_struct); Scale: $(m.scale_Q); Post period partial log-likelihood: $pll")
        end
        return -pll
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
        
        optimum = optimize(obj_func, 
                            lower,upper,
                            θ₀, 
                            Fminbox(NelderMead()),
                            Optim.Options(outer_iterations = 1500,iterations=10000)) #Fminbox(NelderMead()) NelderMead() SimulatedAnnealing()
        
    end
    
    return optimum
end


function calc_std_errors(θ, d, m)
    function log_likelihood(θ)
        if m.scale_est == "No"
            m.κ = θ[1]; m.ϕᶠ = θ[2]; m.Γ = θ[3:end]
        else
            m.ψ = θ[1]; m.κ = θ[2]; m.ϕᶠ = θ[3]; m.Γ = θ[4:end]
        end
        
        ll = calc_log_likelihood(timevar, state, decision, timevar_ent, decision_ent, decision_pe_quit, d, m, Tbar, agg_variables, sMat_data, K, production_parameters, x, demand_parameters_pre, demand_parameters_post, df, n_max, K_space)
        return ll
    end

    Hessian = ForwardDiff.hessian(log_likelihood, θ)
    var_cov = inv(-Hessian)
    if sum(diag(var_cov) .< 0) != 0
        println("There is negative diagonal element")
        println("The diagonal of variance-covariance is $(diag(var_cov))")
        SE2 = sqrt.(abs.(diag(var_cov)))
    else
        SE2 = sqrt.(diag(var_cov))
    end
    #return SE, SE2
    return SE2
end
