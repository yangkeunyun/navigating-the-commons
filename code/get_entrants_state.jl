


function get_entrants_state(
    df::DataFrame,
    x::DataFrame,
    Ω_df::DataFrame,
    K_df::DataFrame,
    state_space::String
    )::Tuple{Int,Array{Float64,1},Int,Int}

    Kᵉ = df.K[df.A .== 1]
    Ωᵉ = df.Ω[df.A .== 1]
    Kᵉ = Int(findnearest(K_df.K, median(Kᵉ))[1])
        
    Ωᵉ = Int(findnearest(Ω_df.Ω, median(skipmissing(Ωᵉ)))[1])
    if state_space == "KAΩ"
        entrants_index = DataFrame(K_index=Kᵉ,Ω_index=Ωᵉ,A_index=1)
        entrants_index = leftjoin(entrants_index,x[:,[:K_index,:Ω_index,:A_index,:state]], on=[:K_index,:Ω_index,:A_index])
    elseif state_space == "KΩ"
        entrants_index = DataFrame(K_index=Kᵉ,Ω_index=Ωᵉ)
        entrants_index = leftjoin(entrants_index,x[:,[:K_index,:Ω_index,:state]], on=[:K_index,:Ω_index])
    end 
    xᵉ = entrants_index.state[1]
    
    Ω_ent_index = sort!(Int.(findnearest(Ω_df.Ω,df.Ω[df.A .== 1])))
    Π_ent = [count(==(i), Ω_ent_index) for i in unique(Ω_ent_index)]./length(Ω_ent_index)
    
    
    Kᵉˣ = df.K[df.exit .== 1]
    Ωᵉˣ = df.Ω[df.exit .== 1]
    Kᵉˣ = Int(findnearest(K_df.K, median(Kᵉˣ))[1])
    Ωᵉˣ = Int(findnearest(Ω_df.Ω, median(skipmissing(Ωᵉˣ)))[1])
    
    return xᵉ, Π_ent, Kᵉˣ, Ωᵉˣ
end