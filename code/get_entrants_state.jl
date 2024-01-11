Kᵉ = df_tmp.K[df_tmp.A .== 1]
Ωᵉ = df_tmp.Ω[df_tmp.A .== 1]
Kᵉ = Int(findnearest(K_df.K, median(Kᵉ))[1])
#if Kᵉ == 1
#    Kᵉ += 1
#end     
Ωᵉ = Int(findnearest(Ω_df.Ω, median(skipmissing(Ωᵉ)))[1])
if m.state_space == "KAΩ"
    entrants_index = DataFrame(K_index=Kᵉ,Ω_index=Ωᵉ,A_index=1)
    entrants_index = leftjoin(entrants_index,x[:,[:K_index,:Ω_index,:A_index,:state]], on=[:K_index,:Ω_index,:A_index])
elseif m.state_space == "KΩ"
    entrants_index = DataFrame(K_index=Kᵉ,Ω_index=Ωᵉ)
    entrants_index = leftjoin(entrants_index,x[:,[:K_index,:Ω_index,:state]], on=[:K_index,:Ω_index])
end 
xᵉ = entrants_index.state[1]

Ω_ent_index = sort!(Int.(findnearest(Ω_df.Ω,df_tmp.Ω[df_tmp.A .== 1])))
Π_ent = [count(==(i), Ω_ent_index) for i in unique(Ω_ent_index)]./length(Ω_ent_index)


Kᵉˣ = df_tmp.K[df_tmp.exit .== 1]
Ωᵉˣ = df_tmp.Ω[df_tmp.exit .== 1]
Kᵉˣ = Int(findnearest(K_df.K, median(Kᵉˣ))[1])
Ωᵉˣ = Int(findnearest(Ω_df.Ω, median(skipmissing(Ωᵉˣ)))[1])

