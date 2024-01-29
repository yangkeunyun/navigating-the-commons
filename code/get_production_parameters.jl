#= Production parameters obtained from static production function estimation =#

struct ProductionParameters
    βₖ::Float64
    βₐ::Float64
    βᵏ::Float64
    βʷ::Float64
    βᵗ::Float64
    β₀::Float64
    λ::Float64
end

function get_production_parameters(state_space::String, PF_estimates::DataFrame)
    println("Production parameters")
    if state_space == "KAΩ"                
        βₖ = PF_estimates[PF_estimates.var.=="bk",:][!,2][1]
        βₐ = PF_estimates[PF_estimates.var.=="ba",:][!,2][1]
        βᵏ = PF_estimates[PF_estimates.var.=="bK",:][!,2][1]
        βʷ = PF_estimates[PF_estimates.var.=="bw",:][!,2][1]
        βᵗ = PF_estimates[PF_estimates.var.=="bt",:][!,2][1]
        β₀ = PF_estimates[PF_estimates.var.=="c",:][!,2][1]
        λ = PF_estimates[PF_estimates.var.=="rho",:][!,2][1]
    else
        βₖ = PF_estimates[PF_estimates.var.=="bk",:][!,2][1]
        βₐ = 9e99
        βᵏ = PF_estimates[PF_estimates.var.=="bK",:][!,2][1]
        βʷ = PF_estimates[PF_estimates.var.=="bw",:][!,2][1]
        βᵗ = PF_estimates[PF_estimates.var.=="bt",:][!,2][1]
        β₀ = PF_estimates[PF_estimates.var.=="c",:][!,2][1]
        λ = PF_estimates[PF_estimates.var.=="rho",:][!,2][1]
    end
    return ProductionParameters(βₖ, βₐ, βᵏ, βʷ, βᵗ, β₀, λ)
end