#= Production parameters obtained from static production function estimation =#
println("Production parameters")
PF_estimates = CSV.read(path*raw"/data/PF_CD.csv", DataFrame)

if m.state_space == "KAΩ"                
    βₖ = PF_estimates[PF_estimates.var.=="bk",:][!,2][1]
    βₐ = PF_estimates[PF_estimates.var.=="ba",:][!,2][1]
    βᵏ = PF_estimates[PF_estimates.var.=="bK",:][!,2][1]
    βʷ = PF_estimates[PF_estimates.var.=="bw",:][!,2][1]
    βᵗ = PF_estimates[PF_estimates.var.=="bt",:][!,2][1]
    β₀ = PF_estimates[PF_estimates.var.=="c",:][!,2][1]
    λ = PF_estimates[PF_estimates.var.=="rho",:][!,2][1]
else
    βₖ = PF_estimates[PF_estimates.var.=="bk",:][!,2][1]
    βᵏ = PF_estimates[PF_estimates.var.=="bK",:][!,2][1]
    βʷ = PF_estimates[PF_estimates.var.=="bw",:][!,2][1]
    βᵗ = PF_estimates[PF_estimates.var.=="bt",:][!,2][1]
    β₀ = PF_estimates[PF_estimates.var.=="c",:][!,2][1]
    λ = PF_estimates[PF_estimates.var.=="rho",:][!,2][1]
end