

function get_demand_parameters(dm_estimates::DataFrame)
    αᵖ = dm_estimates[dm_estimates.var.=="ln_price_whale_catchval_cpi",:][!,2][1]
    αᵖᵒᵖ = dm_estimates[dm_estimates.var.=="ln_resident_population",:][!,2][1]
    αᵍᵈᵖ = dm_estimates[dm_estimates.var.=="ln_real_gdp_per_capita",:][!,2][1]
    αᵗ = dm_estimates[dm_estimates.var.=="t",:][!,2][1]
    α₀ = dm_estimates[dm_estimates.var.=="_cons",:][!,2][1]
    return DemandParameters(αᵖ, αᵖᵒᵖ, αᵍᵈᵖ, αᵗ, 0.0, α₀)
end