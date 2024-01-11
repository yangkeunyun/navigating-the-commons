println("Demand parameters")

dm_estimates = CSV.read(path*raw"/Data/Intermediate/demand_curve.csv", DataFrame)
αᵖ = dm_estimates[dm_estimates.var.=="ln_price_whale_catchval_cpi",:][!,2][1]
αᵖᵒᵖ = dm_estimates[dm_estimates.var.=="ln_resident_population",:][!,2][1]
αᵍᵈᵖ = dm_estimates[dm_estimates.var.=="ln_real_gdp_per_capita",:][!,2][1]
αᵗ = dm_estimates[dm_estimates.var.=="t",:][!,2][1]
αᵖᵒˢᵗ = dm_estimates[dm_estimates.var.=="1.post_1859",:][!,2][1]
αᵖᵉᵗ = dm_estimates[dm_estimates.var.=="1.post_1859#c.ln_petroleum_price",:][!,2][1]
α₀ = dm_estimates[dm_estimates.var.=="_cons",:][!,2][1]

dm_estimates_pre = CSV.read(path*raw"/Data/Intermediate/demand_curve_pre.csv", DataFrame)
αᵖ_pre = dm_estimates_pre[dm_estimates_pre.var.=="ln_price_whale_catchval_cpi",:][!,2][1]
αᵖᵒᵖ_pre = dm_estimates_pre[dm_estimates_pre.var.=="ln_resident_population",:][!,2][1]
αᵍᵈᵖ_pre = dm_estimates_pre[dm_estimates_pre.var.=="ln_real_gdp_per_capita",:][!,2][1]
αᵗ_pre = dm_estimates_pre[dm_estimates_pre.var.=="t",:][!,2][1]
α₀_pre = dm_estimates_pre[dm_estimates_pre.var.=="_cons",:][!,2][1]

dm_estimates_post = CSV.read(path*raw"/Data/Intermediate/demand_curve_post.csv", DataFrame)
αᵖ_post = dm_estimates_post[dm_estimates_post.var.=="ln_price_whale_catchval_cpi",:][!,2][1]*2
αᵖ_post = αᵖ_pre
αᵖᵒᵖ_post = dm_estimates_post[dm_estimates_post.var.=="ln_resident_population",:][!,2][1]
αᵍᵈᵖ_post = dm_estimates_post[dm_estimates_post.var.=="ln_real_gdp_per_capita",:][!,2][1]
αᵗ_post = dm_estimates_post[dm_estimates_post.var.=="t",:][!,2][1]
αᵖᵉᵗ_post = dm_estimates_post[dm_estimates_post.var.=="1.post_1859#c.ln_petroleum_price",:][!,2][1]
α₀_post = dm_estimates_post[dm_estimates_post.var.=="_cons",:][!,2][1]

dm_estimates_pre_NEP = CSV.read(path*raw"/Data/Intermediate/demand_curve_pre_NEP.csv", DataFrame)
αᵖᵒᵖ_pre_NEP = dm_estimates_pre_NEP[dm_estimates_pre_NEP.var.=="ln_resident_population",:][!,2][1]
αᵍᵈᵖ_pre_NEP = dm_estimates_pre_NEP[dm_estimates_pre_NEP.var.=="ln_real_gdp_per_capita",:][!,2][1]
αᵗ_pre_NEP = dm_estimates_pre_NEP[dm_estimates_pre_NEP.var.=="t",:][!,2][1]
α₀_pre_NEP = dm_estimates_pre_NEP[dm_estimates_pre_NEP.var.=="_cons",:][!,2][1]

dm_estimates_post_NEP = CSV.read(path*raw"/Data/Intermediate/demand_curve_post_NEP.csv", DataFrame)
αᵖᵒᵖ_post_NEP = dm_estimates_post_NEP[dm_estimates_post_NEP.var.=="ln_resident_population",:][!,2][1]
αᵍᵈᵖ_post_NEP = dm_estimates_post_NEP[dm_estimates_post_NEP.var.=="ln_real_gdp_per_capita",:][!,2][1]
αᵗ_post_NEP = dm_estimates_post_NEP[dm_estimates_post_NEP.var.=="t",:][!,2][1]
αᵖᵉᵗ_post_NEP = dm_estimates_post_NEP[dm_estimates_post_NEP.var.=="ln_petroleum_price",:][!,2][1]
α₀_post_NEP = dm_estimates_post_NEP[dm_estimates_post_NEP.var.=="_cons",:][!,2][1]
