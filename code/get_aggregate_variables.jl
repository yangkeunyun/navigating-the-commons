struct AggregateVariables
    QData::Vector{Float64}
    KData::Vector{Float64}
    NData::Vector{Int64}
    WData::Vector{Float64}
    W₀::Vector{Float64}
    W_post::Vector{Float64}
    Pop::Vector{Int64}
    GDP::Vector{Int64}
    Pet::Vector{Float64}
    Post1859::Vector{Float64}
    PData::Vector{Float64}
    r_max::Float64
    z::Float64
end

function get_aggregate_variables(
    agg_output::DataFrame,
    agg_capacity::DataFrame,
    agg_incumbents::DataFrame,
    pop_df::DataFrame,
    usgdp_df::DataFrame,
    whale_prices::DataFrame,
    petroleum::DataFrame,
    d::DG_data,
    Tbar::Int64,
    WMin = 100,
    r_max0 = 0.0135,
    z0 = 1.4
    )::AggregateVariables

    println("Load aggregate variables")

    #= Get aggregate variables =#
    QData = agg_output[d.t0 .≤ agg_output.year .≤ d.t1+1,:whale_catch] 
    KData = agg_capacity[d.t0 .≤ agg_capacity.year .≤ d.t1+1,2] 
    NData = agg_incumbents[d.t0 .≤ agg_incumbents.year .≤ d.t1+1,2] 
    WData = pop_df.pop_baleen[d.t0 .≤ pop_df.year .≤ d.t1+1] + pop_df.pop_sperm[d.t0 .≤ pop_df.year .≤ d.t1+1]
    W₀ = pop_df.pop_baleen[pop_df.year .== d.t0] + pop_df.pop_sperm[pop_df.year .== d.t0]
    W_post = pop_df.pop_baleen[pop_df.year .== d.tP] + pop_df.pop_sperm[pop_df.year .== d.tP]
    Pop = usgdp_df.resident_population[d.t0 .≤ usgdp_df.year .≤ d.t1+1]
    GDP = usgdp_df.real_gdp_per_capita[d.t0 .≤ usgdp_df.year .≤ d.t1+1]
    Pet = petroleum.petroleum_price[d.t0 .≤ petroleum.year .≤ d.t1+1]

    Post1859 = zeros(Tbar)
    Post1859[d.tP-d.t0+1:end] .= 1

    #plot(whale_prices.year, whale_prices.p)
    PData = whale_prices[d.t0 .≤ whale_prices.year .<= d.t1, :p]

    function find_whale_pop_params(w_param)
        WVec = zeros(Tbar+1)
        WVec[[1]] = W₀
        for t = 1:(Tbar)
            WVec[[t+1]] .= WVec[t] .+ w_param[1]*WVec[t]*(1 .- (WVec[t]./W₀).^w_param[2]) .- QData[t]*1.71
            if WVec[t+1] < 0
                QData[[t]] .= max((WVec[t] + w_param[1]*WVec[t]*(1 - (WVec[t]/W₀[1])^w_param[2]))/1.71 - WMin, .1)
                WVec[[t+1]] .= WVec[t]
            end   
        end
        gap = sum(abs.(WVec - WData))
        return gap
    end
    
    initial = [r_max0, z0]
    optimum = optimize(find_whale_pop_params, initial)
    w_param_hat = optimum.minimizer

    r_max = w_param_hat[1]; z = w_param_hat[2]

    println("typeof(QData) = ", typeof(QData))
    println("typeof(KData) = ", typeof(KData))
    println("typeof(NData) = ", typeof(NData))
    println("typeof(WData) = ", typeof(WData))
    println("typeof(W₀) = ", typeof(W₀))
    println("typeof(W_post) = ", typeof(W_post))
    println("typeof(Pop) = ", typeof(Pop))
    println("typeof(GDP) = ", typeof(GDP))
    println("typeof(Pet) = ", typeof(Pet))
    println("typeof(Post1859) = ", typeof(Post1859))
    println("typeof(PData) = ", typeof(PData))
    println("typeof(r_max) = ", typeof(r_max))
    println("typeof(z) = ", typeof(z))
    
    return AggregateVariables(QData, KData, NData, WData, W₀, W_post, Pop, GDP, Pet, Post1859, PData, r_max, z)

end