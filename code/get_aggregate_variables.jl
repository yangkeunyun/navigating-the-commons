println("Load aggregate variables")

agg_output = CSV.read(path*raw"/Data/Intermediate/agg_output.csv", DataFrame) 
agg_capacity = CSV.read(path*raw"/Data/Intermediate/agg_capacity.csv", DataFrame) 
agg_incumbents = CSV.read(path*raw"/Data/Intermediate/agg_incumbents.csv", DataFrame)
pop_df = CSV.read(path*raw"/Data/Intermediate/whalepop_estimated.csv", DataFrame) 
usgdp_df = CSV.read(path*raw"/Data/Intermediate/aggregate_gdp_pop_index.csv", DataFrame) 
whale_prices = CSV.read(path*raw"/Data/Intermediate/demand_data.csv", DataFrame)[:,[:year,:price_whale_catchval_cpi]] 
dropmissing!(whale_prices)
rename!(whale_prices, :price_whale_catchval_cpi => :p)
petroleum = CSV.read(path*raw"/Data/Intermediate/petroleum_price.csv", DataFrame)[:,[:year,:petroleum_price_ma]] 
rename!(petroleum, :petroleum_price_ma => :petroleum_price)
petroleum = [DataFrame(year=collect(d.t0:d.tP-1), petroleum_price=ones(d.tP-d.t0)); petroleum]

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

# Find whale population model parameters
WMin = 100; r_max0 = 0.0135; z0 = 1.4; 
function find_whale_pop_params(w_param)
    WVec = zeros(eltype(m.κ),Tbar+1)
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




#=
PData[94:Tbar] = exp.((1/m.α[1]).*log.(QData[94:Tbar]) 
     - (m.α[2]/m.α[1])*log.(Pop[94:Tbar]) 
     - (m.α[3]/m.α[1]*log.(GDP[94:Tbar])) 
     - (m.α[4]/m.α[1])*collect(1909-1800:d.t1-1800) 
     - (m.α[5]/m.α[1])*Post1859[94:Tbar]
     - (m.α[6]/m.α[1])*Post1859[94:Tbar].*Pet[94:Tbar]
    .- (m.α[7]/m.α[1])) =#
#plot(whale_prices[d.t0 .≤ whale_prices.year .≤ 1905, :p],PData[1:90],seriestype=:scatter)
#plot(collect(d.t0:d.t1),PData)
#plot(whale_prices.year,whale_prices.p)

#plot(agg_incumbents.year,agg_incumbents.agg_incumbents)


