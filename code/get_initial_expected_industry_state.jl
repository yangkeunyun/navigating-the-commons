println("Get initial guess of expected industry state in each period directly from data")

df_tmp.Ω_index = findnearest(Ω_df.Ω, df_tmp.Ω) # describe(df.Ω_index)
df_tmp.K_index = findnearest(K_df.K, df_tmp.K) # describe(df.K_index)
df_tmp.A_index = findnearest(A_df.A, df_tmp.A) # describe(df.A_index)
if m.state_space == "KAΩ"
    s_df = df_tmp[:,[:K_index,:Ω_index,:A_index, :tid]]
    s_df = transform(groupby(s_df, [:K_index,:Ω_index,:A_index, :tid]), nrow => :s) # get the number of occurrences for each transition
    s_df = unique(s_df[:,[:tid, :K_index,:Ω_index,:A_index, :s]]) # remove unnecessary rows and columns  
    s_df = rightjoin(s_df,x[:,[:K_index,:Ω_index,:A_index,:state]], on=[:K_index,:Ω_index,:A_index])
elseif m.state_space == "KΩ"
    s_df = df_tmp[:,[:K_index,:Ω_index, :tid]]
    s_df = transform(groupby(s_df, [:K_index,:Ω_index, :tid]), nrow => :s) # get the number of occurrences for each transition
    s_df = unique(s_df[:,[:tid, :K_index,:Ω_index, :s]]) # remove unnecessary rows and columns  
    s_df = rightjoin(s_df,x[:,[:K_index,:Ω_index,:state]], on=[:K_index,:Ω_index])
end   
s_df = s_df[:, [:tid, :state, :s]]
sort!(s_df, [:state]) # Caution! It is needed for Linux
s_df = unstack(s_df, :state, :s)
sort!(s_df, :tid)
dropmissing!(s_df, :tid)
s_df = coalesce.(s_df, 0.0)
s_df = s_df[d.t0 .≤ s_df.tid .≤ d.t1,:]

# Find the missing years in 'tid'
missing_elements = setdiff(d.t0:maximum(s_df.tid), s_df.tid)

# Loop over state columns and fill missing values with zeros
if missing_elements != Int64[]
    # Create a new DataFrame with missing elements and zeros in 'y' columns
    result_df = DataFrame(tid = vcat(s_df.tid, missing_elements))
    for col in names(s_df, Not(:tid))
        missing_values = zeros(length(missing_elements))
        existing_values = s_df[findall(s_df.tid .== missing_elements[1]), col]
        if !isempty(existing_values)
            missing_values[1:length(existing_values)] .= existing_values
        end
        result_df[!, col] = vcat(s_df[!, col], missing_values)
    end
    # Sort the DataFrame based on 'x'
    result_df = sort!(result_df, :tid)
    # Finally, convert the dataframe to a matrix
    sMat_data = Tables.matrix(result_df[d.t0 .≤ result_df.tid .≤ d.t1,2:end])
else
    sMat_data = Tables.matrix(s_df[d.t0 .≤ s_df.tid .≤ d.t1,2:end])
end



