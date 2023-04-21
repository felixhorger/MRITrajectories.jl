
"""
"""
# TODO: this belongs into MRITrajectories bc not directly related to reconstruction, but trajectory design
function swap_duplicates_dynamic!(split_indices::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}, num_dynamic::Integer; maxiter::Integer=0) where N
	if maxiter == 0
		maxiter = sum(length.(split_indices))
	end
	for i = 1:num_dynamic
		indices_dynamic = split_indices[i]
		k = zero(CartesianIndex{N})
		j = 1
		while j ≤ length(indices_dynamic)
			m = indices_dynamic[j]
			if k == m
				found = false
				swap_dynamics = randperm!(collect(1:num_dynamic))
				for i_swap in swap_dynamics
					i_swap == i && continue
					indices_dynamic_swap = split_indices[i_swap]
					swap_indices = randperm!(collect(1:length(indices_dynamic_swap)))
					for j_swap in swap_indices
						k_swap = indices_dynamic_swap[j_swap]
						l = searchsortedfirst(indices_dynamic, k_swap)
						l_swap = searchsortedfirst(indices_dynamic_swap, k)
						l ≤ length(indices_dynamic) && indices_dynamic[l] == k_swap && continue # Is k_swap already in indices_dynamic?
						l_swap ≤ length(indices_dynamic_swap) && indices_dynamic_swap[l_swap] == k && continue # Is k already in indices_dynamic_swap?
						# Swap
						insert!(indices_dynamic, l, k_swap)
						insert!(indices_dynamic_swap, l_swap, k)
						deleteat!(indices_dynamic, (l > j) ? j : (j + 1))
						deleteat!(indices_dynamic_swap, (l_swap > j_swap) ? j_swap : (j_swap + 1))
						found = true
						break
					end
					found && break
				end
				if !found
					error("Not implemented, maybe do +- 1 search or reset to a default value, or do a hardcore search")
				end
			else
				j += 1
			end
			k = m
		end
	end
	return split_indices
end



"""
	floors the floats and adds one
"""
function cartesianise(points::AbstractVector{NTuple{N, Float64}}) where N
	[CartesianIndex(floor.(Int, p) .+ 1) for p in points]
end


function sort_constantL2distance!(split_sampling::AbstractVector{<: AbstractVector{<: CartesianIndex{N}}}, Δ::Real) where N
	@assert length(split_sampling) > 2
	@assert Δ ≥ 0
	num_dynamic = length(split_sampling)
	l = length.(split_sampling)
	for i = 1:sum(l)
		t = mod1(i, num_dynamic)
		j = (i-1) ÷ num_dynamic + 1
		k = Tuple(split_sampling[t][j])
		# Iterate over indices in next dynamic
		r = mod1(t+1, num_dynamic)
		j > l[r] && break
		δ = Inf
		q = 0
		for p = j:l[r]
			δp = abs(Δ - sqrt(sum(abs2, k .- Tuple(split_sampling[r][p]))))
			δp ≥ δ && continue
			δ = δp
			q = p
		end
		# Swap
		split_sampling[r][q], split_sampling[r][j] = split_sampling[r][j], split_sampling[r][q]
	end
	return split_sampling
end

