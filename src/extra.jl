
"""
	Method of Dwork2020

	r is minimum radius, if radius kwarg not given it is assumed constant

	kwarg radius(x::NTuple{N, Real})

	kwarg k determines how many new points are generated per active point

	outputs points in [0, shape[d]]

	Periodic boundaries

	Note: problem for constant radius, it is not guaranteed that every cell is occupied
"""
function poissondisk(shape::NTuple{N, Integer}, r::Real; radius::Function=(_ -> r), k::Integer=30) where N
	# Set up arrays
	rmin = r
	grid_unit = rmin
	grid_shape = ceil.(Int, shape ./ grid_unit)
	points = Array{NTuple{N, Float64}}(undef, grid_shape)
	nan_tuple = ntuple(_ -> NaN, N)
	points .= (nan_tuple,)
	active = let 
		first_point = (shape .- eps(Float64)) .* ntuple(_ -> rand(), N) # shape is excluded
		first_point_grid = floor.(Int, first_point ./ grid_unit) .+ 1
		points[first_point_grid...] = first_point
		Set{CartesianIndex{N}}((CartesianIndex(first_point_grid),))
	end
	num_points = 1
	offset = Float64.(shape .÷ 2)
	# Generate points
	while length(active) > 0
		# Pick random active point
		i = rand(active)
		x = points[i]
		# "Radius" for this point
		r = radius(x .- offset)
		@assert r ≥ rmin
		rsq = r^2
		# Create k new points in the [r, 2r] range and check them
		one_is_valid = false
		for κ = 1:k
			y = nan_tuple
			while !(0.25^2 ≤ sum(abs2, y) ≤ 0.5^2)
				y = ntuple(_ -> rand() - 0.5, N)
			end
			y = mod.(4r .* y .+ x, shape)
			# Note: shape used here as exclusive upper limit because shape - ε would be mapped to shape-1,
			# otherwise maximum would be shape-2
			y_grid = CartesianIndex(floor.(Int, y ./ grid_unit .+ 1))
			# Check
			any(!isnan, points[y_grid]) && continue
			m = ceil(Int, r / rmin)
			lower = m + 1
			ι = y_grid - CartesianIndex(ntuple(_ -> lower, N))
			neighbours = 2m + 1
			is_valid = true
			for j in CartesianIndices(ntuple(_ -> neighbours, N)) # Edges are checked anyways despite not necessary
				l = ι + j
				l = CartesianIndex(mod1.(Tuple(l), grid_shape))
				#For non-periodic boundaries this can be used: !checkbounds(Bool, points, l) && continue
				z = points[l]
				any(isnan, z) && continue
				Δ = periodic_distance.(abs.(y .- z), shape)
				if sum(abs2, Δ) < rsq
					is_valid = false
					break
				end
			end
			if is_valid
				points[y_grid] = y
				push!(active, y_grid)
				one_is_valid = true
				num_points += 1
			end
		end
		# Remove x if no valid points where found
		!one_is_valid && pop!(active, i)
	end
	points = vec(points)
	valid_points = Vector{NTuple{N, Float64}}(undef, num_points)
	j = 1
	for i in eachindex(points)
		point = points[i]
		any(isnan, point) && continue
		valid_points[j] = point
		j += 1
	end
	return valid_points
end
@inline function periodic_distance(Δ::Real, boundary::Real)
	@assert Δ ≥ 0
	if Δ > boundary ÷ 2
		return boundary - Δ
	else
		return Δ
	end
end


function distribute_dynamically(indices::Vector{CartesianIndex{N}}, num_dynamic::Integer) where N
	num_samples = length(indices)
	num_samples_per_dynamic = num_samples ÷ num_dynamic
	upper_extra = mod(length(indices), num_dynamic) + 1
	split_indices = [Vector{CartesianIndex{N}}(undef, 0) for t = 1:num_dynamic]
	lengths = num_samples_per_dynamic .+ (0, 1)
	foreach(t -> sizehint!.(split_indices, lengths[1 + (t < upper_extra)]), 1:num_dynamic)
	sorted_indices = sort(indices)
	k = 1
	min_length = 0
	occupied = zeros(Bool, num_dynamic)
	while k ≤ length(sorted_indices)
		K = sorted_indices[k]
		l = k + 1
		while l ≤ length(sorted_indices)
			sorted_indices[l] != K && break
			l += 1
		end
		n = l - k
		occupied .= false
		while n > 0
			for t in randperm!(collect(1:num_dynamic))
				n == 0 && break
				this_length = length(split_indices[t])
				(occupied[t] || this_length == lengths[1 + (t < upper_extra)] || this_length > min_length) && continue
				push!(split_indices[t], K)
				occupied[t] = true
				n -= 1
			end
			min_length = minimum(length.(split_indices))
		end
		k = l
	end
	return split_indices
end

#function distribute_dynamically(indices::Vector{CartesianIndex{N}}, num_dynamic::Integer) where N
#	num_samples = length(indices)
#	num_samples_per_dynamic = num_samples ÷ num_dynamic
#	upper_extra = mod(length(indices), num_dynamic)
#	split_indices = [Vector{CartesianIndex{N}}(undef, num_samples_per_dynamic + 1 - (t > upper_extra)) for t = 1:num_dynamic]
#	assigned = zeros(Bool, length(indices))
# make assigned integer and remember where you put it, use that to swap
#	for t = 1:num_dynamic
#		inds = split_indices[t]
#		j = 0
#		num_required_samples = num_samples_per_dynamic + 1 - (t > upper_extra)
#		for (i, k) in enumerate(indices)
#			j == num_required_samples && break
#			(assigned[i] || k in inds) && continue
#			push!(inds, k)
#			assigned[i] = true
#			j += 1
#		end
#		j == num_required_samples && continue
#		# Need to swap
#		for n = j+1:num_required_samples
#			for t_swap = 1:t-1
#				inds_swap = indices[t_swap]
#				for (n_swap, k_swap) in inds_swap
#					k_swap in inds && continue
#					for (i, k) in enumerate(indices)
#						(assigned[i] || k in inds_swap) && continue
#						push!(inds, k)
#						assigned[i] = true
#						break
#					end
#				end
#			end
#		end
#	end
#	return split_indices
#end



"""
	sample_density() = (x,y,z,...) where x, y, z, ... ∈ [0, 1]
	samples outside the specified region will be rejected without error
	returns a Vector{NTuple{N, Float64}} 
"""
function rand(Int, sample_density::Function, num::Integer, shape::NTuple{N, Integer}; maxiter=100num) where N
	mask = zeros(Bool, shape)
	valid_points = Vector{CartesianIndex{N}}(undef, num)
	p = 0
	while p < num
		points = rand(sample_density, num - p, shape; maxiter)
		for point in points
			index = CartesianIndex(floor.(Int, point .+ 1))
			if mask[index] == 0
				valid_points[p+1] = index
				mask[index] = 1
				p += 1
			end
		end
	end
	return valid_points
end

# TODO: min max instead of shape
function rand(sample_density::Function, num::Integer, shape::NTuple{N, Integer}; maxiter=100num) where N
	points = Vector{NTuple{N, Float64}}(undef, num)
	i = 1
	j = 0
	upper = CartesianIndex(shape)
	while i ≤ num
		# Check maximum iterations
		j += 1
		if j == maxiter
			error("Maximum number of iterations exceeded")
			return indices
		end
		# Sample
		sample = sample_density()
		!all(0 .≤ sample .< 1) && continue
		points[i] = shape .* sample
		i += 1
	end
	return points
end


"""
	Points sorted along time
	For dynamic sampling:
	if maxi > num_dynamic, the sampling pattern won't be unique (same k-t sampled multiple times).
	You need to apply some randomisation to make it useful for dynamic sampling.
"""
function rand(Int, density::Function, sample_density::Function, num::Integer, shape::NTuple{N, Integer}, mini::Integer, maxi::Integer) where N
	# Determine normalisation sum(density)
	density_sum = 0.0
	for i in CartesianIndices(shape)
		d = density(Tuple(i) ./ shape)
		density_sum += d
	end
	normalisation = num / density_sum
	# Fill in first round of samples
	counts = zeros(Int, shape)
	for i in CartesianIndices(shape)
		counts[i] = min(max(floor.(Int, normalisation * density(Tuple(i) ./ shape)), mini), maxi)
	end
	# Randomly put residual samples
	s = sum(counts)
	while s < num
		for i in rand(Int, sample_density, num - s, shape)
			counts[i] == maxi && continue
			counts[i] += 1
		end
		s = sum(counts)
	end
	# Make indices
	points = Vector{CartesianIndex{N}}(undef, num)
	j = 1
	for i in CartesianIndices(shape)
		n = counts[i]
		while n > 0
			points[j] = i
			j += 1
			n -= 1
		end
	end
	return points, counts
end

function rand(Int, sample_density::Function, num::Integer, shape::NTuple{N, Integer}, mini::Integer, maxi::Integer; maxiter=1000_000) where N
	mask = fill(mini, shape)
	valid_points = Vector{CartesianIndex{N}}(undef, num)
	num_spatial = prod(shape)
	@assert mini * num_spatial ≤ num
	if mini > 0
		valid_points[1:num_spatial] .= shuffle(vec(CartesianIndices(shape)))
	end
	if mini > 1
		@views for i = 2:mini
			j = (i-1) * num_spatial
			valid_points[j+1 : j+num_spatial] = valid_points[1:num_spatial]
			shuffle!(valid_points[1:num_spatial])
		end
	end
	p = num_spatial * mini
	i = 1
	while p < num && i < maxiter
		points = rand(sample_density, num - p, shape; maxiter)
		for point in points
			index = CartesianIndex(floor.(Int, point .+ 1))
			if mask[index] < maxi
				valid_points[p+1] = index
				mask[index] += 1
				p += 1
			end
		end
		i += 1
	end
	return valid_points
end

"""
regular_dynamic(density::Function, shape::NTuple{N, Integer}, num_dynamic::Integer, num_samples::Integer, mini::Integer, maxi::Integer) where N
density must give the pixelwise probability of finding a sample there
"""
#function regular_dynamic(density::Function, shape::NTuple{N, Integer}, num_dynamic::Integer, num_samples::Integer, mini::Integer, maxi::Integer) where N
# TODO: regular dynamic for non-uniform sampling
#	rings of the below
#	sampling = regular_dynamic(shape, num_dynamic, maxi)
#
#	split_sampling = MRIRecon.split_sampling(sampling, num_dynamic)
#	mask = dropdims(MRIRecon.sampling_mask(sampling, num_dynamic); dims=1)
#	for t = 1:num_dynamic
#		for i = 1:length(split_sampling[t])
#			I = split_sampling[t][i]
#			p = density(Tuple(I) ./ shape)
#			if sum(mask[I, :]) / maxi > p
#				# too many samples in I, need to remove some
#			end
#			split_sampling
#		end
#	end
#	return sampling
#end

function disk(shape::NTuple{N, Integer}) where N
	centre = shape .÷ 2 .+ 0.5
	indices = Vector{CartesianIndex{N}}(undef, 0)
	sizehint!(indices, ceil(Int, 0.25π * prod(shape))) 
	for I in CartesianIndices(shape)
		sum(abs2, (Tuple(I) .- centre) ./ shape) > 0.25 && continue
		push!(indices, I)
	end
	return indices
end


