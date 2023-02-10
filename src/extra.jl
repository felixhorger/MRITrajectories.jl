
"""
	Method of Dwork2020

	radius is minimum radius, if radius kwarg not given it is assumed constant

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
	floors the floats and adds one
"""
function cartesianise(points::AbstractVector{NTuple{N, Float64}}) where N
	[CartesianIndex(floor.(Int, p) .+ 1) for p in points]
end

