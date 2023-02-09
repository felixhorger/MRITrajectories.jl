
"""
	Method of Dwork2020

	radius is minimum radius, if radius kwarg not given it is assumed constant

	kwarg radius(x::NTuple{N, Real})

	kwarg k determines how many new points are generated per active point

	outputs points in [0, shape[d]]

	Note: problem for constant radius, it is not guaranteed that every cell is occupied
"""
function poisson(shape::NTuple{N, Integer}, r::Real; radius::Function=(_ -> r), k::Integer=30) where N
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
				#!checkbounds(Bool, points, l) && continue
				if sum(abs2, y .- points[l]) < rsq
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

"""
	floors the floats and adds one
"""
function points2indices(points::AbstractVector{NTuple{N, Float64}}) where N
	[CartesianIndex(floor.(Int, p) .+ 1) for p in points]
end

