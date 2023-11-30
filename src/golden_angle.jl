
#= Convention on the "phis":
	ϕ is a phase
	φ is a polar angle (in 3D θ will be the polar angle and φ the azimuthal angle)
	Φ is the golden angle
=#

function golden_angle(n::Integer)
	@assert n > 0
	return π / (0.5 * (1 + √5) + n - 1)
end

function golden_angle_incremented(num::Integer, angles_per_2π::Integer; n_golden::Integer=1)
	Φ = golden_angle(n_golden)
	# Angle resolution
	resolution = 2π / angles_per_2π
	ξ = Φ / resolution # How many indices to jump if incremented by golden angle
	# Indices of angles in the range
	indices = @. mod(floor(Int, ξ * (0:num-1)), angles_per_2π)
	# Angles
	φ = @. indices * resolution
	# One based indexing
	indices .+= 1
	return φ, indices
end


function golden_means(dim::Integer)
	# Chan2009
	dim += 1
	Fib = zeros(dim, dim)
	Fib[dim+1:dim+1:end] .= 1
	Fib[dim, dim] = 1
	Fib[dim, 1] = 1
	@views begin
		F = eigen(Fib)
		@assert isapprox(imag(F.values[end]), 0.0)
		v = abs.(F.vectors[:, end])
		v = v[1:end-1] ./ v[end]
	end
	return v
end

