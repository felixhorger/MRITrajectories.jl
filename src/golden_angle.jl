
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

function uniquify_indices!(idx::AbstractVector{<: Integer})
	used = zeros(Int, length(idx))

	for i in idx
		used[i] += 1
	end

	for (ii, i) in enumerate(idx)

		used[i] == 1 && continue
		local j_before = 0
		for outer j_before = i-1:-1:1
			used[j_before] != 0 && continue
			break
		end

		local j_after = length(idx) + 1
		for outer j_after = i+1:length(idx)
			used[j_after] != 0 && continue
			break
		end

		if i - j_before < j_after - i && j_before > 0 && used[j_before] == 0
			idx[ii] = j_before
			used[j_before] = 1
		else
			idx[ii] = j_after
			used[j_after] = 1
		end
		used[i] -= 1
	end
	return idx
end

