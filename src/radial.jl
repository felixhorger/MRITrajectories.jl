
function required_num_spokes(num_lines::Integer)
	num_spokes = floor(Int, 0.5π * num_lines)
	if iseven(num_spokes)
		num_spokes += 1
	end
	return num_spokes
end

"""
	Assumes the samples are symmetrically arranged around k = 0,
	i.e. if num_r is even, no sample hits k = 0, if num_r is odd, the centre one does.

"""
macro calculate_spokes(num_φ, φ)
	return esc(quote
		r = reshape(range(-π, π; length=num_r), 1, num_r)
		# Allocate memory
		k = Array{Float64, 3}(undef, 2, num_r, $num_φ)
		e_r = Vector{Float64}(undef, 2) # Unit vector pointing in radial direction
		# Iterate φ
		for spoke = 1:$num_φ
			sine, cosine = sincos($φ)
			e_r[1] = cosine
			e_r[2] = sine
			@. k[:, :, spoke] = r * e_r
		end
		return k
	end)
end
radial_spokes(φ::AbstractVector{<: Real}, num_r::Integer) = @calculate_spokes(length(φ), φ[spoke])
radial_spokes(num_φ::Integer, num_r::Integer) = @calculate_spokes(num_φ, 2π / num_φ * (spoke - 1))


"""
	Sort angles for dynamic imaging to minimise the difference between consecutive angles
	φ[spokes per dynamic, dynamic], indices same

"""
function sort_angles!(
	φ::AbstractVector{<: Real},
	indices::AbstractVector{<: Integer},
	spokes_per_dynamic::Integer,
	num_dynamic::Integer
)
	(φ, indices) = reshape.((φ, indices), spokes_per_dynamic, num_dynamic)
	@views for p = 1:num_dynamic
		perm = sortperm(φ[:, p])
		permute!(φ[:, p], perm)
		permute!(indices[:, p], perm)
	end
	return φ, indices
end

chronological_order(indices::AbstractMatrix) = vec(transpose(indices))

"""
	Stack of stars, fully sampled along partition direction.
	Spoke indices are increased until end is reached, then decreased for the next partition,
	continuing this alternating scheme for all partitions.
	This minimises changes in gradients and thus minimises eddy currents.

"""
function stack_of_stars(spoke_indices::AbstractMatrix{<: Integer}, num_partitions::Integer)
	spokes_per_pulse, num_dynamic = size(spoke_indices)
	incr = 1
	spoke = 1
	k = 1
	sampling = Vector{CartesianIndex{2}}(undef, spokes_per_pulse * num_dynamic * num_partitions)
	for partition = 1:num_partitions
		while (incr > 0 && spoke ≤ spokes_per_pulse) || (incr < 0 && spoke ≥ 1)
			for t = 1:num_dynamic
				sampling[k] = CartesianIndex(spoke_indices[spoke, t], partition)
				k += 1
			end
			spoke += incr
		end
		incr *= -1
		spoke += incr
	end
	return sampling
end

