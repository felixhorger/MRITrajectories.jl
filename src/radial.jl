
function required_num_spokes(num_lines::Integer)
	num_spokes = floor(Int, 0.5π * num_lines)
	if iseven(num_spokes)
		num_spokes -= 1
	end
	return num_spokes
end

oddify(n::Integer) = n + 1 - mod(n, 2)


"""
	Assumes the samples are symmetrically arranged around k = 0,
	i.e. if num_r is even, no sample hits k = 0, if num_r is odd, the centre one does.

	Be aware that if num_φ is even, then the same k are measured twice but in opposed directions

"""
macro calculate_spokes(num_φ, φ)
	return esc(quote
		r = reshape(range(-π + π/num_r, π - π/num_r; length=num_r), 1, num_r)
		#r = reshape(range(-π, π - 2π/num_r; length=num_r), 1, num_r)
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

macro calculate_spokes(num_samples, φ, θ)
	return esc(quote
		r = reshape(range(-π + π/num_r, π - π/num_r; length=num_r), 1, num_r)
		#r = reshape(range(-π, π - 2π/num_r; length=num_r), 1, num_r)
		# Allocate memory
		k = Array{Float64, 3}(undef, 3, num_r, $num_samples)
		e_r = Vector{Float64}(undef, 3) # Unit vector pointing in radial direction
		# Iterate φ
		for spoke = 1:$num_samples
			sineφ, cosineφ = sincos($φ)
			sineθ, cosineθ = sincos($θ)
			e_r[1] = cosineφ * sineθ
			e_r[2] = sineφ   * sineθ
			e_r[3] = cosineθ
			@. k[:, :, spoke] = r * e_r
		end
		return k
	end)
end
function radial_spokes(φ::AbstractVector{<: Real}, θ::AbstractVector{<: Real}, num_r::Integer)
	@assert length(φ) == length(θ)
	@calculate_spokes(length(φ), φ[spoke], θ[spoke])
end
function radial_spokes(num_φ::Integer, num_θ::Integer, num_r::Integer)
	num_samples = num_φ * num_θ
	@calculate_spokes(
		num_samples,
		2π / num_φ * (mod1(spoke, num_φ) - 1),
		π / (num_θ+1) * ((spoke-1) ÷ num_φ + 1)
	)
end



"""
	Sort angles for dynamic imaging to minimise the difference between consecutive angles
	φ[spokes per dynamic, dynamic], indices same

"""
function sort_angles!(
	φ::AbstractVector{<: Real}, # TODO clunky
	indices::AbstractVector{<: Integer},
	spokes_per_dynamic::Integer,
	num_dynamic::Integer
)
	@assert num_dynamic > 1
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
	spoke_indices[spokes_per_pulse, num_dynamic]
	Spoke indices are increased until end is reached, then decreased for the next partition,
	continuing this alternating scheme for all partitions.
	This minimises changes in gradients and thus minimises eddy currents.

	TODO: Would this also work well with spirals?
"""
function stack_of_stars(spoke_indices::AbstractMatrix{<: Integer}, partitions::AbstractVector{<: Integer}; version=:new)
	if version == :old
		@warn "Old stack of stars"
		spokes_per_pulse, num_dynamic = size(spoke_indices)
		incr = 1
		spoke = 1
		k = 1
		sampling = Vector{CartesianIndex{2}}(undef, spokes_per_pulse * num_dynamic * length(partitions))
		for partition in partitions
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
	elseif version == :new
		@warn "New stack of stars"
		spokes_per_pulse, num_dynamic = size(spoke_indices)
		k = 1
		sampling = Vector{CartesianIndex{2}}(undef, spokes_per_pulse * num_dynamic * length(partitions))
		for partition in partitions, spoke = 1:spokes_per_pulse, t = 1:num_dynamic
			sampling[k] = CartesianIndex(spoke_indices[spoke, t], partition)
			k += 1
		end
	else
		error("Unrecognised version")
	end
	return sampling
end


# num_θ referes to the number of angles in 2π, even though θ ∈ [0, π]
function golden_kooshball(num_samples::Integer, num_φ::Integer, num_θ::Integer)

	φ1, φ2 = golden_means(2)

	θ = Vector{Float64}(undef, num_samples)
	φ = similar(θ)
	idx = similar(θ, CartesianIndex{2})

	θ_resolution = 2π / num_θ
	φ_resolution = 2π / num_φ
	for i = 0:num_samples-1
		t = i + 1
		θ[t] = acos(mod(i * φ1, 1))
		φ[t] = 2π * mod(i * φ2, 1)
		@show n = floor(Int, θ[t] / θ_resolution)
		@show m = floor(Int, φ[t] / φ_resolution)
		θ[t] = n * θ_resolution
		φ[t] = m * φ_resolution
		idx[t] = CartesianIndex(n, m)
	end

	return φ, θ, idx
end

