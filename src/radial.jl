
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
		# TODO is this correct?
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
		π / num_θ * ((spoke-1) ÷ num_φ)
	)
end



"""
	Sort angles for dynamic imaging to minimise the difference between consecutive angles
	φ[spokes per dynamic, dynamic], indices same
	works only if φ ∈ [0, 2π] TODO
	TODO: this doesn't work properly, need to manually find the best match, where the diff must be smaller than pi! (circle ...)

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
		perm = sortperm(φ[:, p]; by=β->(β > π ? β - 2π : β))
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


# TODO: convert from continuous angles to gridded ones:
#=
	φ_resolution = 2π / num_φ
	θ_resolution = π / num_θ
	n = mod1(floor(Int, φ[t] / φ_resolution) + 1, num_φ)
	m = mod1(floor(Int, θ[t] / θ_resolution) + 1, num_θ)
	φ[t] = n * φ_resolution
	θ[t] = m * θ_resolution
	idx[t] = CartesianIndex(n, m)
=#

# TODO: doesn't belong here? or 2D radial function in golden.jl has to moved here?
# num_θ referes to the number of angles in π
function golden_kooshball(num_samples::Integer)
	φ1, φ2 = golden_means(2)
	φ = Vector{Float64}(undef, num_samples)
	θ = similar(φ)
	for i = 0:num_samples-1
		t = i + 1
		φ[t] = 2π * mod(i * φ2, 1)
		θ[t] = acos(2 * mod(i * φ1, 1) - 1)
	end
	return φ, θ
end


# TODO some of this could go into an extra module? Since it's not only MRI?
# TODO: this excludes theta = 0 and pi, but if you round/ceil/floor you will end up measuring those angles anyways, how to proceed?
# Could just deal with that when generating the indices from the angles, that would be a better interface
# TODO: better name, how?
function uniform_kooshball(num_φ::Integer, num_θ::Integer)
	θ = Vector{Float64}(undef, 0)
	φ = Vector{Float64}(undef, 0)
	# Determine rough size:
	# Need: ∫_{-π/2}^{π/2} * cos θ dθ ≈ sum_{i=1}^{num_θ} cosθ_i Δθ / num_θ
	# Hence: sum_{i=1}^{num_θ} cosθ_i ≈ (∫_{-π/2}^{π/2} * cos θ dθ) * num_θ / Δθ
	# And: ∫_{-π/2}^{π/2} * cos θ dθ = 2
	# total num ≈ num_φ * sum_{i=1}^{num_θ} cosθ_i = num_φ * 2 * num_θ / π
	sizehint!.((φ, θ), ceil(Int, 2/π * num_φ * num_θ))
	v = Vector{Float64}(undef, num_φ)
	@views for ϑ in range(-π/2, π/2; length=num_θ+2)[2:end-1]
		n = round(Int, num_φ * cos(ϑ))
		n += iseven(n)
		v[1:n] .= ϑ + π/2
		append!(θ, v[1:n])
		append!(φ, range(0, 2π; length=n+1)[1:end-1])
	end
	return φ, θ
end

# See Saff1997
function generalised_spiral_kooshball(num::Integer)
	φ = Vector{Float64}(undef, num)
	θ = Vector{Float64}(undef, num)
	φ[1] = 0
	θ[1] = π
	φ[num] = 0
	θ[num] = 0
	for i = 2:num-1
		h = 2 * (i - 1) / (num - 1) - 1
		θ[i] = acos(h)
		φ[i] = mod2pi( φ[i-1] + 3.6 / √(num * (1 - h^2)) )
	end
	return φ, θ
end

function spherical2cartesian(
	φ::AbstractVector{<: Real},
	θ::AbstractVector{<: Real}
)
	@assert length(φ) == length(θ)
	v = Matrix{Float64}(undef, 3, length(φ))
	for i in eachindex(φ)
		sineφ, cosineφ = sincos(φ[i])
		sineθ, cosineθ = sincos(θ[i])
		v[1, i] = sineθ * cosineφ
		v[2, i] = sineθ * sineφ
		v[3, i] = cosineθ
	end
	return v
end

function project_onto_spherical_grid(
	φ::AbstractVector{<: Real},
	θ::AbstractVector{<: Real},
	v::AbstractMatrix{<: Real} # Needs to be normalised
)
	@assert length(φ) == length(θ)

	nt = Threads.nthreads()
	w = Matrix{Float64}(undef, 3, nt)
	c = Matrix{Float64}(undef, size(v, 2), nt)
	idx = Vector{Int}(undef, length(φ))
	@views @inbounds Threads.@threads :static for i in eachindex(φ)
		tid = Threads.threadid()
		sineφ, cosineφ = sincos(φ[i])
		sineθ, cosineθ = sincos(θ[i])
		w[1, tid] = sineθ * cosineφ
		w[2, tid] = sineθ * sineφ
		w[3, tid] = cosineθ
		mul!(c[:, tid]', w[:, tid]', v)
		idx[i] = argmax(c[:, tid])
	end

	return idx
end

