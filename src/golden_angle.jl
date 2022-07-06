
#= Convention on the "phis":
	ϕ is a phase
	φ is a polar angle (in 3D θ will be the polar angle and φ the azimuthal angle)
	Φ is the golden angle
=#

const Φ = 2π / (1 + √5)

"""
	Niche version of golden_angle_incremented, could be used to sample a semi-circle with radial spokes
	to minimise stimulation in gradient spoiled GRE (overlaying spoiler and rephaser).
	however this requires all corrections to work (mainly B0 I presume), otherwise
	imaging artefacts are to be expected.
"""
function golden_angle_incremented(lower::Real, upper::Real, num::Integer, num_angles_in_range::Integer)
	# Angle resolution
	resolution = (upper - lower) / num_angles_in_range
	ξ = Φ / resolution # How many indices to jump if incremented by golden angle
	# Indices of angles in the range
	indices = @. mod(floor(Int64, ξ * (0:num-1)), num_angles_in_range)
	# Angles
	φ = @. indices * resolution + lower
	# One based indexing
	indices .+= 1
	return φ, indices
end

"""
	Golden angle incremented angles
"""
function golden_angle_incremented(num::Integer, angles_per_2π::Integer, )
	# Angle resolution
	resolution = 2π / angles_per_2π
	ξ = Φ / resolution # How many indices to jump if incremented by golden angle
	# Indices of angles in the range
	indices = @. mod(floor(Int64, ξ * (0:num-1)), angles_per_2π)
	# Angles
	φ = @. indices * resolution
	# One based indexing
	indices .+= 1
	return φ, indices
end

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
	(φ, indices) = reshape.((φ, indices), spokes_per_pulse, num_pulses)
	@views for p = 1:num_pulses
		perm = sortperm(φ[:, p])
		φ[:, p] .= φ[perm, p]
		indices[:, p] .= indices[perm, p]
	end
	return φ, indices
end

