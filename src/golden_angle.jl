
#= Convention on the "phis":
	ϕ is a phase
	φ is a polar angle (in 3D θ will be the polar angle and φ the azimuthal angle)
	Φ is the golden angle
=#

const Φ = 2π / (1 + √5)

"""
	Golden angle incremented angles
"""
function golden_angle_incremented(num::Integer, angles_per_2π::Integer)
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

