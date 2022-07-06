
"""
	Assumes the samples are symmetrically arranged around k = 0,
	i.e. if num_r is even, no sample hits k = 0, if num_r is odd, the centre one does.

"""
function radial_spokes(φ::AbstractVector{<: Real}, num_r::Integer)
	r = reshape(range(-π, π; length=num_r), 1, num_r)
	range(-π, π * (1 - 2*mod(num_r+1, 2) / num_r); length=num_r)
	# Allocate memory
	k = Array{Float64, 3}(undef, 2, num_r, length(φ))
	sinecosine = Vector{Float64}(undef, 2)
	# Iterate φ
	for spoke = 1:length(φ)
		sinecosine .= sincos(φ[spoke])
		k[:, :, spoke] = @. r * sinecosine
	end
	return k
end

