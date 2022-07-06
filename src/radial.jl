
"""
	Assumes the centre sample (or first sample of positive frequencies) is at k = 0

"""
function radial_spokes(φ::AbstractVector{<: Real}, num_r::Integer)
	r = reshape(range(-π, π * (num_r - 2*mod(num_r+1, 2)) / num_r; length=num_r), 1, num_r)
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

