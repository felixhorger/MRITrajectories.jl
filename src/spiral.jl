
"""
	n is number of elements in k-space along one dimension
"""
function required_revolutions(n::Integer, interleaves::Integer)
	return 0.5n / interleaves
end

function spiral(r::Function, Φ::Real, n::Integer, interleaves::Integer)
	k = Array{Float64, 3}(undef, 2, n, interleaves)
	angles = range(0, Φ; length=n+1)[1:n]
	for i = 1:interleaves
		Δφ = 2π * (i-1) / interleaves
		for (j, φ) = enumerate(angles)
			ρ = r(φ)
			sine, cosine = sincos(φ + Δφ)
			k[1, j, i] = ρ * cosine
			k[2, j, i] = ρ * sine
		end
	end
	return k
end

