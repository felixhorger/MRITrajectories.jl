
module MRITrajectories

	using Random
	using LinearAlgebra
	import Base: rand
	import StatsBase

	include("cartesian.jl")
	include("golden_angle.jl")
	include("radial.jl")
	include("spiral.jl")
	include("extra.jl")
	include("utils.jl")

	# TODO: Not sure this is the best place for it
	"""
			dwelltime2bandwidthppx(δt::Real, num_columns::Integer)

		Bandwidth per pixel in units of 1 / (1000 * [δt] * px)
	"""
	dwelltime2bandwidthppx(δt::Real, num_columns::Integer) = 1000. / (δt * num_columns)

end

