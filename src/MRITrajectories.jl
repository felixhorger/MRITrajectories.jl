
module MRITrajectories

	using Random
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
		bandwidth per pixel in units of 1000 / [δt]
	"""
	dwelltime2bandwidthppx(δt::Real, num_columns::Integer) = 1000. / (δt * num_columns)

end

