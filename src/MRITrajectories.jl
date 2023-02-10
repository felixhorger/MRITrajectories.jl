
module MRITrajectories

	using Random
	import Base: rand
	import StatsBase

	include("golden_angle.jl")
	include("radial.jl")
	include("cartesian.jl")
	include("extra.jl")

	# TODO: Not sure this is the best place for it
	"""
		bandwidth per pixel in units of 1000 / [δt]
	"""
	dwelltime2bandwidthppx(δt::Real, num_columns::Integer) = 1000. / (δt * num_columns)

end

