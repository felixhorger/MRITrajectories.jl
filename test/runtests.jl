
using Revise
import MRITrajectories


split_sampling = [
	[CartesianIndex(1, 1), CartesianIndex(1, 2), CartesianIndex(1, 3)],
	[CartesianIndex(1, 3), CartesianIndex(1, 2), CartesianIndex(2, 2)],
	[CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(1, 2)]
]

MRITrajectories.sort_constantL2distance!(split_sampling, 0.0)

