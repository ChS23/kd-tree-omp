import numpy as np


def generate_points(
        num_points: int,
        dimensions: int,
        filename: str,
):
    points = np.random.rand(num_points, dimensions)
    np.savetxt(filename, points, delimiter=',', fmt='%f')


NUM_POINTS = 1000000
NUM_DIMENSIONS = 10

generate_points(NUM_POINTS, NUM_DIMENSIONS, 'data/points.csv')
