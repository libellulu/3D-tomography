import numpy as np


class Line:

    def __init__(self, x0, y0, z0, x1, y1, z1):
        """Line in 3d space

        Parameters
        ----------
        x0, y0, z0, x1, y1, z1: float
            Coordinates of two distinct points in the line
        """
        self.direction = np.array([x1 - x0, y1 - y0, z1 - z0])
        self.direction = self.direction / np.linalg.norm(self.direction)

        self.anchor = np.array([x0, y0, z0])

    def intersect(self, voxel):

        return intersect(voxel=voxel, line=self)


class LineSegment:

    def __init__(self, x0, x1):
        """Line segment in 3d space

        Parameters
        ----------
        x0, x1: Array-like size 3
            Coordinates of the start and end of the line segment
        """
        self.start = np.array(x0)
        self.end = np.array(x1)
        self.length = np.linalg.norm(self.start - self.end)


class Voxel:

    def __init__(self, x0, y0, z0, x1, y1, z1):
        """Parallelepiped voxel in 3d space, aligned with the coordinate axis

        Parameters
        ----------
        x0, y0, z0, x1, y1, z1: float
            Coordinates of opposite corners of the voxel
        """

        self.start = np.array([np.minimum(x0, x1), np.minimum(y0, y1), np.minimum(z0, z1)])
        self.end = np.array([np.maximum(x0, x1), np.maximum(y0, y1), np.maximum(z0, z1)])

    def intersect(self, line):

        return intersect(voxel=self, line=line)


def intersect(voxel: Voxel, line: Line) -> LineSegment:

    intersection = np.zeros((2, 3))

    index = 0

    for i in range(3):
        if not line.direction[i] == 0.0:
            k = (voxel.start[i] - line.anchor[i]) / line.direction[i]
            a = voxel.start[i]
            b = line.anchor[(i+1) % 3] + k * line.direction[(i+1) % 3]
            c = line.anchor[(i+2) % 3] + k * line.direction[(i+2) % 3]

            if (
                    (voxel.start[(i+1) % 3] <= b <= voxel.end[(i+1) % 3]) and
                    (voxel.start[(i+2) % 3] <= c <= voxel.end[(i+2) % 3])
            ):
                intersection[index, i] = a
                intersection[index, (i+1) % 3] = b
                intersection[index, (i+2) % 3] = c
                index += 1
                if index == 2:
                    break

            k = (voxel.end[i] - line.anchor[i]) / line.direction[i]
            a = voxel.end[i]
            b = line.anchor[(i + 1) % 3] + k * line.direction[(i + 1) % 3]
            c = line.anchor[(i + 2) % 3] + k * line.direction[(i + 2) % 3]

            if (
                    (voxel.start[(i + 1) % 3] <= b <= voxel.end[(i + 1) % 3]) and
                    (voxel.start[(i + 2) % 3] <= c <= voxel.end[(i + 2) % 3])
            ):
                intersection[index, i] = a
                intersection[index, (i + 1) % 3] = b
                intersection[index, (i + 2) % 3] = c
                index += 1
                if index == 2:
                    break

    if index == 1:  # Line is only tangent to voxel
        intersection_segment = LineSegment(np.zeros(3), np.zeros(3))

    else:  # Line either crosses the vessel or does not
        intersection_segment = LineSegment(intersection[0], intersection[1])

    return intersection_segment


if __name__ == "__main__":
    # Testing -------------------------------------------------------------------------------

    # my_line = Line(0.01, 0.0, 0.01, 1, 1.01, 1.01)
    my_line = Line(0, 0, 0, 1, 1, 1)
    my_voxel = Voxel(0, 0, 0, 1, 1, 1)

    print(my_line.intersect(my_voxel).length)
