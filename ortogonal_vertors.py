import numpy as np

def basis_vectors(point1, point2):
    #возвращает по заданному направлению 2 перевендикулярных ему вектора
    #point2 - start, point1 - end
    point1 = np.array(point1)
    point2 = np.array(point2)

    direction = point1 - point2
    unit_direction = direction / np.linalg.norm(direction)

    def find_orthogonal_vectors(v):
        if v[0] != 0 or v[1] != 0:
            orthogonal_1 = np.array([-v[1], v[0], 0])
        else:
            orthogonal_1 = np.array([0, v[2], -v[1]])

        orthogonal_1 = orthogonal_1 / np.linalg.norm(orthogonal_1)
        orthogonal_2 = np.cross(v, orthogonal_1)
        orthogonal_2 = orthogonal_2 / np.linalg.norm(orthogonal_2)

        return orthogonal_1, orthogonal_2

    orthogonal_vector_1, orthogonal_vector_2 = find_orthogonal_vectors(unit_direction)

    return unit_direction, orthogonal_vector_1, orthogonal_vector_2


if __name__ == '__main__':
    point1 = np.array([0, 0, 0])
    point2 = np.array([0, 0, 0.5])

    res = basis_vectors(point1, point2)

    print("Единичный вектор направления:", res[0])
    print("Ортогональный вектор 1:", res[1])
    print("Ортогональный вектор 2:", res[2])

