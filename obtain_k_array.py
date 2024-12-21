import re


def interpolate_points(start, end, num_intersections):
    """Generate array (start, end, num_intersections) between two points."""
    return [
        (
            start[0] + (end[0] - start[0]) * i / num_intersections,
            start[1] + (end[1] - start[1]) * i / num_intersections,
            start[2] + (end[2] - start[2]) * i / num_intersections
        )
        for i in range(num_intersections)
    ]


def parse_file(file_path):
    coordinates = []
    num_intersections = 20  # Default value for intersections

    coordinate_pattern = re.compile(r'(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)')
    intersection_pattern = re.compile(r'^(\d{2})(?:\s+!.*)?\s*$')

    with open(file_path, 'r') as file:
        content = file.readlines()

        # Extract number of intersections
        for line in content:
            match = intersection_pattern.match(line.strip())
            if match:
                num_intersections = int(match.groups()[0])
                break

        # Extract coordinates from the content
        for line in content:
            match = coordinate_pattern.match(line.strip())
            if match:
                x, y, z = map(float, match.groups())
                coordinates.append((x, y, z))

    return coordinates, num_intersections


# Main function to process the file and generate all points
def generate_all_points(file_path):
    coordinates, num_intersections = parse_file(file_path)
    all_points = []

    for i in range(0, len(coordinates) - 1, 2):
        start_point = coordinates[i]
        end_point = coordinates[i + 1]
        all_points.extend(interpolate_points(start_point, end_point, num_intersections))

    return all_points, coordinates

