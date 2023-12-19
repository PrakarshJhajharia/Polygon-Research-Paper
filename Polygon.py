#Importing necessary packages
import math
from itertools import combinations
import matplotlib.pyplot as plt


#General Functions
def regular_polygon_vertices(n):
    if n < 3:
        print("A regular polygon must have at least 3 sides.")
        return

    # Angle between consecutive vertices
    angle = 2 * math.pi / n

    # Coordinates of the vertices
    vertices = []
    for i in range(n):
        x = round(math.cos(i * angle), 10)
        y = round(math.sin(i * angle), 10)
        vertices.append((x, y))

    return vertices
def line_equation(point1, point2):
    x1, y1 = point1
    x2, y2 = point2

    # Calculate slope (m) and y-intercept (c)
    if x2 - x1 != 0:
        m = round((y2 - y1) / (x2 - x1), 10)
        c = round(y1 - m * x1, 10)
    else:
        # For vertical lines, slope is infinity, and y-intercept is undefined
        m = math.inf
        c = x1

    return m, c
def intersection_point(line1, line2):
    m1, c1 = line1
    m2, c2 = line2

    # Check if lines are parallel or coincident
    if m1 == m2:
        return None  # Parallel lines, no intersection
    elif m1 == math.inf:
        x = c1
        y = m2 * x + c2
    elif m2 == math.inf:
        x = c2
        y = m1 * x + c1
    else:
        # Calculate intersection point for non-vertical lines
        x = (c2 - c1) / (m1 - m2)
        y = m1 * x + c1

    return round(x, 10), round(y, 10)
def is_point_inside_polygon(point, polygon):
    x, y = point
    count = 0
    for i in range(len(polygon)):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % len(polygon)]
        if ((y1 <= y and y < y2) or (y2 <= y and y < y1)) and (x < x1 + (y - y1) * (x2 - x1) / (y2 - y1)):
            count += 1
    return count % 2 == 1
def find_nearest_adjacent_points(diagonal, intersection_point,vertices,intersection_points):
    # Find the nearest adjacent points above and below the intersection point on the diagonal
    points = list(vertices)+list(intersection_points)
    points_on_line = []
    for z in points:
        if is_point_on_line(z,line_equation(diagonal[0],diagonal[1])):
            points_on_line.append(z)
    a = []
    b = []
    for z in points_on_line:
        distance = math.sqrt((z[0] - intersection_point[0]) ** 2 + (z[1] - intersection_point[1]) ** 2)
        if (a == [] or a[1] > distance) and z[1]>intersection_point[1]:
            a = [z,distance]
        elif (b == [] or b[1] > distance) and z[1]<intersection_point[1]:
            b = [z,distance]
        elif (a == [] or a[1] > distance) and z[0]>intersection_point[0] and round(z[1],7) == round(intersection_point[1],7):
            a = [z, distance]
        elif (b == [] or b[1] > distance) and z[0]<intersection_point[0] and round(z[1],7) == round(intersection_point[1],7):
            b = [z, distance]
    return a[0], b[0]
def is_diagonal_passing_through_points(point1, point2, intersection_points,vertices):
    # Check if any diagonal connecting any 2 vertices passes through both points
    diagonals = list(combinations(vertices, 2))
    for diagonal in diagonals:
        vertex1, vertex2 = diagonal
        line = line_equation(vertex1, vertex2)

        # Check if the line passes through both points
        if is_point_on_line(point1, line) and is_point_on_line(point2, line):

            # Check if no other intersection point lies between these 2 points on the diagonal
            if not is_intersection_between_points(point1, point2, intersection_points,line):
                return True
    return False
def is_point_on_line(point, line):
    x, y = point
    m, c = line
    if m == math.inf:
        if round(x,5) == round(c,5):

            return True
        else:
            return False
    if m == None:
        m = 0
    # Check if the point is on the line
    return round(y, 5) == round(m * x + c, 5)
def is_intersection_between_points(point1, point2, intersection_points,line):
    # Check if any other intersection point lies between the two points on the diagonal
    for intersection_point in intersection_points:
        if is_point_on_line(intersection_point,line):
            # Check if the intersection point is between the two points on the diagonal
            if is_point_between_points(intersection_point, point1, point2):
                return True
    return False
def is_point_between_points(point, start_point, end_point):
    # Check if the point is between the start_point and end_point
    x, y = point
    x1, y1 = start_point
    x2, y2 = end_point

    # Check if the x-coordinate is between the x-coordinates of the start and end points
    x_condition = (round(x1,7) < round(x,7) < round(x2,7)) or (round(x2,7) < round(x,7) < round(x1,7))

    # Check if the y-coordinate is between the y-coordinates of the start and end points
    y_condition = (round(y1,7) < round(y,7) < round(y2,7)) or (round(y2,7) < round(y,7) < round(y1,7))

    if round(x1,7) == round(x,7) == round(x2,7):
        x_condition = True

    if round(y1,7) == round(y,7) == round(y2,7):
        y_condition = True

    return x_condition and y_condition

def plot_lines(vertices, intersection_points, valid_triangles, valid_quad,valid_pent,valid_hex):
    # Plot the vertices
    #plt.scatter(*zip(*vertices), color='red')

    # Plot the lines
    vertex_combinations = list(combinations(vertices, 2))
    for combination in vertex_combinations:
        point1, point2 = combination
        m, c = line_equation(point1, point2)

        if m != math.inf:
            # Not a vertical line
            x_values = [point1[0], point2[0]]
            y_values = [m * x + c for x in x_values]
        else:
            # Vertical line
            x_values = [point1[0], point1[0]]  # Same x value for both points
            y_values = [point1[1], point2[1]]

        plt.plot(x_values, y_values)

    # Plot the intersection points inside the polygon
    #plt.scatter(*zip(*intersection_points), color='blue')

    for triangle in valid_triangles:
        triangle.append(triangle[0])  # Close the loop
        plt.plot(*zip(*triangle), color='green')
        x, y = zip(*triangle)
        plt.fill_between(x, y, color='green', alpha=0.3)

    for triangle in valid_quad:
        triangle.append(triangle[0])  # Close the loop
        x, y = zip(*triangle)
        plt.fill_between(x, y, color='orange', alpha=0.3)

    for triangle in valid_pent:
        triangle.append(triangle[0])  # Close the loop
        plt.plot(*zip(*triangle), color='green')
        x, y = zip(*triangle)
        plt.fill_between(x, y, color='blue', alpha=0.3)

    for triangle in valid_hex:
        triangle.append(triangle[0])  # Close the loop
        plt.plot(*zip(*triangle), color='green')
        x, y = zip(*triangle)
        plt.fill_between(x, y, color='yellow', alpha=0.3)

    # colour = ["Violet","Orange","Red","Pink","Blue","Green","Yellow"]
    # for j in range(6):
    #     plt.scatter(valid_triangles[0][j][0], valid_triangles[0][j][1],color=colour[j])

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.figure(figsize=(150, 100))
    plt.show()

#Triangle
def calculate_valid_triangles(vertices, intersection_points):
    valid_triangles = []
    # Iterate over each intersection point
    for intersection_point in intersection_points:
        # Find the two diagonals that intersect at the current intersection point
        diagonals = []
        for i in range(len(vertices) - 1):
            for j in range(i + 1, len(vertices)):
                line1 = line_equation(vertices[i], intersection_point)
                line2 = line_equation(vertices[j], intersection_point)
                a,b = line1
                c,d = line2
                if round(a,3)==round(c,3) and round(b,3)==round(d,3):
                    diagonals.append([vertices[i], vertices[j]])

        # Find the 4 nearest adjacent points on the diagonals
        adjacent_points = []
        for diagonal in diagonals:
            nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, intersection_point,vertices,intersection_points)
            adjacent_points.append([nearest_up, nearest_down])
        adjacent_pairs = [[adjacent_points[0][0],adjacent_points[1][0]],[adjacent_points[0][1],adjacent_points[1][1]],[adjacent_points[0][0],adjacent_points[1][1]],[adjacent_points[0][1],adjacent_points[1][0]]]

        # Check if a diagonal connecting any 2 vertices passes through both points
        for pair in adjacent_pairs:
            if is_diagonal_passing_through_points(pair[0], pair[1], intersection_points,vertices):
                valid_triangles.append(list(pair) + [intersection_point])
    return valid_triangles
def remove_duplicate_triangle(valid_triangles):
    no_duplicate_valid_triangles = []
    for x in valid_triangles:
        flag = 0
        for y in no_duplicate_valid_triangles:
            if x[0] in y and x[1] in y and x[2] in y:
                flag = 1
        if flag == 0:
            no_duplicate_valid_triangles.append(x)
    return no_duplicate_valid_triangles

#Quadrilateral
def quad_check(pair1, pair2, intersection_points, vertices, intersection_point):
    diagonals = []
    for i in range(len(vertices) - 1):
        for j in range(i + 1, len(vertices)):
            line1 = line_equation(vertices[i], pair1)
            line2 = line_equation(vertices[j], pair1)
            a, b = line1
            c, d = line2
            if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                diagonals.append([vertices[i], vertices[j]])

    # Find the 4 nearest adjacent points on the diagonals
    adjacent_points1 = []
    for diagonal in diagonals:
        nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, pair1, vertices, intersection_points)
        adjacent_points1.append(nearest_up)
        adjacent_points1.append(nearest_down)

    diagonals = []
    for i in range(len(vertices) - 1):
        for j in range(i + 1, len(vertices)):
            line1 = line_equation(vertices[i], pair2)
            line2 = line_equation(vertices[j], pair2)
            a, b = line1
            c, d = line2
            if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                diagonals.append([vertices[i], vertices[j]])

    # Find the 4 nearest adjacent points on the diagonals
    adjacent_points2 = []
    for diagonal in diagonals:
        nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, pair2, vertices, intersection_points)
        adjacent_points2.append(nearest_up)
        adjacent_points2.append(nearest_down)

    for z in adjacent_points1:
        if z in adjacent_points2 and z != intersection_point:
            return z, True
    return 0, False
def calculate_valid_quad(vertices, intersection_points):
    valid_quad = []
    # Iterate over each intersection point
    for intersection_point in intersection_points:
        # Find the two diagonals that intersect at the current intersection point
        diagonals = []
        for i in range(len(vertices) - 1):
            for j in range(i + 1, len(vertices)):
                line1 = line_equation(vertices[i], intersection_point)
                line2 = line_equation(vertices[j], intersection_point)
                a, b = line1
                c, d = line2
                if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                    diagonals.append([vertices[i], vertices[j]])

        # Find the 4 nearest adjacent points on the diagonals
        adjacent_points = []
        for diagonal in diagonals:
            nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, intersection_point, vertices,
                                                                    intersection_points)
            adjacent_points.append([nearest_up, nearest_down])

        # Find all 4 adjacent pairs of the 4 points
        adjacent_pairs = [[adjacent_points[0][0], adjacent_points[1][0]],
                          [adjacent_points[0][1], adjacent_points[1][1]],
                          [adjacent_points[0][0], adjacent_points[1][1]],
                          [adjacent_points[0][1], adjacent_points[1][0]]]


        # Check if a diagonal connecting any 2 vertices passes through both points
        for pair in adjacent_pairs:
            if pair[0] not in vertices and pair[1] not in vertices:
                common_point, check = quad_check(pair[0], pair[1], intersection_points, vertices, intersection_point)
                if check:
                    valid_quad.append([pair[0]] + [common_point] + [pair[1]]+ [intersection_point])
    return valid_quad
def remove_duplicate_quad(valid_quad):
    no_duplicate_valid_quad = []
    for x in valid_quad:
        flag = 0
        for y in no_duplicate_valid_quad:
            if x[0] in y and x[1] in y and x[2] in y and x[3] in y:
                flag = 1
        if flag == 0:
            no_duplicate_valid_quad.append(x)
    return no_duplicate_valid_quad

#Pentagon
def pent_check(pair1, pair2, intersection_points,vertices,intersection_point,adjacent):
    diagonals = []
    for i in range(len(vertices) - 1):
        for j in range(i + 1, len(vertices)):
            line1 = line_equation(vertices[i], pair1)
            line2 = line_equation(vertices[j], pair1)
            a,b = line1
            c,d = line2
            if round(a,3)==round(c,3) and round(b,3)==round(d,3):
                diagonals.append([vertices[i], vertices[j]])

    # Find the 4 nearest adjacent points on the diagonals
    adjacent_points1 = []
    for diagonal in diagonals:
        nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, pair1,vertices,intersection_points)
        adjacent_points1.append(nearest_up)
        adjacent_points1.append(nearest_down)

    diagonals = []
    for i in range(len(vertices) - 1):
        for j in range(i + 1, len(vertices)):
            line1 = line_equation(vertices[i], pair2)
            line2 = line_equation(vertices[j], pair2)
            a, b = line1
            c, d = line2
            if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                diagonals.append([vertices[i], vertices[j]])

    # Find the 4 nearest adjacent points on the diagonals
    adjacent_points2 = []
    for diagonal in diagonals:
        nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, pair2, vertices, intersection_points)
        adjacent_points2.append(nearest_up)
        adjacent_points2.append(nearest_down)

    # for z in adjacent_points1:
    for point1 in adjacent_points1:
        if point1 != intersection_point and point1 not in adjacent and point1 not in adjacent_points2 and point1 not in vertices:
            for point2 in adjacent_points2:
                if point2 != intersection_point and point2 not in adjacent and point2 not in adjacent_points1 and point2 not in vertices:
                    if is_diagonal_passing_through_points(point1, point2, intersection_points,vertices) and not is_diagonal_passing_through_points(point1, pair2, intersection_points,vertices) and not is_diagonal_passing_through_points(point2, pair1, intersection_points,vertices):
                        return point1, point2, True
    return 0, 0,False
def calculate_valid_pent(vertices, intersection_points):
    valid_pent = []
    # Iterate over each intersection point
    for intersection_point in intersection_points:
        # Find the two diagonals that intersect at the current intersection point
        diagonals = []
        for i in range(len(vertices) - 1):
            for j in range(i + 1, len(vertices)):
                line1 = line_equation(vertices[i], intersection_point)
                line2 = line_equation(vertices[j], intersection_point)
                a,b = line1
                c,d = line2
                if round(a,3)==round(c,3) and round(b,3)==round(d,3):
                    diagonals.append([vertices[i], vertices[j]])


        # Find the 4 nearest adjacent points on the diagonals
        adjacent_points = []
        adjacent = []
        for diagonal in diagonals:
            nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, intersection_point,vertices,intersection_points)
            adjacent_points.append([nearest_up, nearest_down])
            adjacent.append(nearest_up)
            adjacent.append(nearest_down)


        # Find all 4 adjacent pairs of the 4 points
        adjacent_pairs = [[adjacent_points[0][0],adjacent_points[1][0]],[adjacent_points[0][1],adjacent_points[1][1]],[adjacent_points[0][0],adjacent_points[1][1]],[adjacent_points[0][1],adjacent_points[1][0]]]

        # Check if a diagonal connecting any 2 vertices passes through both points
        for pair in adjacent_pairs:
            if pair[0] not in vertices and pair[1] not in vertices:
                common_point1, common_point2, check = pent_check(pair[0], pair[1], intersection_points,vertices,intersection_point,adjacent)
                if check and not is_diagonal_passing_through_points(pair[0], pair[1], intersection_points,vertices):
                    valid_pent.append([pair[0]]+[intersection_point]+[pair[1]] + [common_point2] + [common_point1])

    return valid_pent
def remove_duplicate_pent(valid_pent):
    no_duplicate_valid_pent = []
    for x in valid_pent:
        flag = 0
        for y in no_duplicate_valid_pent:
            if x[0] in y and x[1] in y and x[2] in y and x[3] in y and x[4] in y:
                flag = 1
        if flag == 0:
            no_duplicate_valid_pent.append(x)
    return no_duplicate_valid_pent

#Hexagon
def hex_check(pair1, pair2, intersection_points, vertices, intersection_point, adjacent):
    diagonals = []
    for i in range(len(vertices) - 1):
        for j in range(i + 1, len(vertices)):
            line1 = line_equation(vertices[i], pair1)
            line2 = line_equation(vertices[j], pair1)
            a, b = line1
            c, d = line2
            if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                diagonals.append([vertices[i], vertices[j]])

    # Find the 4 nearest adjacent points on the diagonals
    adjacent_points1 = []
    for diagonal in diagonals:
        nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, pair1, vertices, intersection_points)
        adjacent_points1.append(nearest_up)
        adjacent_points1.append(nearest_down)

    diagonals = []
    for i in range(len(vertices) - 1):
        for j in range(i + 1, len(vertices)):
            line1 = line_equation(vertices[i], pair2)
            line2 = line_equation(vertices[j], pair2)
            a, b = line1
            c, d = line2
            if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                diagonals.append([vertices[i], vertices[j]])

    # Find the 4 nearest adjacent points on the diagonals
    adjacent_points2 = []
    for diagonal in diagonals:
        nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, pair2, vertices, intersection_points)
        adjacent_points2.append(nearest_up)
        adjacent_points2.append(nearest_down)

    # for z in adjacent_points1:
    for point1 in adjacent_points1:
        if point1 != intersection_point and point1 not in adjacent and point1 not in adjacent_points2 and point1 not in vertices:
            for point2 in adjacent_points2:
                if point2 != intersection_point and point2 not in adjacent and point2 not in adjacent_points1 and point2 not in vertices:
                    if not is_diagonal_passing_through_points(point1, point2, intersection_points, vertices) and not is_diagonal_passing_through_points(point1,pair2,intersection_points,vertices) and not is_diagonal_passing_through_points(point2, pair1, intersection_points, vertices) and not is_diagonal_passing_through_points(point1,intersection_point,intersection_points,vertices) and not is_diagonal_passing_through_points(point2,intersection_point,intersection_points,vertices):

                        diagonals = []
                        for i in range(len(vertices) - 1):
                            for j in range(i + 1, len(vertices)):
                                line1 = line_equation(vertices[i], point1)
                                line2 = line_equation(vertices[j], point1)
                                a, b = line1
                                c, d = line2
                                if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                                    diagonals.append([vertices[i], vertices[j]])

                        # Find the 4 nearest adjacent points on the diagonals
                        adjacent_p1 = []
                        for diagonal in diagonals:
                            nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, point1, vertices,
                                                                                    intersection_points)
                            adjacent_p1.append(nearest_up)
                            adjacent_p1.append(nearest_down)

                        diagonals = []
                        for i in range(len(vertices) - 1):
                            for j in range(i + 1, len(vertices)):
                                line1 = line_equation(vertices[i], point2)
                                line2 = line_equation(vertices[j], point2)
                                a, b = line1
                                c, d = line2
                                if round(a, 3) == round(c, 3) and round(b, 3) == round(d, 3):
                                    diagonals.append([vertices[i], vertices[j]])

                        # Find the 4 nearest adjacent points on the diagonals
                        adjacent_p2 = []
                        for diagonal in diagonals:
                            nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, point2, vertices,
                                                                                    intersection_points)
                            adjacent_p2.append(nearest_up)
                            adjacent_p2.append(nearest_down)

                        for point3 in adjacent_p1:
                            if point3 in adjacent_p2 and point3 != intersection_point and point3 not in adjacent and point3 not in adjacent_points1 and point3 not in vertices and point3 not in adjacent_points2:
                                if not is_diagonal_passing_through_points(point3,pair1,intersection_points,vertices) and not is_diagonal_passing_through_points(point3,pair2,intersection_points,vertices) and not is_diagonal_passing_through_points(point3,intersection_point,intersection_points,vertices):
                                    return point1,point2,point3, True

    return 0, 0, 0, False
def calculate_valid_hex(vertices, intersection_points):
    valid_hex = []
    # Iterate over each intersection point
    for intersection_point in intersection_points:
        # Find the two diagonals that intersect at the current intersection point
        diagonals = []
        for i in range(len(vertices) - 1):
            for j in range(i + 1, len(vertices)):
                line1 = line_equation(vertices[i], intersection_point)
                line2 = line_equation(vertices[j], intersection_point)
                a,b = line1
                c,d = line2
                if round(a,3)==round(c,3) and round(b,3)==round(d,3):
                    diagonals.append([vertices[i], vertices[j]])

        # Find the 4 nearest adjacent points on the diagonals
        adjacent_points = []
        adjacent = []
        for diagonal in diagonals:
            nearest_up, nearest_down = find_nearest_adjacent_points(diagonal, intersection_point,vertices,intersection_points)
            adjacent_points.append([nearest_up, nearest_down])
            adjacent.append(nearest_up)
            adjacent.append(nearest_down)

        # Find all 4 adjacent pairs of the 4 points
        adjacent_pairs = [[adjacent_points[0][0],adjacent_points[1][0]],[adjacent_points[0][1],adjacent_points[1][1]],[adjacent_points[0][0],adjacent_points[1][1]],[adjacent_points[0][1],adjacent_points[1][0]]]


        # Check if a diagonal connecting any 2 vertices passes through both points
        for pair in adjacent_pairs:
            if pair[0] not in vertices and pair[1] not in vertices:
                common_point1, common_point2, common_point3, check = hex_check(pair[0], pair[1], intersection_points,vertices,intersection_point,adjacent)
                if check and not is_diagonal_passing_through_points(pair[0], pair[1], intersection_points,vertices):
                    valid_hex.append([pair[0]]+[intersection_point]+[pair[1]] + [common_point2] + [common_point3] + [common_point1])

    return valid_hex
def remove_duplicate_hex(valid_hex):
    no_duplicate_valid_hex = []
    for x in valid_hex:
        flag = 0
        for y in no_duplicate_valid_hex:
            if x[0] in y and x[1] in y and x[2] in y and x[3] in y and x[4] in y and x[5] in y:
                flag = 1
        if flag == 0:
            no_duplicate_valid_hex.append(x)
    return no_duplicate_valid_hex

# Running the code
def calculate(n):

    # Calculate the rounded vertices
    vertices = regular_polygon_vertices(n)

    # Calculate and store the intersection points
    intersection_points = set()
    diagonals_combinations = list(combinations(vertices, 2))
    for i in range(len(diagonals_combinations) - 1):
        for j in range(i + 2, len(diagonals_combinations)):
            line1 = line_equation(*diagonals_combinations[i])
            line2 = line_equation(*diagonals_combinations[j])
            intersection = intersection_point(line1, line2)
            if intersection and is_point_inside_polygon(intersection, vertices):
                intersection_points.add(intersection)
    i =[]
    for x in vertices:
        for z in intersection_points:
            if round(x[0],5)==round(z[0],5) and round(x[1],5)==round(z[1],5):
                i.append(z)
    for z in i:
        intersection_points.remove(z)

    # Calculate and store the valid polygons
    valid_triangles = calculate_valid_triangles(vertices, intersection_points)
    valid_quad = calculate_valid_quad(vertices, intersection_points)
    valid_pent = calculate_valid_pent(vertices, intersection_points)
    valid_hex = calculate_valid_hex(vertices, intersection_points)
    no_duplicate_valid_triangles = remove_duplicate_triangle(valid_triangles)
    no_duplicate_valid_quad = remove_duplicate_quad(valid_quad)
    no_duplicate_valid_pent = remove_duplicate_pent(valid_pent)
    no_duplicate_valid_hex = remove_duplicate_hex(valid_hex)

    # Print the number of valid polygons inside the polygon
    num_valid_triangles = len(no_duplicate_valid_triangles)
    num_valid_quad = len(no_duplicate_valid_quad)
    num_valid_pent = len(no_duplicate_valid_pent)
    num_valid_hex = len(no_duplicate_valid_hex)

    plot_lines(vertices, intersection_points, no_duplicate_valid_triangles,no_duplicate_valid_quad,no_duplicate_valid_pent,no_duplicate_valid_hex)

    if n!=5 and n!=7:
        print(n,"\t\t",num_valid_triangles,"\t\t",num_valid_quad,"\t\t",num_valid_pent,"\t\t",num_valid_hex)
    else:
        print(n,"\t\t",num_valid_triangles,"\t\t",num_valid_quad,"\t\t\t",num_valid_pent,"\t\t\t",num_valid_hex)
def run(start,end):
    print("N \t\t Triangles \t Quad \t\t Pentagons \t Hexagons")

    for i in range(start,end,2):
        calculate(i)

start = 17
end = 21
run(start,end)
