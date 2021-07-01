# %%
import geopandas as gpd
from scipy.spatial import ConvexHull
from shapely.geometry import mapping, LineString
import numpy as np
import math
import matplotlib.pyplot as plt


def get_contour(in_file='./data/Oahu_10m_nonoaa_contours.shp', contour=-2.0):
    contours = gpd.read_file(in_file)
    contours['shape_length'] = contours['geometry'].length
    max_lengths = contours.groupby('ELEV')['shape_length'].idxmax()
    contours = contours.set_index(['ID']).loc[max_lengths]
    contours = contours[contours['ELEV'] == contour]
    return contours['geometry'].values[0]


def get_convex_points(contour_geom, plot=False):
    points = contour_geom.coords
    hull = ConvexHull(points)
    if plot:
        points = mapping(contour_geom)['coordinates']
        plt.plot([points[i][0] for i in hull.vertices], [points[i][1] for i in hull.vertices], 'r--', lw=2)
        plt.show()
    return hull.vertices


def slice_shape_by_idxs(shape, slice_start, slice_end):
    points = shape.coords
    if slice_end < slice_start:
        points = list(zip(*points.xy))
        return points[slice_start: ] + points[: slice_end+1]
    return points[slice_start: slice_end+1]


def get_contour_convex_slices(contour_geom, convex_idxs):
    convex_idxs = sorted(convex_idxs)
    convex_idxs = np.append(convex_idxs, convex_idxs[0])
    buff = []
    for i in range(len(convex_idxs)-1):
        slice_start = convex_idxs[i]
        slice_end = convex_idxs[i+1]
        buff.append(slice_shape_by_idxs(contour_geom, slice_start, slice_end))
    return buff


def pt_to_pt_distance(pt, other):
    return math.sqrt(math.pow(pt[0]-other[0], 2) + math.pow(pt[1]-other[1], 2))


def calculate_angle(pt, other):
    rads = math.atan2((other[0] - pt[0]), (other[1] - pt[1]))
    rads = rads*180/math.pi
    return -rads + 180


def edge_to_edge_intersection(edge1, edge2, tol=0.001):
    xdiff = (edge1[0][0] - edge1[1][0], edge2[0][0] - edge2[1][0])
    ydiff = (edge1[0][1] - edge1[1][1], edge2[0][1] - edge2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       return False

    d = (det(*edge1), det(*edge2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    xmin = min(edge2[0][0], edge2[1][0])
    xmax = max(edge2[0][0], edge2[1][0])
    ymin = min(edge2[0][1], edge2[1][1])
    ymax = max(edge2[0][1], edge2[1][1])
    if x >= xmin-tol and x <= xmax+tol and y >= ymin-tol and y <= ymax+tol:
        return True
    return False


def get_next_point(p0, points, index, distance=2000, plot=False):

    def filter_by_distance(i):
        return pt_to_pt_distance(p0, i) < distance

    def filter_by_intersection(i):
        i_idx = points.index(i)
        if i_idx == index +1:
            return True
        intermediate_edges = [(points[j], points[j+1]) for j in range(index+1, i_idx-1)]
        line_to_test = (p0, i)
        intersections = [edge_to_edge_intersection(line_to_test, j) for j in intermediate_edges]
        if len(intersections) == 0:
            return True
        return not max(intersections)

    def filter_by_angle(i, angle0, angle1):
        angle = calculate_angle(p0, i)
        if angle1 < angle0:
            if angle < angle1 or angle >= angle0:
                return True
            else:
                return False
        if angle >= angle0 and angle < angle1:
            return True
        else:
            return False

    def get_index_of_max_angle():
        angle0 = calculate_angle(p0, points[index+1])
        if points.index(p0) == 0:
            angle1 = calculate_angle(p0, points[-1])
        else:
            angle1 = calculate_angle(p0, points[index-1])
        filtered_points = [i for i in points[index+1:] if filter_by_distance(i) and filter_by_angle(i, angle0, angle1) and filter_by_intersection(i)]
        # filtered_points = [i for i in filtered_points if filter_by_intersection(i)]
        if len(filtered_points) == 0:
            return index+1
        return points.index(filtered_points[-1])

    return get_index_of_max_angle()


def simplify_single_slice(points):
    buff = []
    i = 0
    p0 = points[i]
    buff.append(p0)
    while i < len(points)-1:
        print(i)
        next_idx = get_next_point(p0, points, i)
        p0 = points[next_idx]
        if next_idx == i:
            break
        i = next_idx
        buff.append(p0)
    buff.append(points[-1])
    return buff


def simplify_convex_slices(slices):
    buff = []
    print('Slices: {}'.format(len(slices)))
    slice_number = 0
    for points in slices:
        print('Slice {}'.format(slice_number))
        buff.append(simplify_single_slice(points))
        slice_number += 1
    return np.concatenate(buff)


def test_is_clockwise(points):
    def test_index(i):
        x = points[i+1][0]-points[i][0]
        y = points[i+1][1]+points[1][0]
        return x*y
            
    idxs = range(len(points)-1)
    clockwise_test = list(map(lambda i: test_index(i), idxs))
    return sum(clockwise_test) > 0
    

contour_2m = get_contour()
assert(test_is_clockwise(contour_2m.coords))
coords = list(zip(*contour_2m.coords.xy))
convex_idxs = get_convex_points(contour_2m)
slices = get_contour_convex_slices(contour_2m, convex_idxs)
points = simplify_convex_slices(slices)

with open('output.txt', 'w') as out:
    out.write(LineString(points).wkt)

# %%
plt.plot([i[0] for i in coords], [i[1] for i in coords], 'b', lw=0.1)
plt.plot([i[0] for i in points], [i[1] for i in points], 'r--', lw=0.2)
plt.show()

# %%

u=[0,1]
v=[1,1]
l=math.sqrt((math.pow(u[0], 2) +  math.pow(u[1],2)) * (math.pow(v[0], 2) + math.pow(v[1],2)))
theta = math.acos(u[0]*v[1] - u[1]*v[0] / l)

# %%

print(calculate_angle((0,0), (1,0))) #Right
print(calculate_angle((0,0), (0,1))) #Top
print(calculate_angle((0,0), (-1,0))) #Left
print(calculate_angle((0,0), (0,-1))) #Bottom
# %%
