# Based on https://gist.github.com/tixxit/252229
# Chan's Convex Hull O(n log h) - Tom Switzer <thomas.switzer@gmail.com>
# Modify to work with points in numpy arrays and python 3.x
# Points of type array([x, y])
# Carlos Sepúlveda nov-2023
# carsepmo@gmail.com

import numpy as np

TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

def cross_product(o, a, b):
    """Cross product between vectors OA and OB."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

def turn(p, q, r):
    """Returns -1, 0, 1 if p, q, r forms a right, straight or left turn."""
    cross = cross_product(p, q, r)
    return bool(cross > 0) - bool(cross < 0)

def _dist(a, b):
    return np.linalg.norm(a - b)

def _keep_left(hull, r):
    """Check if r point should be added to convex hull."""
    while len(hull) > 1 and turn(hull[-2], hull[-1], r) != TURN_LEFT:
        hull = np.delete(hull, -1, 0)
    r = np.expand_dims(r, axis=0)
    hull = np.vstack((hull, r))
    return hull

def _graham_scan(points):
    """Returns points on convex hull of an array of points in CCW order."""
    points = points[np.lexsort((points[:, 1], points[:, 0]))]
    lh = np.empty((0, points.shape[1]), dtype=points.dtype)
    uh = np.empty((0, points.shape[1]), dtype=points.dtype)

    for point in points:
        lh = _keep_left(lh, point)

    for point in points[::-1]:
        uh = _keep_left(uh, point)

    uh = uh[1:-1]
    convex_hull = np.vstack((lh, uh))
    return convex_hull

def _rtangent(hull, p):
    """Return the index of the point in hull that the right tangent line from p to hull touches."""
    l, r = 0, len(hull)
    l_prev = turn(p, hull[0], hull[-1])
    l_next = turn(p, hull[0], hull[(l + 1) % r])
    while l < r:
        c = (l + r) // 2
        c_prev = turn(p, hull[c], hull[(c - 1) % len(hull)])
        c_next = turn(p, hull[c], hull[(c + 1) % len(hull)])
        c_side = turn(p, hull[l], hull[c])

        if c_prev != TURN_RIGHT and c_next != TURN_RIGHT:
            return c
        elif c_side == TURN_LEFT and (l_next == TURN_RIGHT or l_prev == l_next) or c_side == TURN_RIGHT and c_prev == TURN_RIGHT:
            r = c
        else:
            l = c + 1
            l_prev = -c_next
            l_next = turn(p, hull[l], hull[(l + 1) % len(hull)])

    return l

def _min_hull_pt_pair(hulls):
    """Returns the hull, point index pair that is minimal."""
    h, p = 0, 0
    for i in range(len(hulls)):
        j = min(range(len(hulls[i])), key=lambda j: (hulls[i][j][0], hulls[i][j][1]))
        if (hulls[i][j][0], hulls[i][j][1]) < (hulls[h][p][0], hulls[h][p][1]):
            h, p = i, j
    return (h, p)

def _next_hull_pt_pair(hulls, pair):
    """Return the (hull, point) index pair of the next point in the convex hull."""
    p = hulls[pair[0]][pair[1]]
    next = (pair[0], (pair[1] + 1) % len(hulls[pair[0]]))
    for h in (i for i in range(len(hulls)) if i != pair[0]):
        s = _rtangent(hulls[h], p)
        q, r = hulls[next[0]][next[1]], hulls[h][s]
        t = turn(p, q, r)
        if t == TURN_RIGHT or t == TURN_NONE and _dist(p, r) > _dist(p, q):
            next = (h, s)
    return next

def convex_hull(pts):
    """Return the points on the convex hull of pts in CCW order."""
    for m in (1 << (1 << t) for t in range(len(pts))):
        hulls = [_graham_scan(pts[i:i + m]) for i in range(0, len(pts), m)]
        hull = [_min_hull_pt_pair(hulls)]
        for _ in range(m):
            p = _next_hull_pt_pair(hulls, hull[-1])
            if p == hull[0]:
                return np.array([hulls[h][i] for h, i in hull])
            hull.append(p)

def obb(convexhull):
    """
    Return the minimum area oriented bounding box of convex hull
    in CCW order.
    """
    num_points = len(convexhull)
    
    #start minimun obb
    min_area = float('inf')
    obb = None
    
    #iterate trough convexhull points
    for i in range(num_points):
        j = (i + 1) % num_points
        
        side_vec = convexhull[j] - convexhull[i]
        side_lenght = np.linalg.norm(side_vec)
        
        #normalize vector direction
        side_direction = side_vec/side_lenght
        
        #Get perpendicular vector direction
        perpen_direction = np.array([-side_direction[1], side_direction[0]])
        
        #Project all points of convex onto current side
        projected_points = np.dot(convexhull -convexhull[i], side_direction)
        
        #Get coordinares (min and max)
        min_projection = np.min(projected_points)
        max_projection = np.max(projected_points)
        
        #Get OBB area
        width = max_projection - min_projection
        height = side_lenght
        
        area = width * height
        
        #update min area if we found it
        if area < min_area:
            min_area = area
            obb = {
                'center': convex_hull[i] + 0.5 * side_vec,
                'width': width,
                'height': height,
                'angle': np.arctan2(side_direction[1], side_direction[0]),
                'area': area
            }
    return obb
        
        
        
    