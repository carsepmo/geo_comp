# Based on https://gist.github.com/tixxit/252229
# Chan's Convex Hull O(n log h) - Tom Switzer <thomas.switzer@gmail.com>
# Modify to work with points in numpy arrays and python 3.x
# Points of type array([x, y])
# Carlos Sepúlveda nov-2023
# carsepmo@gmail.com
"""
agregar np.ascontiguosarray()
"""

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

def _get_bbox_vertices(pts, angle):
    mean = np.float32([pts[:, 0].mean(), pts[:, 0].mean()])
    c, s = np.cos(angle[0]), np.sin(angle[0])
    #Rotation matrix
    R = np.float32([c, -s, s, c]).reshape(2, 2)
    pts = (pts.astype(np.float32) - mean) @ R
    x0, y0 = pts[:, 0].min(), pts[:, 1].min()
    x1, y1 = pts[:, 0].max(), pts[:, 1].max()
    corners = np.float32([x0, y0, x0, y1, x1, y1, x1, y0])
    corners = corners.reshape(-1, 2) @ R.T + mean
    return corners

def _compute_area(pts, caliper_angles):
    """Uses fact that inv(R) = R.T"""
    c = np.cos(caliper_angles[0])
    s = np.sin(caliper_angles[0])
    R = np.float32([c, -s, s, c]).reshape(2, 2)
    pts = pts @ R
    x0, y0 = pts[:, 0].min(), pts[:, 1].min()
    x1, y1 = pts[:, 0].max(), pts[:, 1].max()
    return (x1 - x0) * (y1 - y0)


def minimun_obb(convexhull):
    """
    Return the minimum area oriented bounding box of convex hull.
    """
    num_points = len(convexhull)
    #get left and bottom
    i,l = [convexhull[:, i].argmin()for i in range(2)]
    #get right and upper
    k,j = [convexhull[:, i].argmax()for i in range(2)]

    #start minimun obb
    min_area = np.inf
    obb = None
    
    calipers = np.int32([i, j, k, l]) #first BBX
    caliper_angles = np.float32([0.5, 0, -0.5, 1]) * np.pi


    #iterate trough convexhull points
    for i in range(num_points):
      #roll vertices cw
      calipers_advanced = (calipers + 1 ) % num_points 
      
      #Vectors from previous calipers to candidates
      vec = convexhull[calipers_advanced] - convexhull[calipers]
      
      #Find angles of candidate edgelines
      angles = np.arctan2(vec[:,1], vec[:,0]) 
      
      #Find candidate angle deltas
      angle_deltas = caliper_angles -angles
      
      #Select pivot with smaller rotation
      pivot = np.abs(angle_deltas).argmin()
      calipers[pivot] = calipers_advanced[pivot]
      caliper_angles -= angle_deltas[pivot]

      #check if better set of calipers
      area = _compute_area(convexhull[calipers], caliper_angles)
      if area < min_area:
        min_area = area
        best_calipers = calipers.copy()
        best_caliper_angles = caliper_angles.copy()
    
    vertices_obb = _get_bbox_vertices(convexhull[best_calipers], best_caliper_angles)
    obb ={
            'vertices': vertices_obb,
            'calipers': best_calipers,
            'angles' : best_caliper_angles,
            'area' : min_area
          }
    return obb
        
        
     
    