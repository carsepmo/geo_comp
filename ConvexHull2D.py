import numpy as np
import warnings

class ConvexHull2D():
    """
    comentarios
    """
    
    
    def __init__(self, points: np.ndarray):
        if np.ma.isMaskedArray(points):
            raise ValueError('Input points cannot be a masked array')
        
        points = points[np.lexsort((points[:, 1], points[:, 0]))]
        points = np.ascontiguousarray(points, dtype=np.double)
        try: 
            graham_scan_ch2d = _GrahamScanConvexHull2D(points)
        except GrahamScanConvexHullError as e:
            print(f'error while trying to get convex hull: {e}')
        
        self.convex_hull = graham_scan_ch2d.convx_hull
        self._obbx = None
    
    @property
    def obbx(self):
        if self.convex_hull is None:
            msg = "There is no convex hull calculated"
            warnings.warn(msg, Warning)
            return None
        
        if self._obbx is None:
            self._obbx = self.get_obbx(self.convex_hull)
        
        return self._obbx
    
    def __get_bbox_vertices(self, pts, angle):
        mean = np.mean(pts, axis=0)
        c, s = np.cos(angle[0]), np.sin(angle[0])
        R = np.array([c, -s, s, c]).reshape(2, 2)
        pts_transformed = (pts - mean) @ R
        x0, y0 = pts_transformed.min(axis=0)
        x1, y1 = pts_transformed.max(axis=0)
        corners = np.array([x0, y0, x0, y1, x1, y1, x1, y0]).reshape(-1, 2)
        return corners @ R.T + mean
    
    def __compute_area(self, pts, caliper_angles):
        c, s = np.cos(caliper_angles[0]), np.sin(caliper_angles[0])
        R = np.array([c, -s, s, c]).reshape(2, 2)
        pts_transformed = pts @ R
        x0, y0 = pts_transformed.min(axis=0)
        x1, y1 = pts_transformed.max(axis=0)
        w = x1 - x0
        h = y1 - y0
        return w*h, w, h
    
    def get_obbx(self, convx_hull):
        num_points = len(convx_hull)
        i, l = convx_hull[:, 0].argmin(), convx_hull[:, 1].argmin()
        k, j = convx_hull[:, 0].argmax(), convx_hull[:, 1].argmax()

        min_area = np.inf
        calipers = np.array([i, j, k, l], dtype=np.int32)
        caliper_angles = np.array([0.5, 0, -0.5, 1]) * np.pi

        for _ in range(num_points):
            calipers_advanced = (calipers + 1) % num_points
            vec = convx_hull[calipers_advanced] - convx_hull[calipers]
            angles = np.arctan2(vec[:, 1], vec[:, 0])
            angle_deltas = caliper_angles - angles
            pivot = np.abs(angle_deltas).argmin()
            calipers[pivot] = calipers_advanced[pivot]
            caliper_angles -= angle_deltas[pivot]

            area, width, height = self.__compute_area(convx_hull[calipers], 
                                                      caliper_angles)
            if area < min_area:
                min_area = area
                w, h = width, height
                best_calipers = calipers.copy()
                best_caliper_angles = caliper_angles.copy()

        vertices_obb = self.__get_bbox_vertices(convx_hull[best_calipers], best_caliper_angles)
        return {
            'vertices': vertices_obb,
            'calipers': best_calipers,
            'angles': best_caliper_angles,
            'area': min_area,
            'width': w,
            'height': h
        }

class _GrahamScanConvexHull2D():
    TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)
    
    def __init__(self, points):
        self.convx_hull = self._convex_hull(points)

    def __cross_product(self, o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    
    def __turn(self, p, q, r):
        cross = self.__cross_product(p, q, r)
        return np.sign(cross)
    
    def __dist(self, a, b):
        return np.linalg.norm(a - b)
    
    def __keep_left(self, hull, r):
        while len(hull) > 1 and \
            self.__turn(hull[-2], hull[-1], r)  != self.TURN_LEFT:
                
            hull = np.delete(hull, -1, 0)
            
        r = np.expand_dims(r, axis=0)
        return np.vstack((hull, r))
    
    def _graham_scan(self, points):
        lower_hull = np.empty((0, 2), dtype=points.dtype)
        upper_hull = np.empty((0, 2), dtype=points.dtype)

        for point in points:
            lower_hull = self.__keep_left(lower_hull, point)
        for point in points[::-1]:
            upper_hull = self.__keep_left(upper_hull, point)

        return np.vstack((lower_hull, upper_hull[1:-1]))
    
    def __rtangent(self, hull, p):
        l, r = 0, len(hull)
        while l < r:
            c = (l + r) // 2
            c_prev_turn = self.__turn(p, hull[c], hull[(c - 1) % len(hull)])
            c_next_turn = self.__turn(p, hull[c], hull[(c + 1) % len(hull)])
            c_side = self.__turn(p, hull[l], hull[c])
            

            if c_prev_turn != self.TURN_RIGHT and c_next_turn != self.TURN_RIGHT:
                return c
            elif c_side == self.TURN_LEFT or \
                (c_side == self.TURN_RIGHT and c_prev_turn == self.TURN_RIGHT):
                    
                r = c
            else:
                l = c + 1

        return l
    
    def __min_hull_pt_pair(self, hulls):
        min_points = np.array([np.min(hull, axis=0) for hull in hulls])
        min_point_idx = np.argmin(min_points, axis=0)[0]
        min_hull_idx = np.argmin(hulls[min_point_idx][: , 0])
        return min_point_idx, min_hull_idx
        
    def __next_hull_pt_pair(self, hulls, pair):
        p = hulls[pair[0]][pair[1]]
        next = (pair[0], (pair[1] + 1) % len(hulls[pair[0]]))
        
        for h in (i for i in range(len(hulls)) if i != pair[0]):
            s = self.__rtangent(hulls[h], p)
            q, r = hulls[next[0]][next[1]], hulls[h][s]
            t = self.__turn(p, q, r)
            if t == self.TURN_RIGHT or t == self.TURN_NONE and \
                self.__dist(p, r) > self.__dist(p, q):
                next = (h, s)
        return next
    
    def _convex_hull(self, points):
          for m in (1 << (1 << t) for t in range(len(points))):
            hulls = [self._graham_scan(points[i:i + m]) for i in range(0, len(points), m)]
            hull = [self.__min_hull_pt_pair(hulls)]
            for _ in range(m):
                p = self.__next_hull_pt_pair(hulls, hull[-1])
                # return in CW order to be compatible with obb calipers method
                if p == hull[0]:
                    return np.array([hulls[h][i] for h, i in hull])[::-1]
                hull.append(p)
            
class GrahamScanConvexHullError(Exception):
    def __init__(self, message):
        super().__init__(message)


""" points = np.array([(2, 6), (4, 9), (5, 7), (9, 10), (7, 6), (13, 6), (11, 1), (5, 4), (7,11), (11,3), (19,7), (12.5,15)])
cnvex_hull = ConvexHull2D(points)
hull = cnvex_hull.convex_hull
obb = cnvex_hull.obbx

print(f'convex hull: {hull} \n obb: {obb}') """