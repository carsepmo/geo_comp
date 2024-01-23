import numpy as np
import ConvexHull2D as cvxh

class OrientedHexGrid():
    SQRT3 = np.sqrt(3)
    
    def __init__(self, corners, angle, hex_size):
        self.grid_corners = corners
        self.orientation = angle
        self.hex_size = hex_size
        self.__hex_centers = None
        
    def _rotate_pts_to_hor(self, points, angle):
        # Calcula la matriz de rotación inversa
        c = np.cos(angle)
        s = np.sin(angle)
        rotation_matrix_inverse = np.array([[c, -s],
                                            [s, c]])        
        return points @ rotation_matrix_inverse
    
    @property
    def hex_centers(self):
        if self.__hex_centers is None:
            self.__hex_centers = self.__get_hex_centers()
        
        return self.__hex_centers
    
    def __get_hex_centers(self):
        # rotated to horizontal plane
        rotated_corners = \
            self._rotate_pts_to_hor(self.grid_corners, self.orientation)
        
        grid_width = rotated_corners[:,0].max() - \
            rotated_corners[:,0].min()
        
        grid_height = rotated_corners[:,1].max() - \
            rotated_corners[:,1].min()
        
        n_hex_w = np.ceil((grid_width + 0.5 * self.hex_size) 
                          / (1.5 * self.hex_size)).astype(int)
        
        n_hex_h = np.ceil((grid_height * self.SQRT3) / 
                          (3 * self.hex_size) + 0.5).astype(int)
        
        x_coords = self.hex_size * (0.5 + 1.5 * np.arange(n_hex_w))
        
        y_coords_even = self.SQRT3 * self.hex_size * (0.5 + np.arange(n_hex_h))
        y_coords_odd = self.SQRT3 * self.hex_size * np.arange(n_hex_h)
        
        xx, yy_even = np.meshgrid(x_coords, y_coords_even, indexing='ij')
        xx, yy_odd = np.meshgrid(x_coords, y_coords_odd, indexing='ij')
        yy = np.zeros_like(yy_even)
        yy[0::2, :] = yy_even[0::2, :]
        yy[1::2, :] = yy_odd[1::2, :]

        centers = np.dstack([xx, yy])

        centers[:,:,0] = centers[:,:,0] + rotated_corners[:,0].min()
        centers[:,:,1] = rotated_corners[:,1].max() - centers[:,:,1]
        hex_rot_centers = self._rotate_pts_to_hor(centers, -self.orientation)
        return  hex_rot_centers

"""         
######## test
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Polygon
rng = np.random.default_rng()
points = rng.random((300, 2)) *10 
#points = np.array([(2, 6), (4, 9), (5, 7), (9, 10), (7, 6), (13, 6), (11, 1), (5, 4), (7,11), (11,3), (19,7), (12.5,15)])
cnvex_hull = cvxh.ConvexHull2D(points)
hull = cnvex_hull.convex_hull
obb = cnvex_hull.obbx
hex_size = 1.8

hexgrid = OrientedHexGrid(obb['vertices'], obb['angles'][1], hex_size)
centers = hexgrid.hex_centers
fig, ax = plt.subplots()
ax.set_aspect('equal', adjustable='datalim')
for i in np.ndindex(centers.shape[:2]):
    x, y = centers[i]
    hexagon = RegularPolygon((x, y), numVertices = 6,
                                 radius =hex_size, orientation=np.radians(30)+obb['angles'][1],
                             facecolor='none', edgecolor='black'
                             )
    ax.add_patch(hexagon)
    ax.scatter(x, y)
rect_closed = np.vstack((obb['vertices'], obb['vertices'][0]))  # Asegurar que el OBB esté cerrado
rect_patch = Polygon(rect_closed, closed=True, edgecolor='g', facecolor='none', lw=2, label='rect')
ax.add_patch(rect_patch)
points = points[np.lexsort((points[:, 1], points[:, 0]))]

for p in points:
  ax.scatter(p[0],p[1], marker='*' )

plt.xlim(obb['vertices'][:, 0].min() -5, obb['vertices'][:, 0].max() +5)
plt.ylim(obb['vertices'][:, 1].min() -5, obb['vertices'][:, 1].max() +5)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Hexagonal Grid in Rectangle')
plt.show() """