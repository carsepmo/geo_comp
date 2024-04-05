import numpy as np
import ConvexHull2D as cvxh
from shapely.geometry import Polygon
import hexy as hx
import hextools as hxt

class OrientedHexGrid():
    SQRT3 = hxt.SQRT3
    
    def __init__(self, corners, angle, hex_size):
        self.grid_corners = corners
        self.orientation = angle
        self.hex_size = hex_size
        self.__hex_centers = None
        
    def _rotate_pts(self, points, angle):
        return hxt.rotate_pts(points, angle)
    
    @property
    def hex_centers(self):
        if self.__hex_centers is None:
            self.__hex_centers = self.__get_hex_centers()
        
        return self.__hex_centers
    
    def __get_hex_centers(self):
        
        # Oriented boundingbox corners rotated to horizonal plane
        rotated_corners = \
            self._rotate_pts(self.grid_corners, self.orientation)
        
        grid_width = rotated_corners[:,0].max() - \
            rotated_corners[:,0].min()
        
        grid_height = rotated_corners[:,1].max() - \
            rotated_corners[:,1].min()
        
        n_hex_w = np.ceil((grid_width + 0.5 * self.hex_size) 
                          / (1.5 * self.hex_size)).astype(int)
        
        n_hex_h = np.ceil((grid_height * self.SQRT3) / 
                          (3 * self.hex_size) + 0.5).astype(int)
        x0 = rotated_corners[:,0].min()
        y0 = rotated_corners[:,1].max()
        
        ij_matrix = np.transpose(np.indices((n_hex_w, n_hex_h)), axes=(1,2,0))
        centers = hxt.evenq2pixel(ij_matrix.reshape(-1,2), 
                    (x0,y0), self.hex_size).reshape(n_hex_w, n_hex_h,-1)
        
        hex_rot_centers = self._rotate_pts(centers, -self.orientation)
        return  hex_rot_centers

class OHGFromPolygon(OrientedHexGrid):
    def __init__(self, polygon_points, obb, hex_size):
        self.points = polygon_points
        corners = obb['vertices']
        if obb['width'] > obb['height']:
            angle = obb['angles'][0]
        else:
            angle = obb['angles'][1]
        
        super().__init__(corners, angle, hex_size)
        
        self.__hex_centers = None
        self.__valid_hex_tile = None
        self._x0 = None
        self._y0 = None
    
    @property
    def hex_centers(self):
        if self.__hex_centers is None:
            self.__hex_centers, self.__valid_hex_tile = self.__get_hex_centers()        
        return self.__hex_centers
    
    @property
    def valid_hex_tiles(self):
        if self.__valid_hex_tile is None:
            self.__hex_centers, self.__valid_hex_tile = self.__get_hex_centers()
        
        return self.__valid_hex_tile
    
    def _get_hex_tiles_vertices(self, hexcenters, size):
        angles = np.linspace(0, 2*np.pi, 7)[:-1]
        vertices_x = hexcenters[..., 0, None] + size * np.cos(angles)
        vertices_y = hexcenters[..., 1, None] + size * np.sin(angles)

        # Combina las coordenadas para formar los vértices de todos los hexágonos
        # si hex centers es de tamaño MxNx2 devuelve un array
        # de tamaño MxNx2x6
        vertices = np.stack((vertices_x, vertices_y), axis=-2)
        return vertices
    
    def __get_valid_hex(self, hex_vertices, poly_vertices):
        polygon = Polygon(poly_vertices)
        #hex_vertex_reshaped = hex_vertices.reshape((-1,6,2))
        intersecting_hex_tiles = np.zeros((hex_vertices.shape[0],
                                           hex_vertices.shape[1]),
                                           dtype=bool)
        
        #for i in range(hex_vertex_reshaped.shape[0]):
        for i in np.ndindex(hex_vertices.shape[:2]):
            hex_polygon = Polygon(np.column_stack(hex_vertices[i]))
            intersecting_hex_tiles[i] = polygon.intersects(hex_polygon)
        return intersecting_hex_tiles
    
    def __get_hex_centers(self):
        # rotated to horizontal plane
        rotated_corners = \
            self._rotate_pts(self.grid_corners, self.orientation)
        
        
        grid_width = rotated_corners[:,0].max() - \
            rotated_corners[:,0].min()
        
        grid_height = rotated_corners[:,1].max() - \
            rotated_corners[:,1].min()
        
        n_hex_w = np.ceil((grid_width + 0.5 * self.hex_size) 
                          / (1.5 * self.hex_size)).astype(int)
        
        n_hex_h = np.ceil((grid_height * self.SQRT3) / 
                          (3 * self.hex_size) + 0.5).astype(int)
        
        self._x0 = rotated_corners[:,0].min()
        self._y0 = rotated_corners[:,1].max()
        
        ij_matrix = np.transpose(np.indices((n_hex_w, n_hex_h)), axes=(1,2,0))
        centers = hxt.evenq2pixel(ij_matrix.reshape(-1,2), 
                            (self._x0,self._y0), 
                            hex_size).reshape(n_hex_w, n_hex_h,-1)
        
       
        # hex vertex in local coordinate (no rotation)
        hex_vertices = self._get_hex_tiles_vertices(centers, self.hex_size)
        
        # polygon/hull points in local coord
        rot_poly_points = self._rotate_pts(self.points, self.orientation)
        
        self._notrot_hex_cent = centers.copy()
        
        # hex that intersect polygon/hull
        valid_hex_tile = self.__get_valid_hex(hex_vertices, rot_poly_points) 
        
        # hex centers rotated to angle of major OBB axis (global coords) 
        hex_centers = self._rotate_pts(centers, -self.orientation)
        
        return  hex_centers, valid_hex_tile 
    
    def evenq_to_pixel(self, col, row):
        x,y = hx.evenq_to_pixel_local(col, row)
        x = x +self._x0
        y = self._y0 -y
        return self._rotate_pts(np.array([x,y]), -self.orientation)
       
  
class GraphFromHexGrid():
    pass

#E NE NW W SW SE
# []

######## test
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.patches import Polygon as Pol
rng = np.random.default_rng()
points = rng.random((75, 2)) *2000

#points = np.array([(2, 6), (4, 9), (5, 7), (9, 10), (7, 6), (13, 6), (11, 1), (5, 4), (7,11), (11,3), (19,7), (12.5,15)])
cnvex_hull = cvxh.ConvexHull2D(points)
hull = cnvex_hull.convex_hull

obb = cnvex_hull.obbx
hex_size = 32

if obb['width'] > obb['height']:
    angle2 = obb['angles'][0]
else:
    angle2 = obb['angles'][1]


#hexgrid = OrientedHexGrid(obb['vertices'], angle, hex_size)
hexgrid = OHGFromPolygon(hull, obb, hex_size)

centers = hexgrid.hex_centers
fig, ax = plt.subplots()
ax.set_aspect('equal', adjustable='datalim')
print(f'angle:{np.rad2deg(angle2)}')
print(f'width:{np.rad2deg(obb["width"])}')
print(f'height:{np.rad2deg(obb["height"])}')
print(f'GRidShape:{hexgrid._grid_size}')
for i in np.ndindex(centers.shape[:2]):
    facecol = 'gray'
    if not hexgrid.valid_hex_tiles[i]:
        facecol = 'cyan'
    
    if i == (0,0):
        facecol = 'black'
    elif i ==(0,1):
        facecol = 'blue'

    x, y = centers[i]
    hexagon = RegularPolygon((x, y), numVertices = 6,
                                    radius =hex_size, orientation=np.radians(30)+angle2,
                                facecolor=facecol, edgecolor='black'
                                )
    ax.add_patch(hexagon)
    #ax.text(x, y, i)
rect_closed = np.vstack((obb['vertices'], obb['vertices'][0]))  # Asegurar que el OBB esté cerrado
rect_patch = Pol(rect_closed, closed=True, edgecolor='g', facecolor='none', lw=2, label='rect')
ax.add_patch(rect_patch)
cvx_hull = Pol(hull, closed=True, edgecolor='r', facecolor='red', 
                   alpha=0.6, lw=2, label='hull')
ax.add_patch(cvx_hull)
plt.xlim(obb['vertices'][:, 0].min() -5, obb['vertices'][:, 0].max() +5)
plt.ylim(obb['vertices'][:, 1].min() -5, obb['vertices'][:, 1].max() +5)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Hexagonal Grid in Rectangle')
plt.show()