import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Polygon
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
            self.__hex_centers = self.__get_hex_centers
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
        
        xx, yy_even = np.meshgrid(x_coords, y_coords_even)
        xx, yy_odd = np.meshgrid(x_coords, y_coords_odd)
        
        yy = np.where(np.arange(n_hex_h) % 2 == 0, yy_even, yy_odd)
        
        centers = np.stack((xx, yy), axis=1).reshape(-1, 2)
        centers[:,0] = centers[:,0] + rotated_corners[:,0].min()
        centers[:,1] = rotated_corners[:,1].max() - centers[:,1]
        hex_rot_centers = self._rotate_pts_to_hor(centers, -self.orientation)
        return  hex_rot_centers
        