import numpy as np
import hexy as hx
from shapely.geometry import Polygon, Point
from shapely.affinity import rotate

## "even-q" offset coordiantes

class HexTile():
   
    def __init__(self, offsetCoords):
        if not len(offsetCoords) == 2:
            raise ValueError('offsetCoord shall be an (col, row) pair')
        
        if not all(isinstance(coord, int) for coord in offsetCoords):
            raise ValueError('Elements of offsetCoord shall be int, due they are index ')
        
        self.offsetCoords = offsetCoords
        self.cube_coordinates = np.array(hx.qoffset_to_cube(offsetCoords))


