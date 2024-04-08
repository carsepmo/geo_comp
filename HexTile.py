import numpy as np
import hextools as hxt
from shapely.geometry import Polygon, Point
from shapely.affinity import rotate

## "even-q" offset coordiantes

class HexTile():
   
    def __init__(self, offsetCoords, position):
        if not len(offsetCoords) == 2:
            raise ValueError('offsetCoord shall be an (col, row) pair')
        
        if not all(isinstance(coord, int) for coord in offsetCoords):
            raise ValueError('Elements of offsetCoord shall be int, due they are index ')
        
        self.evenq = offsetCoords
        self.cube = hxt.evenq2cube(offsetCoords)
        self.position = position
        self.connectedTo = {}
