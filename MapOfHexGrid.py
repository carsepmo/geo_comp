from errors import  IncorrectCoordinatesError, HexExistsError, MistmatchError
import numpy as np

def make_key_from_coordinates(indexes):
    """
    Converts indexes to string for hashing
    : param indexes: the indexes of a hex. nx2, n=number of index pairs
    : return: key for hashing
    """
    if isinstance(indexes, tuple) and len(indexes) ==2:
        return str(int(indexes[0])) + ',' + str(int(indexes[1]))
    else:
        return[str(int(index[0])) + ',' + str(int(index[1])) for index in indexes]

class MapOfHexGrid:
    def __init__(self) -> None:
        self._map = {}
    
    def keys(self):
        yield from self._map.key()
    
    def values(self):
        yield from self._map.values()
    
    def items(self):
        yield from self._map.items()
        
    def __len__(self):
        return self._map.__len__()
    
    def __iter__(self):
        yield from self._map
    
    def __setitem__(self, coordinates, hex_tile_object):
        """
        Assigns hex tile objects as values to coordinates as key (evenq coords).
        The nuber of coordinae and hex tile objects should be equal.
        :param coordinates: evenq coordinates of hex tile
        :param: hex_tile: the hex tile objects themselves
        :return: None
        """
        
        if len(coordinates) != len(hex_tile_object):
            raise MistmatchError("Number of coordinates does not match number of hex tile objects.")
        
        keys = make_key_from_coordinates(coordinates)
        for key, hex in zip(keys, hex_tile_object):
            if key in self._map.keys():
                raise HexExistsError('key ' + key + ' already exists.')
            
            self._map[key] = hex
    
    def setitem_direct(self, key, value):
        if key in self._map.keys():
            raise HexExistsError("key " + key + "already exists.")
        self._map[key] = value
    
  
    def overwrite_entries(self, coordinates, hex):
        keys = make_key_from_coordinates(coordinates)
        for key in keys:
            self._map[key] = hex
    
    def __delitem__(self, coordinates):
        if len(coordinates.shape) == 1:
            coordinates = np.array([coordinates])
        
        keys = make_key_from_coordinates(coordinates)
        for key in self.keys():
            del self._map[key]
    
    def __getitem__(self, coordinates):
        """
        Retrives hexes stores at 'coordinates'
        :param coordinate: the locations used as keys for hexes. You can pass more then one coordinate
        :return: list of hexes tiles mapped to using 'coordinates'
        """
        if len(coordinates.shape) == 1:
            coordinates = np.array([coordinates])
        keys = make_key_from_coordinates(coordinates)
        return [self._map.get(k) for k in keys if k in self._map.keys()]
           
        
    
