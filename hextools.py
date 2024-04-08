import numpy as np
"""
Different tools to convert between hex coordinates systems, hex selection
and other. To work with flat top hex version

     ___
 ___/ N \___ 
/NW \___/NE \   Flat Top neighbors
\___/ x \___/  
/SW \___/SE \ 
\___/ S \___/
    \___/
"""


SQRT3 = np.sqrt(3)

# coords transformations functions

def evenq2cube(evenq, fractional = False):
    """
    Convert evenq to cube coordinates
    :param evenq: A coordinate tuple in evenq form (q,r)
    :return: `evenq` in cube form (q,r,s).
    """
    if isinstance(evenq, np.ndarray):
        q = evenq[:,0]
        if not fractional:
            r = evenq[:,1] - ((q + (q.astype(int) & 1))/2).astype(int)
        else:
            r = evenq[:,1] - ((q + (q.astype(int) & 1))/2)
        s = -q -r
        return np.vstack((q,r,s)).T
 
    else:          
        if not len(evenq) == 2:
            raise ValueError('evenq shall be a (col, row) pair')
        
        if not all(isinstance(coord, int) for coord in evenq) \
                                                and not fractional:
            raise ValueError('evenq elements are index (shall be int).')
        col, row = evenq
        q = col
        r = row - ((col + (int(col) & 1)) / 2)
        if not fractional:
            r = int(r)
        s = -q -r
        return (q, r, s)

def cube2evenq(cube):
    """
    Convert cube to evenq coordinates
    :param cube: A coordinate tuple in evenq form (q,r,s)
    :return: `cube` in evenq form (q,r).
    """
    if isinstance(cube, np.ndarray):
        q = cube[:,0]
        r = cube[:,1]
        s = cube[:,2]
        if np.any((q+r+s)) != 0:
            raise ValueError('some cube coordinates are inconsistent. Pleas check')
        row = r + ((q + (q & 1))/2).astype(int)
        return np.vstack((q,row)).T    
        
    else:   
        if not len(cube) == 3:
            raise ValueError('cube shall be an (q, r, s) triplet')
        if not all(isinstance(coord, int) for coord in cube):
            raise ValueError('cube elements are index (shall be int).')
        q, r, s = cube
        if q + r + s !=0:
            raise ValueError('Inconsistent cube coordinates')
        col = q
        row = r + int((q + (q & 1)) /2)
        return (col, row)

def evenq2axial(evenq):
    """
    Convert evenq to axial coordinates
    :param evenq: A coordinate tuple in evenq form (q,r)
    :return: `evenq` in axial form (q,r).
    """
    if isinstance(evenq, np.ndarray):
        q = evenq[:,0]
        r = evenq[:,1] - ((q + (q & 1))/2).astype(int)
        
        return np.vstack((q,r).T)
    
    else:      
        if not len(evenq) == 2:
            raise ValueError('evenq shall be a (col, row) pair')
        
        if not all(isinstance(coord, int) for coord in evenq):
            raise ValueError('evenq elements are index (shall be int).')
        col, row = evenq
        q = int(col)
        r = int(row - (col + (col&1)) / 2)
        return (q, r)

def axial2evenq(axial):
    """
    Convert axial to evenq coordinates
    :param axial: A coordinate tuple in axial form (q,r)
    :return: `axial` in evenq form (q,r).
    """
    if isinstance(axial, np.ndarray):
        q = axial[:,0]
        r = axial[:,1] + ((q + (q & 1))/2).astype(int)
        
        return np.vstack((q,r).T)
    
    else:
        if not len(axial) == 2:
            raise ValueError('axial shall be a (q, r) pair')
        
        if not all(isinstance(coord, int) for coord in axial):
            raise ValueError('axial elements are index (shall be int).')
        
        q, r = axial
        col = int(q)
        row = int(r + (q + (q&1)) / 2)
        return (col, row)

def cube2axial(cube):
    """
    Convert cube to axial coordinates
    :param cube: A coordinate tuple in axial form (q, r, s)
    :return: `axial` in evenq form (q, r).
    """
    if isinstance(cube, np.ndarray):
        return np.vstack((cube[:, 0], cube[:, 1])).T
    
    else:
            
        if not len(cube) == 3:
            raise ValueError('cube shall be an (q, r, s) triplet')
        if not all(isinstance(coord, int) for coord in cube):
            raise ValueError('cube elements are index (shall be int).')
        q,r,_ = cube
        return (q, r)

def axial2cube(axial):
    """
    Convert axial to cube coordinates
    :param axial: A coordinate tuple in axial form (q,r)
    :return: `cube` in evenq form (q, r, s).
    """
    if isinstance(axial, np.ndarray):
        x = axial[:, 0]
        z = axial[:, 1]
        y = -x - z
        cube_coords = np.vstack((x, y, z)).T
        return cube_coords
            
    else:
            
        if not len(axial) == 2:
            raise ValueError('axial shall be a (q, r) pair')
        
        if not all(isinstance(coord, int) for coord in axial):
            raise ValueError('axial elements are index (shall be int).')
        q,r = axial
        s = -q -r
        return (q,r,s)

def evenq2pixel(evenq, offset, hexsize = 1, angle=0):
    """
    Convert evenq to pixel coordinates (position of hex center)
    :param evenq: A coordinate tuple in axial form (q,r)
    :param offset: coordinates of the upper left bbox corner (x0,y0)
    :param angle: rotatition angle of hex from flat top position
    :param hexsize: hez radius
    :return: `evenq` in pixel form (x,y).
    """
    x0 = offset[0]
    y0 = offset[1]
    if isinstance(evenq, np.ndarray):
        q = evenq[:,0]
        r = evenq[:,1]
        x = (0.5 + 3/2*q)*hexsize
        y = (SQRT3*(0.5*((q + 1) &1) + r))*hexsize
        x = x + x0
        y = y0 - y
        pixel = np.vstack((x,y)).T
        
        if angle !=0:
            return rotate_pts(pixel,angle)
        else:
            return pixel
    else:
        q,r = evenq
        x = (0.5 + 3/2*q)*hexsize
        y = (SQRT3*(0.5*((q + 1) &1) + r))*hexsize
        x = x + x0
        y = y0 - y
        if angle !=0:
            return rotate_pts((x,y), angle)
        else:
            return (x,y)

def pixel_to_evenq(pixel, offset, hexsize = 1, angle=0):
    x0 = offset[0]
    y0 = offset[1]
    #rotar cuando angulo no sea cero
    x,y =pixel
    x = x - x0
    y = y0- y
    q = (2/3)*((x/hexsize) - 0.5)
    r = (y/(SQRT3*hexsize)) - (1/2)*(int(q + 1)&1)
    ## rounding to neareast hex
    fracq, fracr, fracs = evenq2cube((q,r), fractional=True)
    cq, cr, cs = np.round((fracq, fracr, fracs))
    q_diff = np.abs(cq - fracq)
    r_diff = np.abs(cr - fracr)
    s_diff = np.abs(cs - fracs)
    if q_diff > r_diff and q_diff > s_diff:
        cq = -cr - cs
    elif r_diff > s_diff:
        cr = -cq -cs
    else:
        cs = -cq -cr
    eq, er = cube2evenq((int(cq),int(cr),int(cs)))
    return eq, er

def rotate_pts(points, angle):
    """
    Rotate points according to given angle
    :param points: points in 2d in format (x,y)
    :param angle: angle rotation in rad
    :return: `points` rotated according angle.
    """
    c = np.cos(angle)
    s = np.sin(angle)
    if isinstance(points, np.ndarray):
        rotation_matrix = np.array([[c, -s],
                                    [s,c]])
        return points @ rotation_matrix
    else:
        x = points[0]*c + points[1]*s
        y = points[0]*-s + points[1]*c
        return x,y


## Hex selection
# These are the vectors for moving from any hex to one of its neighbors
# flat top
# These are the vectors for moving from any hex to one of its neighbors
# flat top and cube coordinates system
S = np.array((0, +1, -1))
SE = np.array((+1, 0, -1))
SW = np.array((-1, +1, 0))
N = np.array((0, -1, +1))
NW = np.array((-1, 0, +1))
NE = np.array((+1, -1, 0))

ALL_DIRECTIONS_CUBE_COORDS = np.array([SE, NE, N, NW, SW, S])
# same coordinates order
ALL_DIR_EVENQ_EVENCOLS = np.array([[+1, +1], [+1,  0], [ 0, -1], 
     [-1,  0], [-1, +1], [ 0, +1]])
ALL_DIR_EVENQ_ODDCOLS = np.array([[+1,  0], [+1, -1], [ 0, -1], 
     [-1, -1], [-1,  0], [ 0, +1]])

dir={'SE':0, 'NE':1, 'N':2, 'NW':3, 'SW':4, 'S':5 }

def get_neighbor(currentHex, r_dir, cubecoords =False):
    """
    Return de hex neighbor in the selected direction
    :param currenttHex: coodinates of hex in evenq or cube
    :param r_dir: the relative direction of movement 1 to 5 according to
                  SE, NE, N, NW, SW, S
    :return: hex neighbor in evenq or cube form.
    """
    
    if cubecoords:
        return currentHex + ALL_DIRECTIONS_CUBE_COORDS[r_dir]
    else:
        q,_ = currentHex
        #odd col
        if q & 1:
            return currentHex + ALL_DIR_EVENQ_ODDCOLS[r_dir]
        else:
            return currentHex + ALL_DIR_EVENQ_EVENCOLS[r_dir]
        
        
def get_neighbors(currentHex, cubecoords = False):
    if cubecoords:
        return currentHex + ALL_DIRECTIONS_CUBE_COORDS
    else:
        q,_ = currentHex
        
        #odd col
        if q & 1:
            return currentHex + ALL_DIR_EVENQ_ODDCOLS
        else:
            return currentHex + ALL_DIR_EVENQ_EVENCOLS


coords = [(2,1), (3,1), (0,0), (97,25), (1193487283,32872321)]

def test():
    cubecoords =[]
    for pto in coords:        
        cubecoords.append(evenq2cube(pto))
    
    axialcoords = []
    for pto in coords:
        axialcoords.append(evenq2axial(pto))
    
    ax2cube = []
    for ax in axialcoords:
        ax2cube.append(axial2cube(ax))
    
    res = []
    for c3 in ax2cube:
        res.append(cube2evenq(c3))
    assert np.array_equal(np.array(coords),np.array(res))
    
    qr1 = []
    for  ax  in axialcoords:
       qr1.append(axial2evenq(ax))
    assert np.array_equal(np.array(qr1),np.array(coords))
    
    qr2 = []
    for  ax  in cubecoords:
       qr2.append(cube2evenq(ax))
    assert np.array_equal(np.array(coords),np.array(qr2))
    
    
    #print(coords)
    #print(cubecoords)
    #print(axialcoords)
    #print(ax2cube)
    #print(qr1)
    #print(qr2)
    #print(res)

test()
    
    