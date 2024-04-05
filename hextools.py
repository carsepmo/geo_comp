import numpy as np
SQRT3 = np.sqrt(3)

# coords transformations functions

def evenq2cube(evenq):
    """
    Convert evenq to cube coordinates
    :param evenq: A coordinate tuple in evenq form (q,r)
    :return: `evenq` in cube form (q,r,s).
    """
    if isinstance(evenq, np.ndarray):
        q = evenq[:,0]
        r = evenq[:,1] - ((q + (q & 1))/2).astype(int)
        s = -q -r
        return np.vstack((q,r,s)).T
 
    else:          
        if not len(evenq) == 2:
            raise ValueError('evenq shall be a (col, row) pair')
        
        if not all(isinstance(coord, int) for coord in evenq):
            raise ValueError('evenq elements are index (shall be int).')
        col, row = evenq
        q = col
        r = row - int((col + (col & 1)) / 2)
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

def pixel_to_evenq(evenq, offset, angle, hexsize = 1):
    x0 = offset[0]
    y0 = offset[1]
    pass

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
    
    