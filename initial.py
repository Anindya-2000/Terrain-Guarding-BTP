import random, math
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from sympy import Point, Line, Segment

MAXLEN = 20
INF = 100000000
OFFSET = 2

# class Point:
#     def __init__(self, pts):
#         self.x = pts[0]
#         self.y = pts[1]
  
# # Given three collinear points p, q, r, the function checks if 
# # point q lies on line segment 'pr' 
# def onSegment(p, q, r):
#     if ( (q.x <= max(p.x, r.x)) and (q.x >= min(p.x, r.x)) and 
#            (q.y <= max(p.y, r.y)) and (q.y >= min(p.y, r.y))):
#         return True
#     return False
  
# def orientation(p, q, r):
#     # to find the orientation of an ordered triplet (p,q,r)
#     # function returns the following values:
#     # 0 : Collinear points
#     # 1 : Clockwise points
#     # 2 : Counterclockwise
      
#     # See https://www.geeksforgeeks.org/orientation-3-ordered-points/amp/ 
#     # for details of below formula. 
      
#     val = (float(q.y - p.y) * (r.x - q.x)) - (float(q.x - p.x) * (r.y - q.y))
#     if (val > 0):
          
#         # Clockwise orientation
#         return 1
#     elif (val < 0):
          
#         # Counterclockwise orientation
#         return 2
#     else:
          
#         # Collinear orientation
#         return 0
  
# # The main function that returns true if 
# # the line segment 'p1q1' and 'p2q2' intersect.
# def doIntersect(p1,q1,p2,q2):
      
#     # Find the 4 orientations required for 
#     # the general and special cases
#     o1 = orientation(p1, q1, p2)
#     o2 = orientation(p1, q1, q2)
#     o3 = orientation(p2, q2, p1)
#     o4 = orientation(p2, q2, q1)
  
#     # General case
#     if ((o1 != o2) and (o3 != o4)):
#         return True
  
#     # Special Cases
  
#     # p1 , q1 and p2 are collinear and p2 lies on segment p1q1
#     if ((o1 == 0) and onSegment(p1, p2, q1)):
#         return True
  
#     # p1 , q1 and q2 are collinear and q2 lies on segment p1q1
#     if ((o2 == 0) and onSegment(p1, q2, q1)):
#         return True
  
#     # p2 , q2 and p1 are collinear and p1 lies on segment p2q2
#     if ((o3 == 0) and onSegment(p2, p1, q2)):
#         return True
  
#     # p2 , q2 and q1 are collinear and q1 lies on segment p2q2
#     if ((o4 == 0) and onSegment(p2, q1, q2)):
#         return True
  
#     # If none of the cases
#     return False

def toRadian(angle):
    return angle * math.pi / 180

def generateTerrainPoints(A0, A1, A2, numTerrains, numPts):
    A0 = toRadian(A0)
    A1 = toRadian(A1)
    A2 = toRadian(A2)
    terrain = []
    dev = A0
    valid = True
    for i in range(numTerrains):
        if dev <= 0:
            valid = False
            break

        subTerrain = []

        if len(terrain) == 0: subTerrain.append((0, 0))
        else: subTerrain.append(terrain[-1][-1])
        
        for j in range(numPts[i] + 1):
            l = int(MAXLEN * random.random())
            if i == numTerrains - 2 and j == numPts[i]:
                l = MAXLEN
            x = subTerrain[-1][0] + l * math.sin(dev)
            y = subTerrain[-1][1] - l * math.cos(dev)
            subTerrain.append((x, y))
            if j < numPts[i]:
                dev += math.pi - A1

            if dev >= math.pi:
                valid = False
                break

        terrain.append(subTerrain)
        if not valid:
            break
        dev += math.pi - A2

    return terrain

def bruteForce(T):
    maxX, maxY, _ = maxMinCoordinates(T)
    maxY += 1000
    maxX /= 2

    l = len(T)
    sol = []
    solLen = INF

    def checkIfVisible(p1, p2, start, end):
        valid = True
        for t in range(start, end + 1):
            if not valid: break
            for ind in range(len(T[t]) - 1):
                if not valid: break
                if doIntersect(p1, p2, Point(T[t][ind]), Point(T[t][ind + 1])):
                    valid = False
                    break
        return valid

    for msk in range(1, 1 << (l + 1)):
        candidate = [i for i in range(l + 1) if (msk >> i) & 1]
        valid = True
        for subterrain in range(l):
            if not valid: break
            if subterrain in candidate or subterrain + 1 in candidate: continue
            left = right = [-1 for _ in range(len(candidate))]
            rleft = rright = [False for _ in range(len(candidate))]
            for i in range(len(candidate)):
                g = candidate[i]
                g_pt = T[g][0] if g < l else T[l - 1][-1]
                p1 = Point(g_pt)
                if g < subterrain:
                    for ind in range(1, len(T[subterrain])):
                        pt = T[subterrain][ind]
                        p2 = Point(pt)
                        valid = checkIfVisible(p1, p2, 0, subterrain)
                        if valid:
                            left[i] = ind
                            break
                    if checkIfVisible(p1, Point(T[subterrain][0]), 0, subterrain):
                        rleft[i] = True
                else:
                    for ind_ in range(1, len(T[subterrain])):
                        ind = len(T[subterrain]) - ind_ - 1
                        pt = T[subterrain][ind]
                        p2 = Point(pt)
                        valid = checkIfVisible(p1, p2, subterrain, l - 1)
                        if valid:
                            right[i] = ind
                            break
                    if checkIfVisible(p1, Point(T[subterrain][-1]), subterrain, l - 1):
                        rright[i] = True
            valid2 = False
            for g in range(len(candidate)):
                if (left[g] == 1 and rleft[g]) or (right[g] == len(T[subterrain]) - 2 and rright[g]):
                    valid2 = True
            if valid2: continue
            for i in range(len(candidate)):
                if valid2: break
                for j in range(i + 1, len(candidate)):
                    if left[i] >= 0 and right[j] >= 0 and left[i] <= right[j]:
                        valid2 = True
                        break
            
            if not valid2:
                valid = False

        if valid:
            if len(candidate) < solLen:
                solLen = len(candidate)
                sol = candidate
           
    for i in range(len(sol)):
        if sol[i] == l:
            sol[i] = T[l - 1][-1]
        else:
            sol[i] = T[sol[i]][0]
    return sol

def maxMinCoordinates(T):
    maxX = 0
    minY = INF
    maxY = -INF
    for i in range(len(T)):
        for j in range(len(T[i])):
            maxX = max(maxX, T[i][j][0])
            maxY = max(maxY, T[i][j][1])
            minY = min(minY, T[i][j][1])
    return [maxX, maxY, minY]

def plotTerrain(T, sol):
    verts = []
    maxX, maxY, minY = maxMinCoordinates(T)
    for i in range(len(T)):
        for j in range(len(T[i])):
            if i:
                if j:
                    verts.append(T[i][j])
            else:
                verts.append(T[i][j])
    
    codes = [Path.MOVETO if not i else Path.LINETO for i in range(len(verts))]
    path = Path(verts, codes)

    fig, ax = plt.subplots()
    patch = patches.PathPatch(path, facecolor='none', lw=2)
    ax.add_patch(patch)

    xs, ys = zip(*verts)
    ax.plot(xs, ys, 'bo', lw=2, ms=5)
    xs, ys = zip(*sol)
    ax.plot(xs, ys, 'ro', lw=2, ms=5)

    ax.set_xlim(-OFFSET, maxX + OFFSET)
    ax.set_ylim(minY - OFFSET, maxY + OFFSET)
    plt.show()

if __name__ == '__main__':
    T = generateTerrainPoints(20, 150, 271, 3, [3, 5, 5])
    # sol = bruteForce(T)
    # plotTerrain(T, sol)
    p1, p2, p3, p4 = Point(0, 0), Point(1, 1), Point(0, 0), Point(0, 1.333)
    l1 = Line(p1, p2)
    s1 = Segment(p3, p4)
    
    
    # using intersection() method
    showIntersection = l1.intersection(s1)
    
    print(showIntersection)
