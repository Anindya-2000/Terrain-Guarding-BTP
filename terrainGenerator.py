import random, math
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from sympy import Point, Line, Segment
import dotenv
import os, sys

MAXLEN = 20
INF = 100000000
OFFSET = 2
MAX_ITER_PER_TERRAIN = 10

PATH_FIGURES = os.getcwd() + '/Figures'
PATH_DATA = os.getcwd() + '/Data'

dotenv_file = dotenv.find_dotenv()
dotenv.load_dotenv(dotenv_file)
run_seq = os.environ['RUN']

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
            if not l: l += 1
            # if i == numTerrains - 2 and j == numPts[i]:
            #     l = MAXLEN
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

    def FindPOI(line, terrain, pt_to_ignore, msk = 0):
        for i in range(len(terrain) - 1):
            if i == pt_to_ignore: continue
            s = Segment(terrain[i], terrain[i + 1])
            poi = line.intersection(s)
            if len(poi) == 1:
                return poi
        return []

    def checkIfSameSide(coeff, point1, point2 = [0, maxY]):
        a, b, c = coeff
        val1 = a * point1[0] + b * point1[1] + c
        val2 = a * point2[0] + b * point2[1] + c
        if (val1 >= 0 and val2 >= 0) or (val1 < 0 and val2 < 0):
            return True
        return False
  
    for msk in range(1, 1 << (l + 1)):
        candidate = [i for i in range(l + 1) if (msk >> i) & 1]
        valid = True
        for subterrain in range(l):
            if not valid: break
            if subterrain in candidate or subterrain + 1 in candidate: continue
            leftLine = rightLine = [-1 for _ in range(len(candidate))]
            leftPoint = rightPoint = [[] for _ in range(len(candidate))]
            valid2 = False
            for i in range(len(candidate)):
                g = candidate[i]
                g_pt = T[g][0] if g < l else T[l - 1][-1]
                p1 = Point(g_pt)
                if g < subterrain:
                    p2 = Point(T[subterrain][0])
                    leftLine[i] = Line(p1, p2)
                    if Point.is_collinear(p1, p2, Point(T[subterrain][1])) or checkIfSameSide(leftLine[i].coefficients, T[subterrain][1]):
                        valid2 = True
                        break
                    leftPoint[i] = FindPOI(leftLine[i], T[subterrain], 0)

                else:
                    p2 = Point(T[subterrain][-1])
                    rightLine[i] = Line(p1, p2)
                    if Point.is_collinear(p1, p2, Point(T[subterrain][-2])) or checkIfSameSide(rightLine[i].coefficients, T[subterrain][-2]):
                        valid2 = True
                        break
                    rightPoint[i] = FindPOI(rightLine[i], T[subterrain], len(T[subterrain]) - 2)
            
            if valid2: continue
            for i in range(len(candidate)):
                if  valid2: break
                for j in range(i + 1, len(candidate)):
                    if not len(leftPoint[i]) or not len(rightPoint[j]): continue
                    L = Line(leftPoint[i][0], rightPoint[j][0])
                    poi = leftLine[i].intersection(rightLine[j])[0]
                    if not checkIfSameSide(L.coefficients, poi.coordinates):
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

def plotTerrain(T, sol, run_seq):
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
    
    plt.savefig(PATH_FIGURES + '/Figure-' + run_seq, bbox_inches='tight')
    # plt.show()
    

def writeDataToFile(run_seq, T):
    f = open(PATH_DATA + "/Data-" + run_seq + '.txt', "w")
    f.write(f"{len(T)}\n")
    for i in range(len(T)):
        f.write(f"{len(T[i])}\n")
        for pt in T[i]:
            f.write(f"{pt[0]} {pt[1]}\n")
    f.close()

def readTerrainFromFile(run_seq):
    try:
        f = open(PATH_DATA + "/Data-" + run_seq + '.txt', "r")
    except Exception as _:
        raise "File doesn't Exist!!\n"
    numTerrains = int(f.readline())
    T = []
    for i in range(numTerrains):
        T.append([])
        numPts = int(f.readline())
        for _ in range(numPts):
            x, y = [int(val) for val in f.readline().split(' ')]
            T[-1].append((x, y))
    return T

def recursivelyBuildTerrains(cur, numVerts, A0, A1, A2, minVerts, maxVerts, numTerrains):
    if cur == numTerrains:
        for _ in range(MAX_ITER_PER_TERRAIN):
            T = generateTerrainPoints(A0, A1, A2, numTerrains, numVerts)
            if len(T) != numTerrains:
                continue
            sol = bruteForce(T)
            if len(sol) < int((numTerrains + 1)/ 2):
                global run_seq
                writeDataToFile(run_seq, T)
                plotTerrain(T, sol, run_seq)
                run_seq = str(int(run_seq) + 1)
        return

    for v in range(minVerts, maxVerts + 1):
        numVerts.append(v)
        cur += 1
        recursivelyBuildTerrains(cur, numVerts, A0, A1, A2, minVerts, maxVerts, numTerrains)
        numVerts.pop()
        cur -= 1

if __name__ == '__main__':
    try:
        os.mkdir(PATH_FIGURES)
    except Exception as _:
        pass

    try:
        os.mkdir(PATH_DATA)
    except Exception as _:
        pass
    
    
    if len(sys.argv) > 1:
        numTerrains = int(sys.argv[1])
        for A0 in range(20, 60):
            for A1 in range(150, 180):
                for A2 in range(280, 360):
                    # minVerts = math.floor((A2 - A0 - math.pi) / (math.pi - A1)) + 1
                    # maxVerts = math.ceil(math.pi / (math.pi - A1)) - 1
                    # if minVerts < 0:
                    #     minVerts = 1
                    # if maxVerts < 0:
                    #     maxVerts = 1
                    minVerts = 3
                    maxVerts = 20
                    recursivelyBuildTerrains(0, [], A0, A1, A2, minVerts, maxVerts, numTerrains)
    else:
        T = generateTerrainPoints(20, 150, 271, 3, [3, 5, 5])
        writeDataToFile(run_seq, T)
        sol = bruteForce(T)
        plotTerrain(T, sol, run_seq)
        run_seq = str(int(run_seq) + 1)
    
    os.environ['RUN'] = str(int(run_seq))
    dotenv.set_key(dotenv_file, "RUN", os.environ["RUN"])

