# distanceCalculator.py
# ---------------------
# Licensing Information: Please do not distribute or publish solutions to this
# project. You are free to use and extend these projects for educational
# purposes. The Pacman AI projects were developed at UC Berkeley, primarily by
# John DeNero (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# For more info, see http://inst.eecs.berkeley.edu/~cs188/sp09/pacman.html

"""
This file contains a Distancer object which computes and 
caches the shortest path between any two points in the maze. 

Example:
distancer = Distancer(gameState.data.layout)
distancer.getDistance( (1,1), (10,10) )
"""

import sys, time, random

class Distancer:
  def __init__(self, layout, default = 10000):
    """
    Initialize with Distancer(layout).  Changing default is unnecessary.    
    """
    self._distances = None
    self.default = default    
    self.dc = DistanceCalculator(layout, self, default)

  def getMazeDistances(self):
    self.dc.run()
    
  def getDistance(self, pos1, pos2):
    """
    The getDistance function is the only one you'll need after you create the object.
    """
    if self._distances == None:
      return manhattanDistance(pos1, pos2)
    if isInt(pos1) and isInt(pos2):
      return self.getDistanceOnGrid(pos1, pos2)
    pos1Grids = getGrids2D(pos1)
    pos2Grids = getGrids2D(pos2)
    bestDistance = self.default
    for pos1Snap, snap1Distance in pos1Grids:
      for pos2Snap, snap2Distance in pos2Grids:
        gridDistance = self.getDistanceOnGrid(pos1Snap, pos2Snap)
        distance = gridDistance + snap1Distance + snap2Distance
        if bestDistance > distance:
          bestDistance = distance
    return bestDistance

  def getDistanceOnGrid(self, pos1, pos2):
    key = (pos1, pos2)
    if key in self._distances:
      return self._distances[key]
    else:
      raise Exception("Positions not in grid: " + str(key))

  def isReadyForMazeDistance(self):
    return self._distances != None

def manhattanDistance(x, y ):
  return abs( x[0] - y[0] ) + abs( x[1] - y[1] )

def isInt(pos):
  x, y = pos
  return x == int(x) and y == int(y)

def getGrids2D(pos):
  grids = []
  for x, xDistance in getGrids1D(pos[0]):
    for y, yDistance in getGrids1D(pos[1]):
      grids.append(((x, y), xDistance + yDistance))
  return grids
  
def getGrids1D(x):
  intX = int(x)
  if x == int(x):
    return [(x, 0)]
  return [(intX, x-intX), (intX+1, intX+1-x)]
  
##########################################
# MACHINERY FOR COMPUTING MAZE DISTANCES #
##########################################

distanceMap = {}

class DistanceCalculator:
  def __init__(self, layout, distancer, default = 10000):
    self.layout = layout
    self.distancer = distancer
    self.default = default
  
  def run(self):
    global distanceMap

    if self.layout.walls not in distanceMap:
      distances = computeDistances(self.layout)
      distanceMap[self.layout.walls] = distances
    else:
      distances = distanceMap[self.layout.walls]

    self.distancer._distances = distances  

def computeDistances(layout):
    "Runs UCS to all other positions from each position"
    distances = {}
    allNodes = layout.walls.asList(False)
    for source in allNodes:
        dist = {}
        closed = {}
        for node in allNodes:
            dist[node] = sys.maxint
        import util
        queue = util.PriorityQueue()
        queue.push(source, 0)
        dist[source] = 0
        while not queue.isEmpty():
            node = queue.pop()
            if node in closed:
                continue
            closed[node] = True
            nodeDist = dist[node]
            adjacent = []
            x, y = node
            if not layout.isWall((x,y+1)):
                adjacent.append((x,y+1))
            if not layout.isWall((x,y-1)):
                adjacent.append((x,y-1) )
            if not layout.isWall((x+1,y)):
                adjacent.append((x+1,y) )
            if not layout.isWall((x-1,y)):
                adjacent.append((x-1,y))
            for other in adjacent:
                if not other in dist:
                    continue
                oldDist = dist[other]
                newDist = nodeDist+1
                if newDist < oldDist:
                    dist[other] = newDist
                    queue.push(other, newDist)
        for target in allNodes:
            distances[(target, source)] = dist[target]
    return distances
    

def getDistanceOnGrid(distances, pos1, pos2):
    key = (pos1, pos2)
    if key in distances:
      return distances[key]
    return 100000
  
