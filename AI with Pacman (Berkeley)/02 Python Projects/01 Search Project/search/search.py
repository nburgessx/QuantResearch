# search.py
# ---------
# Licensing Information: Please do not distribute or publish solutions to this
# project. You are free to use and extend these projects for educational
# purposes. The Pacman AI projects were developed at UC Berkeley, primarily by
# John DeNero (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# Student side autograding was added by Brad Miller, Nick Hay, and Pieter
# Abbeel in Spring 2013.
# For more info, see http://inst.eecs.berkeley.edu/~cs188/pacman/pacman.html
"""
In search.py, you will implement generic search algorithms which are called
by Pacman agents (in searchAgents.py).
"""
import util
class SearchProblem:
    """
    This class outlines the structure of a search problem, but doesn't implement
    any of the methods (in object-oriented terminology: an abstract class).
    You do not need to change anything in this class, ever.
    """
    def getStartState(self):
        """
        Returns the start state for the search problem
        """
        util.raiseNotDefined()
    def isGoalState(self, state):
        """
          state: Search state
        Returns True if and only if the state is a valid goal state
        """
        util.raiseNotDefined()
    def getSuccessors(self, state):
        """
          state: Search state
        For a given state, this should return a list of triples,
        (successor, action, stepCost), where 'successor' is a
        successor to the current state, 'action' is the action
        required to get there, and 'stepCost' is the incremental
        cost of expanding to that successor
        """
        util.raiseNotDefined()
    def getCostOfActions(self, actions):
        """
         actions: A list of actions to take
        This method returns the total cost of a particular sequence of actions.  The sequence must
       be composed of legal moves
        """
        util.raiseNotDefined()
 
def tinyMazeSearch(problem):
    """
    Returns a sequence of moves that solves tinyMaze.  For any other
    maze, the sequence of moves will be incorrect, so only use this for tinyMaze
    """
    from game import Directions
    s = Directions.SOUTH
    w = Directions.WEST
    return  [s,s,w,s,w,w,s,w]
def depthFirstSearch(problem):
    """
    Search the deepest nodes in the search tree first
    Your search algorithm needs to return a list of actions that reaches
    the goal.  Make sure to implement a graph search algorithm
    To get started, you might want to try some of these simple commands to
    understand the search problem that is being passed in:
    print "Start:", problem.getStartState()
    print "Is the start a goal?", problem.isGoalState(problem.getStartState())
    print "Start's successors:", problem.getSuccessors(problem.getStartState())
    """
    "*** YOUR CODE HERE ***"
    # Initialise Variables
    closed = set()
    fringe = util.Stack() #LIFO
    node = util.Queue() #FIFO
    foundGoal = False
    #Initialise Start Fringe, which consists of Position, Action and StepCost
    startNode = ( problem.getStartState(), list(), 0 )
    fringe.push( startNode )
    counter = 0
    while fringe.isEmpty() == False:
        counter = counter + 1
      
        # While Loop
        node = fringe.pop()
        parent_position, node_path, parent_stepCost = node
        # Check for Goal State and if found return the Path to the Goal
        if problem.isGoalState( parent_position ):
            foundGoal = True
            break
        # Expansion Routine
        if parent_position not in closed:
            # Add Node to Closed Set
            closed.add( parent_position )
            # Expand Node
            for child in problem.getSuccessors( parent_position ):
                # Reset Child Path
                child_path = list()
                for element in node_path:
                    child_path.append( element )
              
                # Get Child Details
                child_position, child_action, child_stepCost = child
                # Update Path to Child
                child_path.append( child_action )
              
                # Create an Updated Child Node
                childNode = ( child_position, child_path, parent_stepCost + child_stepCost )
                # Add Updated Child to Fringe
                fringe.push( childNode )
    # Return the Result
    if foundGoal == False:
  
        # No Solution Found
        print 'Failure: No Solution Found'
        return fringe # empty fringe
  
    else:
        # Print Successful Path Queue
        print 'Solution:', node_path
        return node_path
        #tinyMaze Solution
        #return ['South','South','West','South','West','West','South','West']
  
def breadthFirstSearch(problem):
    """
    Search the shallowest nodes in the search tree first.
    """
    "*** YOUR CODE HERE ***"
 
    # Initialise Variables
    closed = set()
    fringe = util.Queue() #FIFO
    node = util.Queue() #FIFO
    foundGoal = False
    #Initialise Start Fringe, which consists of Position, Action and StepCost
    startNode = ( problem.getStartState(), list(), 0 )
    fringe.push( startNode )
    counter = 0
    while fringe.isEmpty() == False:
        counter = counter + 1
      
        # While Loop
        node = fringe.pop()
        parent_position, node_path, parent_stepCost = node
        # Check for Goal State and if found return the Path to the Goal
        if problem.isGoalState( parent_position ):
            foundGoal = True
            break
        # Expansion Routine
        if parent_position not in closed:
            # Add Node to Closed Set
            closed.add( parent_position )
            # Expand Node
            for child in problem.getSuccessors( parent_position ):
                # Reset Child Path
                child_path = list()
                for element in node_path:
                    child_path.append( element )
              
                # Get Child Details
                child_position, child_action, child_stepCost = child
                # Update Path to Child
                child_path.append( child_action )
              
                # Create an Updated Child Node
                childNode = ( child_position, child_path, parent_stepCost + child_stepCost )
                # Add Updated Child to Fringe
                fringe.push( childNode )
    # Return the Result
    if foundGoal == False:
  
        # No Solution Found
        print 'Failure: No Solution Found'
        return fringe # empty fringe
  
    else:
        # Print Successful Path Queue
        print 'Solution:', node_path
        return node_path
        #tinyMaze Solution
        #return ['South','South','West','South','West','West','South','West']
 
def uniformCostSearch(problem):
    "Search the node of least total cost first. "
    "*** YOUR CODE HERE ***"
 
    # Initialise Variables
    closed = set()
    fringe = util.PriorityQueue() #Priority Queue
    node = util.Queue() #FIFO
    foundGoal = False

    #Initialise Start Fringe, which consists of Position, Action and StepCost
    startNode = ( problem.getStartState(), list(), 0 )
    fringe.push( startNode, 0 )

    while fringe.isEmpty() == False:
      
        # While Loop
        node = fringe.pop()
        parent_position, node_path, parent_stepCost = node

        # Check for Goal State and if found return the Path to the Goal
        if problem.isGoalState( parent_position ):
            foundGoal = True
            break

        # Expansion Routine
        if parent_position not in closed:

            # Add Node to Closed Set
            closed.add( parent_position )

            # Expand Node
            for child in problem.getSuccessors( parent_position ):

                # Reset Child Path
                child_path = list()
                for element in node_path:
                    child_path.append( element )
              
                # Get Child Details
                child_position, child_action, child_stepCost = child

                # Update Path to Child
                child_path.append( child_action )
              
                # Create an Updated Child Node
                childNode = ( child_position, child_path, parent_stepCost + child_stepCost )

                # Add Updated Child to Fringe
                fringe.push( childNode, parent_stepCost + child_stepCost )

    # Return the Result
    if foundGoal == False:
  
        # No Solution Found
        print 'Failure: No Solution Found'
        return fringe # empty fringe
  
    else:

        # Print Successful Path Queue
        print 'Solution:', node_path
        return node_path
        #tinyMaze Solution
        #return ['South','South','West','South','West','West','South','West']

def nullHeuristic(state, problem=None):
    """
    A heuristic function estimates the cost from the current state to the nearest
    goal in the provided SearchProblem.  This heuristic is trivial.
    """
    return 0

def aStarSearch(problem, heuristic=nullHeuristic):
    "Search the node that has the lowest combined cost and heuristic first."
    "*** YOUR CODE HERE ***"
 
    # Initialise Variables
    foundGoal                   = False
    closed                      = {}
    path                        = {}
    pathgCost                   = {}
    pathfCost                   = {}
    moreInfo                    = False

    # State: ( Position, Cost )
    position                    = problem.getStartState()
    gCost                       = 0
    hCost                       = heuristic( position, problem )
    fCost                       = gCost + hCost
    startNode                   = ( position, gCost )
    if moreInfo == True: print 'startNode:', startNode
    
    path[ position ]            = list()
    pathgCost[ position ]       = gCost
    pathfCost[ position ]       = fCost

    ################################################################################################################    
    # Fringe Set-Up                                                                                                #
    # IMPORTANT: Finge Queue is given a function to prioritise the 'fringe.pop()' selection function               #
    # For an aStarSearch this function should represent the totalCost = path cost + heuristic cost i.e. f = g + h  #
    ################################################################################################################
    fringe = util.PriorityQueueWithFunction( lambda node: pathgCost[ node[0] ] + heuristic( node[0], problem ) )
    fringe.push( startNode )
    
    # Fringe Loop
    counter = 0
    while fringe.isEmpty() == False:

        counter += 1
        if moreInfo == True: print 'count:', counter
        
        # Update Parent Node
        parentNode              = fringe.pop()
        
        parentPosition          = parentNode[0]
        parentgCost             = parentNode[1] 
        parentfCost             = parentgCost + heuristic( parentPosition, problem )
        if moreInfo == True: print '#####################################'
        if moreInfo == True: print 'expansion:',counter
        if moreInfo == True: print 'selected:', parentNode

        # Don't Allow Costs to be Negative
        if parentgCost < 0:
            parentgCost = 0

        parentPath                  = path[ parentPosition ]
        pathgCost[ parentPosition ] = parentgCost
        if moreInfo == True: print 'parentPath:', parentPath
        if moreInfo == True: print 'parentgCost:', parentgCost

        # Goal Test
        if problem.isGoalState( parentPosition ) == True:
            foundGoal = True
            break
        
        # Check if Closed Dictionary contains Node and switch closed node if new node cost is smaller
        isInClosed = closed.has_key( parentPosition )
        if moreInfo == True: print 'IsInClosedSet:',isInClosed
        
        # MEMORY EFFICIENT GRAPH SEARCH FEATURE
        #################################################################
        isNewCostSmaller = False
        if isInClosed:
            isNewCostSmaller = parentgCost < closed[ parentPosition ]
        #################################################################
                
        # Expand Node
        if isInClosed == False or isNewCostSmaller == True:
        
                closed[ parentPosition ] = parentgCost
                if moreInfo == True: print 'closed:', closed
        
                for child in problem.getSuccessors( parentPosition ):

                    # Extract Child Information
                    childPosition, childAction, childgCost = child
                    if moreInfo == True: print 'childPosition:', childPosition
                    if moreInfo == True: print 'childAction:', childAction

                    # Update Child Costs                
                    childgCost                  = pathgCost[ parentPosition ] + childgCost
                    if childgCost < 0:
                        childgCost = 0
                        
                    childfCost                  = childgCost + heuristic( childPosition, problem )
                    if childfCost < 0:
                        childfCost = 0

                    ############################################################################
                    # Update Path, Costs and Dictionaries                                      #
                    # IMPORTANT                                                                #
                    # Only Update paths if they are f-cost cheaper than any existing path data #
                    ############################################################################

                    isfCostCheaper = True
                    doesPathAlreadyExist = path.has_key( childPosition )
                    if doesPathAlreadyExist:
                        isfCostCheaper = childfCost < pathfCost[ childPosition ]

                    if moreInfo == True: print 'doesPathAlreadyExist:', doesPathAlreadyExist
                    if moreInfo == True: print 'isfCostCheaper:', isfCostCheaper

                    if isfCostCheaper == True:
                    
                        # Update Child Path
                        childPath = list()
                        for action in path[ parentPosition ]:
                            childPath.append( action )
                            if moreInfo == True: print'action:',action
                        if moreInfo == True: print 'parentPath:', childPath

                        childPath.append( childAction )                    
                        if moreInfo == True: print 'childPath:', childPath

                        # Update Path and Cost Dictionaries
                        path[ childPosition ]       = childPath                                    
                        pathgCost[ childPosition ]  = childgCost                                   
                        pathfCost[ childPosition ]  = childfCost
                        if moreInfo == True: print 'path:', path[ childPosition ]
                        if moreInfo == True: print 'pathgCost:', pathgCost[ childPosition ]
                        if moreInfo == True: print 'pathfCost:', pathfCost[ childPosition ]

                    #############################################################################
                    #############################################################################
                        
                    # Create Child Node and Add to Fringe
                    childNode = ( childPosition, childgCost )
                    #print 'childNode:', childNode

                    fringe.push( childNode )
                    if moreInfo == True: print 'childNode:', childNode
                    if moreInfo == True: print '---------------------------------'
            
    # Return Failure or Solution
    if foundGoal == False:
        if moreInfo == True: print 'Failure: No Solution Found'
        if moreInfo == True: print 'last attempted path:', parentPath
        return list()
    else:
        if moreInfo == True: print 'Solution:', parentPath
        return parentPath

 
# Abbreviations
bfs = breadthFirstSearch
dfs = depthFirstSearch
astar = aStarSearch
ucs = uniformCostSearch

