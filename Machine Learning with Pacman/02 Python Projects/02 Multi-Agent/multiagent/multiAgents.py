# multiAgents.py
# --------------
# Licensing Information: Please do not distribute or publish solutions to this
# project. You are free to use and extend these projects for educational
# purposes. The Pacman AI projects were developed at UC Berkeley, primarily by
# John DeNero (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# Student side autograding was added by Brad Miller, Nick Hay, and Pieter
# Abbeel in Spring 2013.
# For more info, see http://inst.eecs.berkeley.edu/~cs188/pacman/pacman.html
 

from util import manhattanDistance
from game import Directions
import random, util
from game import Agent
 

class ReflexAgent(Agent):
    """
      A reflex agent chooses an action at each choice point by examining
      its alternatives via a state evaluation function.
      The code below is provided as a guide.  You are welcome to change
      it in any way you see fit, so long as you don't touch our method
      headers.
    """

    def getAction(self, gameState):
        """
        You do not need to change this method, but you're welcome to.
        getAction chooses among the best options according to the evaluation function.
        Just like in the previous project, getAction takes a GameState and returns
        some Directions.X for some X in the set {North, South, West, East, Stop}
        """
        # Collect legal moves and successor states
        legalMoves = gameState.getLegalActions()

        # Choose one of the best actions
        scores = [self.evaluationFunction(gameState, action) for action in legalMoves ]
        bestScore = max(scores)
        bestIndices = [index for index in range(len(scores)) if scores[index] == bestScore]
        chosenIndex = random.choice(bestIndices) # Pick randomly among the best

        "Add more of your code here if you want to"
        return legalMoves[chosenIndex]

    def evaluationFunction(self, currentGameState, action):
        """
        Design a better evaluation function here.
        The evaluation function takes in the current and proposed successor
        GameStates (pacman.py) and returns a number, where higher numbers are better.
        The code below extracts some useful information from the state, like the
        remaining food (newFood) and Pacman position after moving (newPos).
        newScaredTimes holds the number of moves that each ghost will remain
        scared because of Pacman having eaten a power pellet.
        Print out these variables to see what you're getting, then combine them
        to create a masterful evaluation function.
        """

        # Useful information you can extract from a GameState (pacman.py)
        successorGameState = currentGameState.generatePacmanSuccessor(action)
        newPos = successorGameState.getPacmanPosition()
        newFood = successorGameState.getFood()
        newGhostStates = successorGameState.getGhostStates()
        newScaredTimes = [ghostState.scaredTimer for ghostState in newGhostStates]
 
        "*** YOUR CODE HERE ***"
 
        # Initialize Variables
        pacmanPosition = successorGameState.getPacmanPosition()
        ghostStates = successorGameState.getGhostStates()
        ghostPositions = successorGameState.getGhostPositions()
        scaredTimes = [ghostState.scaredTimer for ghostState in newGhostStates]
        foodPositions = successorGameState.getFood().asList()
        foodRemaining = successorGameState.getNumFood()
 
        # Closest Ghost and Ghost Index
        closestGhost = 999999
        closestGhostIndex = 0
        ghostCount = 0
 
        for ghostPosition in ghostPositions:
            ghostCount += 1
            ghostDistance = util.manhattanDistance( pacmanPosition, ghostPosition )

            if ghostDistance < closestGhost:
                closestGhostIndex = ghostCount
                closestGhost = ghostDistance

        # Are there any Ghosts?
        if closestGhostIndex == 0:
           ghostsInMaze = False
        else:
           ghostsInMaze = True

        # Is the Closest Ghost Scared? - with 5 steps on the clock
        isGhostScared = False
        if ghostsInMaze == True:
            if successorGameState.getGhostState( closestGhostIndex )\
            .scaredTimer > 5:

                isGhostScared = True
 
        # Closest Ghost Scared Time
        ghostScaredTime = successorGameState.getGhostState( closestGhostIndex )

        # Ignore Ghosts if they are more than 5 steps away
        if closestGhost >=5 and isGhostScared == False:
            closestGhost = 999999.9
 
        # Food Distances
        foodDistances = list()
        for foodPosition in foodPositions:
          foodDistance = util.manhattanDistance( pacmanPosition, foodPosition )
          foodDistances.append( foodDistance )
 
        # Closest Food
        closestFood = 999999.9
        for foodDistance in foodDistances:
          closestFood = min( foodDistance, closestFood )
 
        # Inverse Closest Food
        if closestFood <= 0:
            inverseClosestFood = 999999
        else:         
            inverseClosestFood = 1.0 / closestFood
 
        # Chase Ghosts if they are scared
        if ghostsInMaze == True and isGhostScared == True:
            if closestGhost == 0:
                # Eat Ghost for big points!
                closestGhost = 999999
            else:             
                # Can't allow closest Ghost to decrease when ghosts are scared
                closestGhost += ( 999999 / closestGhost )
         
        # Score
        score = successorGameState.getScore() + closestGhost \
        + inverseClosestFood + foodRemaining
     
        return score

def scoreEvaluationFunction(currentGameState):
    """
      This default evaluation function just returns the score of the state.
      The score is the same one displayed in the Pacman GUI.
      This evaluation function is meant for use with adversarial search agents
      (not reflex agents).
    """
    return currentGameState.getScore()
 

class MultiAgentSearchAgent(Agent):
    """
      This class provides some common elements to all of your
      multi-agent searchers.  Any methods defined here will be available
      to the MinimaxPacmanAgent, AlphaBetaPacmanAgent & ExpectimaxPacmanAgent.
      You *do not* need to make any changes here, but you can if you want to
      add functionality to all your adversarial search agents.  Please do not
      remove anything, however.
      Note: this is an abstract class: one that should not be instantiated.  It's
      only partially specified, and designed to be extended.  Agent (game.py)
      is another abstract class.
    """

    def __init__(self, evalFn = 'scoreEvaluationFunction', depth = '2'):
        self.index = 0 # Pacman is always agent index 0
        self.evaluationFunction = util.lookup(evalFn, globals())
        self.depth = int(depth)
 
class MinimaxAgent(MultiAgentSearchAgent):
    """
      Your minimax agent (question 2)
    """
 
    def getAction(self, gameState):
        """
          Returns the minimax action from the current gameState using self.depth
          and self.evaluationFunction.
          Here are some method calls that might be useful when implementing minimax.
          gameState.getLegalActions(agentIndex):
          Returns a list of legal actions for an agent
            agentIndex=0 means Pacman, ghosts are >= 1
          gameState.generateSuccessor(agentIndex, action):
            Returns the successor game state after an agent takes an action
          gameState.getNumAgents():
            Returns the total number of agents in the game
        """
 
        "*** YOUR CODE HERE ***"

        # Constant Definitions
        PACMAN              = 0
        MAX_DEPTH           = self.depth
        MAX                 = 'MAX'
        MIN                 = 'MIN'
        TERMINUS            = 'TERMINUS'
        RESULT              = 'RESULT'

        # Initialize Variables
        rootPosition        = 0
        childPosition       = 0
        currentDepth        = 0
        currentAgent        = 0
        currentState        = gameState
        nAgents             = gameState.getNumAgents()
        resultAction        = Directions.STOP

        # Debug Information
        #############################
        DEBUG_SETUP         = False
        DEBUG_NODES         = False
        DEBUG_SEARCH        = False
        DEBUG_DICT          = False
        DEBUG_RESULT        = False
        #############################

        # Initialize Placeholder Objects
        nodeStack           = util.Stack()  # ( rootAction, operator, nodeState, parent, depth, agent )
        searchTree          = dict()        # ( parent, operator, legalActions, score, rootAction, nodesToExpand )

        # Initialize the Root
        parent              = RESULT
        currentOperator     = MAX
        nodesToExpand       = []
        score               = -999999.9
        rootAction          = Directions.STOP

        legalActions        = []
        if gameState.isWin() == False and gameState.isLose() == False:
            legalActions    = gameState.getLegalActions( PACMAN )

        # Problem Set-Up
        ################

        if DEBUG_SETUP:
            print ''
            print 'PROBLEM SET-UP'
            print '--------------'
            print 'nAgents:',nAgents
            print 'nGhosts',nAgents - 1
            print 'currentAgent:',currentAgent
            print 'currentDepth:',currentDepth
            print 'MAX_DEPTH:',MAX_DEPTH
            print ''
            
        # Push Search Tree Dictionary: Add Root to Search Tree
        # Thereafter ... Only Add to Tree when Popping Nodes
        searchTree[ rootPosition ] = ( parent, currentOperator, legalActions, score, rootAction, nodesToExpand )                

        if DEBUG_SETUP:
            print 'PUSH SEARCH TREE (ROOT)'
            print 'node:', rootPosition
            print 'parent:', parent
            print 'operator:', currentOperator
            print 'legalActions:', legalActions
            print 'score:', score
            print 'rootAction:', rootAction
            print 'nodesToExpand:',nodesToExpand
            print ''

        # Generate Successors and Push to Node Stack
        # ------------------------------------------

        # Successor: currentAgent, currentDepth
        if nAgents == 1:
            currentDepth +=1
            
        elif nAgents > 1:
            currentAgent += 1

            if currentAgent >= nAgents:
                currentAgent = PACMAN
                currentDepth += 1

        if DEBUG_SETUP:
            print 'SUCCESSORS'
            print 'nextAgent:',currentAgent
            print 'nextDepth:',currentDepth
            print ''

        # Iterate through Legal Actions                    
        if legalActions == []:
            return []
        
        for legalAction in legalActions:

            # Generate Successor
            successorState  = currentState.generateSuccessor( PACMAN, legalAction )

            # Evaluate Successor
            currentOperator = MIN
            if successorState.isWin() == True or successorState.isLose() == True or ( currentDepth >= MAX_DEPTH and currentAgent == PACMAN ):
                currentOperator = TERMINUS

            elif currentAgent == PACMAN:
                currentOperator = MAX

            # Push Node
            node = ( legalAction, currentOperator, successorState, rootPosition, currentDepth, currentAgent )
            nodeStack.push ( node )

            if DEBUG_SETUP:
                print 'PUSH NODE (ROOT CHILD)'
                print 'legalAction:',legalAction
                print 'operator:',currentOperator
                print 'parent:',rootPosition
                print 'depth:',currentDepth
                print 'agent:',currentAgent
                print ''














        # Node Stack Loop
        #################
        if DEBUG_NODES:
            print 'NODE LOOP'
            print '---------'
            print ''

        nodeLoop = 0
        while nodeStack.isEmpty() == False:
            nodeLoop += 1

            # Pop Node Stack
            rootAction, currentOperator, currentState, nodeParent, nodeDepth, nodeAgent = nodeStack.pop()
            
            # Update Positions for Search Tree
            parentPosition = nodeParent
            childPosition += 1
            
            if DEBUG_NODES:
                print 'NODE LOOP ',nodeLoop
                print '-------------'
                print 'currentAgent:',currentAgent
                print 'currentDepth:',currentDepth
                print ''
                print 'POP NODE'
                print 'legalAction:',legalAction
                print 'operator:',currentOperator
                print 'nodeParent:',nodeParent
                print 'nodeDepth:',nodeDepth
                print 'nodeAgent:',nodeAgent
                print ''

            # Update currentAgent and currentDepth
            currentDepth = nodeDepth
            currentAgent = nodeAgent

            if nAgents == 1:
                currentDepth += 1

            elif nAgents > 1:
                currentAgent += 1

                if currentAgent >= nAgents:
                    currentAgent = PACMAN
                    currentDepth += 1

            # Check for Terminal States
            if currentOperator == TERMINUS or ( currentDepth >= MAX_DEPTH and currentAgent != PACMAN ) or currentState.isWin() or currentState.isLose() or legalActions == []:

                
                # Update Search Tree on Every Node Popped
                score = self.evaluationFunction( currentState )
                
                nodesToExpand = []
                searchTree[ childPosition ] = ( parentPosition, currentOperator, [] , score, rootAction, nodesToExpand )

                # Update SearchTree Parent Nodes To Expand
                p, o, l, s, a, e = searchTree[ parentPosition ]
                e.append( childPosition )
                searchTree[ parentPosition ] = ( p, o, l, s, a, e )

                if DEBUG_NODES:
                    print 'PUSH SEARCH TREE (TERMINUS)'
                    print '---------------------------'
                    print 'node:',childPosition
                    print 'parentPosition:',parentPosition
                    print 'operator',currentOperator
                    print 'legalActions:',Directions.STOP
                    print 'score:',score
                    print 'rootAction:',rootAction
                    print 'nodesToExpand:',nodesToExpand
                    print ''
                
            else:

                # Update Search Tree on Every Node Popped
                score = 999999.9
                if currentOperator == MAX:
                    score *= -1.0

                nodesToExpand = []                    
                searchTree[ childPosition ] = ( parentPosition, currentOperator, legalActions, score, rootAction, nodesToExpand )

                # Update SearchTree Parent Nodes To Expand
                p, o, l, s, a, e = searchTree[ parentPosition ]
                e.append( childPosition )
                searchTree[ parentPosition ] = ( p, o, l, s, a, e )
                
                if DEBUG_NODES:
                    print 'PUSH SEARCH TREE (NON_TERMINUS)'
                    print '-------------------------------'
                    print 'node:',childPosition
                    print 'parentPosition:',parentPosition
                    print 'operator',currentOperator
                    print 'legalActions:',legalActions
                    print 'score:',score
                    print 'rootAction:',rootAction
                    print 'nodesToExpand:',nodesToExpand
                    print ''


                # Successor: nextAgent, nextAgent
                nextAgent = currentAgent
                nextDepth = currentDepth

                # Update the Successor Operator and Score
                nextOperator = MIN
                score = 999999.9

                if nextAgent == PACMAN:
                    nextOperator = MAX 
                    score *= -1.0

                ##################################################################
                ##################################################################
                # CHECK
                # Don't generate next successors if current Depth Breached
                
                ##################################################################
                ##################################################################
                    
                # Iterate through Legal Actions
                nextLegalActions = []
                if currentState.isWin() == False and currentState.isLose() == False:
                    nextLegalActions = currentState.getLegalActions( PACMAN )
                    
                for nextLegalAction in nextLegalActions:

                    # Generate Next Successor if no breach of current depth
                    if ( currentDepth >= MAX_DEPTH and currentAgent != PACMAN ):
                        nextOperator = TERMINUS

                    else:
                        
                        successorState  = currentState.generateSuccessor( PACMAN, nextLegalAction )

                        nextOperator = MIN

                        if successorState.isWin() or successorState.isLose() or ( currentDepth >= MAX_DEPTH and currentAgent != PACMAN ) or currentOperator == TERMINUS or legalActions == []:
                            nextOperator = TERMINUS

                        elif nextAgent == PACMAN:
                            nextOperator = MAX
                        
                        # Push Node
                        parentPosition = childPosition
                        node = ( nextLegalAction, nextOperator, successorState, parentPosition, currentDepth, currentAgent )
                        nodeStack.push ( node )

                        if DEBUG_NODES:
                            print 'PUSH NODE'
                            print 'legalAction:',nextLegalAction
                            print 'operator:',nextOperator
                            print ''




















        # Search Tree Loop
        # ----------------
        if DEBUG_SEARCH:
            print 'SEARCH TREE LOOP'
            print '----------------'
            print ''

        for key in searchTree:
            parent, operator, legalActions, score, rootAction, nodesToExpand = searchTree[key]

            if DEBUG_SEARCH:
                print 'key:',key
                print '---------'
                print 'parent:',parent
                print 'operator:',operator
                print 'legalActions',legalActions
                print 'score:',score
                print 'rootAction:',rootAction
                print 'nodesToExpand:',nodesToExpand
                print ''


        if DEBUG_SEARCH:
            print 'OPERATE'
            print '#######'
            print ''
            
        # Navigate Search Tree
        # --------------------
        
        # KEY
        # (n)ode = ( (p)arent, (o)perator, (l)egalActions, (s)core, root(A)ction, nodesTo(E)xpand )
        
        n = 0 
        p, o, l, s, a, e = searchTree[ n ]

        if DEBUG_SEARCH:
            print 'LOOP 0'
            print 'n:',n
            print 'p:',p
            print 'o:',o
            print 'l:',l
            print 's:',s
            print 'a:',a
            print 'e:',e
            print ''

        searchCount = 0            
        while not( n == 0 and len( e ) == 0 ):
            searchCount += 1

            if DEBUG_SEARCH:
                print 'LOOP ', searchCount
            
            # update start node, popping e
            copyn = n
            copyp = p
            copyo = o
            copyl = l
            copys = s
            copya = a
            copye = e
            
            # If Terminal Node or no more nodes to expand
            if copyo == TERMINUS or len( e ) == 0:

                # move up a node
                if DEBUG_SEARCH:
                    print 'MOVE UP'
                n = copyp
                p, o, l, s, a, e = searchTree[ n ]
                
                # Apply Operator, Update Score and Update Root Action on Root Nodes
                if DEBUG_SEARCH:
                    print 'o:',o
                    print 'copys:',copys
                    print 's:',s
                if o == MAX:
                    if DEBUG_SEARCH:
                        print 'APPLY MAX'
                    if copys > s:
                        s = copys
                        if n == 0:
                            a = copya
                        if DEBUG_SEARCH:
                            print 'new s:', s
                            print 'a:', a

                elif o == MIN:
                    if DEBUG_SEARCH:
                        print 'APPLY MIN'
                    if copys < s:
                        s = copys
                        if n == 0:
                            a = copya
                        if DEBUG_SEARCH:
                            print 'new s:', s
                            print 'a:', a

                # update node
                searchTree[ n ] = ( p, o, l, s, a, e )
                if DEBUG_SEARCH:
                    print 'UP TREE:',searchTree[ n ]
                
            elif len( e ) > 0:

                # Pop the expanded node
                lowern = e.pop()
                copye = e
            
                # Remove expanded node
                copye = e # after pop
                searchTree[ copyn ] = ( copyp, copyo, copyl, copys, copya, copye )

                # move down a node
                if DEBUG_SEARCH:
                    print 'MOVE DOWN'
                n = lowern
                p, o, l, s, a, e = searchTree[ n ]


            if DEBUG_SEARCH:            
                print 'n:',n
                print 'p:',p
                print 'o:',o
                print 'l:',l
                print 's:',s
                print 'a:',a
                print 'e:',e
                print ''
        
        
        if DEBUG_DICT or DEBUG_RESULT:
            print 'ROOT'
            print '----'
            print 'rootParent:',p
            print 'rootOperator:',o
            print 'rootLegalActions:',l
            print 'nRootLegalActions:',len(l)
            print 'rootScore:',s
            print 'rootAction:',a
            print 'nodesToExpand:',e
            print ''

        if DEBUG_DICT:
            print 'SEARCH TREE'
            print '-----------'
            for key in searchTree:
                print 'key',key
                print searchTree[ key ]
                print ''
                    
        resultAction = a
        resultScore = s

        if DEBUG_DICT or DEBUG_RESULT:
            print ''
            print 'resultAction:',resultAction
            print 'resultScore:',resultScore
            print ''

        return resultAction


 
 
class AlphaBetaAgent(MultiAgentSearchAgent):
    """
      Your minimax agent with alpha-beta pruning (question 3)
    """
 
    def getAction(self, gameState):
        """
          Returns the minimax action using self.depth and self.evaluationFunction
        """
        "*** YOUR CODE HERE ***"
        util.raiseNotDefined()
 
class ExpectimaxAgent(MultiAgentSearchAgent):
    """
      Your expectimax agent (question 4)
    """
 
    def getAction(self, gameState):
        """
          Returns the expectimax action using self.depth and self.evaluationFunction
          All ghosts should be modeled as choosing uniformly at random from their
          legal moves.
        """
        "*** YOUR CODE HERE ***"
        util.raiseNotDefined()
 
def betterEvaluationFunction(currentGameState):
    """
      Your extreme ghost-hunting, pellet-nabbing, food-gobbling, unstoppable
      evaluation function (question 5).
      DESCRIPTION: <write something here so we know what you did>
    """
    "*** YOUR CODE HERE ***"
    util.raiseNotDefined()
 
# Abbreviation
better = betterEvaluationFunction
class ContestAgent(MultiAgentSearchAgent):
    """
      Your agent for the mini-contest
    """
 
    def getAction(self, gameState):
        """
          Returns an action.  You can use any method you want and search to any depth you want.
          Just remember that the mini-contest is timed, so you have to trade off speed and computation.
          Ghosts don't behave randomly anymore, but they aren't perfect either -- they'll usually
          just make a beeline straight towards Pacman (or away from him if they're scared!)
        """
        "*** YOUR CODE HERE ***"
        util.raiseNotDefined()
