# tutorialTestClasses.py
# ----------------------
# Licensing Information: Please do not distribute or publish solutions to this
# project. You are free to use and extend these projects for educational
# purposes. The Pacman AI projects were developed at UC Berkeley, primarily by
# John DeNero (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# Student side autograding was added by Brad Miller, Nick Hay, and Pieter 
# Abbeel in Spring 2013.
# For more info, see http://inst.eecs.berkeley.edu/~cs188/pacman/pacman.html

import testClasses

# Simple test case which evals an arbitrary piece of python code.
# The test is correct if the output of the code given the student's
# solution matches that of the instructor's.
class EvalTest(testClasses.TestCase):

    def __init__(self, testDict):
        super(EvalTest, self).__init__(testDict)
        self.preamble = compile(testDict.get('preamble', ""), "%s.preamble" % self.getPath(), 'exec')
        self.test = compile(testDict['test'], "%s.test" % self.getPath(), 'eval')
        self.success = testDict['success']
        self.failure = testDict['failure']

    def evalCode(self, moduleDict):
        bindings = dict(moduleDict)
        exec self.preamble in bindings
        return str(eval(self.test, bindings))

    def execute(self, grades, moduleDict, solutionDict):
        result = self.evalCode(moduleDict)
        if result == solutionDict['result']:
            grades.addMessage('PASS: %s' % self.path) 
            grades.addMessage('\t%s' % self.success)
            return True 
        else:
            grades.addMessage('FAIL: %s' % self.path)
            grades.addMessage('\t%s' % self.failure)
            grades.addMessage('\tstudent result: "%s"' % result)
            grades.addMessage('\tcorrect result: "%s"' % solutionDict['result'])
            
        return False            

    def writeSolution(self, moduleDict, filePath):
        handle = open(filePath, 'w')
        handle.write('# This is the solution file for %s.\n' % self.path)
        handle.write('# The result of evaluating the test must equal the below when cast to a string.\n')
        
        handle.write('result: "%s"\n' % self.evalCode(moduleDict))
        handle.close()
        return True

