# unpack.py
# ---------
# Licensing Information: Please do not distribute or publish solutions to this
# project. You are free to use and extend these projects for educational
# purposes. The Pacman AI projects were developed at UC Berkeley, primarily by
# John DeNero (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# For more info, see http://inst.eecs.berkeley.edu/~cs188/sp09/pacman.html

import os, cPickle, sys

if len(sys.argv) != 3: 
  print 'Usage: %s stats_file team_name' % sys.argv[0]
  print 'Unpacks the stats file of a server into a bunch of replay files.'
  if len(sys.argv) == 2: 
    d = cPickle.load(open(sys.argv[1]))
    print 'Team names:', d.keys()
  sys.exit(2)

d = cPickle.load(open(sys.argv[1]))
user = sys.argv[2]
k = 0
print 'Unpacking games for', user
for g, w in d[user]['gameHistory']:
    k += 1
    t = {'layout': g.state.data.layout, 'agents': g.agents, 'actions': g.moveHistory, 'length': g.length}
    fname = 'replay_' + user + '_' + str(k)
    print 'Game:', fname
    cPickle.dump(t,file(fname, 'w'))