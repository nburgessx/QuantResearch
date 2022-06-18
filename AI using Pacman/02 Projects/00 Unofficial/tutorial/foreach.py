# foreach.py
# ----------
# Licensing Information: Please do not distribute or publish solutions to this
# project. You are free to use and extend these projects for educational
# purposes. The Pacman AI projects were developed at UC Berkeley, primarily by
# John DeNero (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# For more info, see http://inst.eecs.berkeley.edu/~cs188/sp09/pacman.html

fruits = ['apples','oranges','pears','bananas']
for fruit in fruits:
    print fruit + ' for sale'
    
fruitPrices = {'apples':2.00, 'oranges': 1.50, 'pears': 1.75}
for fruit, price in fruitPrices.items():
    print '%s cost %f a pound' % (fruit, price)