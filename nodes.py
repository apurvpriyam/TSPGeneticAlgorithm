import numpy as np

class Node(object):
	def __init__(self, i, x, y):
		self.id = i # id of the ndoe
		self.x = x # x coordinate of node
		self.y = y # y coordinate of node

	def distance(self, node):
		"""
		Calculates distance between two nodes. Input:
		node: node from which the distance needs to be calculated
		"""
		xDis = abs(self.x - node.x)
		yDis = abs(self.y - node.y)
		return(np.sqrt((xDis ** 2) + (yDis ** 2)))