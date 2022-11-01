#!/usr/bin/python3

# from os import SCHED_OTHER
from enum import Enum, auto
from re import X
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

# minVal = 0
# maxVal = 1

# class Cell:
# 	def __init__(self, x:int, y:int):
# 		self.x = x
# 		self.y = y

# 	def __hash__(self):
# 		return (self.y * (maxVal - minVal)) + self.x

	# def __eq__(self, other):
	# 	return (self.x, self.y) == (other.x, other.y)

class Direction(Enum):
	RIGHT = auto()
	DOWN = auto()
	DIAGONAL = auto()
	
class Cost:
	def __init__(self, cost:int, prev:tuple, dir:Direction):
		self.costVal = cost
		self.prev = prev
		self.direction = dir

	def __eq__(self, other):
		return (self.costVal, self.prev, self.direction) == (other.costVal, other.prev, self.direction)


class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def printDict(self, dict, seq1, seq2):
		table_data = [[]]
		# Put the first row aka listing the seq1
		for i in range(-1, self.numColumns):
			if i == -1:
				table_data[0].append(" ")
			elif i == 0:
				table_data[0].append("_")
			else:
				table_data[0].append(seq1[i-1])
		
		for i in range(0, self.numRows):
			table_data.append([])
			last = len(table_data)-1
			for j in range(-1, self.numColumns):
				if j == -1 and i == 0:
					table_data[last].append("_")
				elif j == -1:
					table_data[last].append(seq2[i-1])
				else:
					string = "" #f"({i},{j}):"
					item = dict.get(tuple((i, j)))
					if type(item) == Cost:
						string += f"{item.costVal}"
					else:
						string += " "
					table_data[last].append(string)

		formatString = "{: >5} "
		for row in table_data:
			for item in row:
				print(formatString.format(item), end="")
			print()


	def findBestCost(self, costs):
		bestCost = Cost(None, None, None)
		for cost in costs:
			if bestCost.costVal is None:
				bestCost = cost
			elif cost.costVal < bestCost.costVal:
				bestCost = cost
		return bestCost

	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.numColumns = len(seq1)+1
		self.numRows = len(seq2)+1
		#matrix = seq1 * seq2
		#0 row and 0 column == _
		self.costDict = dict() # Key-Coordinate[row, column] : Value-[Cost, PrevCoordinate]
		
		self.costDict[tuple((0,0))] = Cost(0, None, None)
		# Fill first row
		for i in range(1, self.numColumns):
			prevCell = tuple((0, i-1))
			prevCost = self.costDict.get(prevCell)
			self.costDict[tuple((0,i))] = Cost(prevCost.costVal+INDEL, prevCell, Direction.RIGHT)

		#Fill first col
		for i in range(1, self.numRows):
			prevCell = tuple((i-1, 0))
			prevCost = self.costDict.get(prevCell)
			self.costDict[tuple((i,0))] = Cost(prevCost.costVal+INDEL, prevCell, Direction.DOWN)

		self.printDict(self.costDict, seq1, seq2)

		#Run algorithm for all cells
		for rowIndex in range(1, self.numRows):
			for colIndex in range(1, self.numColumns):
				rightPrevCell = tuple((rowIndex, colIndex-1))
				rightCost = Cost(self.costDict.get(rightPrevCell).costVal+INDEL, rightPrevCell, Direction.RIGHT)
				
				downPrevCell = tuple((rowIndex-1, colIndex))
				downCost = Cost(self.costDict.get(downPrevCell).costVal+INDEL, downPrevCell, Direction.DOWN)
				
				# Calculate diagonal cost if match/sub
				diagonalPrevCell = tuple((rowIndex-1, colIndex-1))
				diagonalCost = Cost(self.costDict.get(diagonalPrevCell).costVal, diagonalPrevCell, Direction.DIAGONAL)
				if (seq1[colIndex-1] == seq2[rowIndex-1]):
					diagonalCost.costVal += MATCH
				else:
					diagonalCost.costVal += SUB
				
				self.costDict[tuple((rowIndex,colIndex))] = self.findBestCost([diagonalCost, rightCost, downCost])
				self.printDict(self.costDict, seq1, seq2)
		
		# cell = tuple((seq1, seq2))
		# finalCost = costDict.get(cell)
		# score = finalCost.costVal

		# while not (cell[0] == 0 and cell[1] == 0):
		

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = random.random()*100
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}


