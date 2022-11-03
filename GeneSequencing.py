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

# Used to compute the bandwidth for banded version
MAXINDELS = 3
SEQRETURNLEN = 100

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

# Used when building up the matrix of which direction the previous cell came from
class Direction(Enum):
	LEFT = auto()
	TOP = auto()
	DIAGONAL = auto()
	
# Contains the costVal, previous cell which helped calculate that, and the direction the previous cell came from
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

	# Print the cost dictionary 
	def printDict(self, seq1, seq2):
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
					item = self.costDict.get(tuple((i, j)))
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
		print()

	# Traverse through the matrix and build up the alignment strings
	def getAlignmentStrings(self, seq1, seq2, cell:tuple):
		alignment1, alignment2 = "", ""
		while not (cell[0] == 0 and cell[1] == 0):
			cost = self.costDict.get(cell)
			if cost.direction == Direction.LEFT:
				alignment1 += seq1[cell[1]-1]
				alignment2 += "-"
			elif cost.direction == Direction.TOP:
				alignment1 += "-"
				alignment2 += seq2[cell[0]-1]
			elif cost.direction == Direction.DIAGONAL:
				alignment1 += seq1[cell[1]-1]
				alignment2 += seq2[cell[0]-1]			
			cell = cost.prev

		#Reverse strings and cut to 100
		return alignment1[::-1][:100], alignment2[::-1][:100]

	# Check the banded setting and return the value
	def checkBanded(self, banded:bool, trueVal, falseVal):
		if banded:
			return trueVal
		else:
			return falseVal
		
	# Given two sequences, find the optimal alignment of Insert/Deletion, Substitution or Match
	#    @required _seq1_(top of matrix) and _seq2_(side of matrix) are two sequences to be aligned
	#    @required _banded_ is a boolean for computing banded alignment or full alignment
	# 	 @required align_length_ = how many base pairs to use in computing the alignment
	def align(self, seq1, seq2, banded, align_length):
		if len(seq1) > align_length:
			seq1 = seq1[:align_length]
		if len(seq2) > align_length:
			seq2 = seq2[:align_length]
		if seq1 == seq2:
			return {'align_cost':len(seq1)*MATCH, 'seqi_first100':seq1, 'seqj_first100':seq2}
		
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.numColumns = len(seq1)+1
		self.numRows = len(seq2)+1
		self.costDict = dict() # Key-Coordinate[row, column] : Value-[Cost, PrevCoordinate]
		
		self.costDict[tuple((0,0))] = Cost(0, None, None)
	
		# Fill first row
		numFirstRow = self.checkBanded(banded, min(1+MAXINDELS,self.numColumns), self.numColumns)
		for i in range(1, numFirstRow):
			prevCell = tuple((0, i-1))
			self.costDict[tuple((0,i))] = Cost(self.costDict.get(prevCell).costVal+INDEL, prevCell, Direction.LEFT)

		#Fill first col
		numFirstCol = self.checkBanded(banded, min(1+MAXINDELS,self.numRows), self.numRows)
		for i in range(1, numFirstCol):
			prevCell = tuple((i-1, 0))
			self.costDict[tuple((i,0))] = Cost(self.costDict.get(prevCell).costVal+INDEL, prevCell, Direction.TOP)

		# self.printDict(seq1, seq2)

		#Run algorithm for all cells
		for rowIndex in range(1, self.numRows):
			colStart = self.checkBanded(banded, max(rowIndex-MAXINDELS,1), 1)
			colEnd = self.checkBanded(banded, min(rowIndex+MAXINDELS+1,self.numColumns), self.numColumns)
			# Calculate costs with priority: left > top > diagonal
			for colIndex in range(colStart, colEnd):
				minCost = Cost(math.inf, None, None)

				# Calculate diagonal cost if match/sub
				diagonalPrevCell = tuple((rowIndex-1, colIndex-1))
				diagonalPrevCost = self.costDict.get(diagonalPrevCell)
				if diagonalPrevCost is not None:
					diagonalCost = Cost(diagonalPrevCost.costVal, diagonalPrevCell, Direction.DIAGONAL)
					if (seq1[colIndex-1] == seq2[rowIndex-1]):
						diagonalCost.costVal += MATCH
					else:
						diagonalCost.costVal += SUB
					if diagonalCost.costVal < minCost.costVal:
						minCost = diagonalCost

				# Calculate top cost if indel
				topPrevCell = tuple((rowIndex-1, colIndex))
				topPrevCost = self.costDict.get(topPrevCell)
				if topPrevCost is not None:
					topCost = Cost(topPrevCost.costVal+INDEL, topPrevCell, Direction.TOP)
					if topCost.costVal <= minCost.costVal:
						minCost = topCost

				# Calculate left cost if indel
				leftPrevCell = tuple((rowIndex, colIndex-1))
				leftPrevCost = self.costDict.get(leftPrevCell)
				if leftPrevCost is not None:
					leftCost = Cost(leftPrevCost.costVal+INDEL, leftPrevCell, Direction.LEFT)
					if leftCost.costVal <= minCost.costVal:
						minCost = leftCost
				
				self.costDict[tuple((rowIndex,colIndex))] = minCost
				# self.printDict(seq1, seq2)
		
		cell = tuple((len(seq2), len(seq1)))
		cost = self.costDict.get(cell)
		score = math.inf
		if cost is not None: # This would happen if it is outside of the band
			score = cost.costVal
		if score != math.inf:
			alignment1, alignment2 = self.getAlignmentStrings(seq1, seq2, cell)
		else:
			alignment1 = alignment2 = "No Alignment Possible"
		# print(alignment1)
		# print(alignment2)

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		# score = random.random()*100
		# alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq1), align_length, ',BANDED' if banded else '')
		# alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}


