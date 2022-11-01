from GeneSequencing import *

def test_print():
  solver = GeneSequencing()
  solver.align("hi", "hello", False, 0)
  
  cell = solver.costDict[tuple((0,0))]
  assert(cell.costVal == 0)
  assert(cell.prev == None)
  assert(cell.direction == None)
  
  cell = solver.costDict[tuple((0,1))]
  assert(cell.costVal == 5)
  assert(cell.prev == tuple((0,0)))
  assert(cell.direction == Direction.RIGHT)
  
  cell = solver.costDict[tuple((0,2))]
  assert(cell.costVal == 10)
  assert(cell.prev == tuple((0,1)))
  assert(cell.direction == Direction.RIGHT)