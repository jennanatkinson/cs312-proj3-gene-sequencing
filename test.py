from GeneSequencing import *

def test_should_solve_mismatchedLen():
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

  cell = solver.costDict[tuple((1,0))]
  assert(cell.costVal == 5)
  assert(cell.prev == tuple((0,0)))
  assert(cell.direction == Direction.DOWN)

  cell = solver.costDict[tuple((1,1))]
  assert(cell.costVal == -3)
  assert(cell.prev == tuple((0,0)))
  assert(cell.direction == Direction.DIAGONAL)

  cell = solver.costDict[tuple((1,2))]
  assert(cell.costVal == 2)
  assert(cell.prev == tuple((1,1)))
  assert(cell.direction == Direction.RIGHT)

  #Final Cell
  cell = solver.costDict[tuple((5,2))]
  assert(cell.costVal == 13)
  assert(cell.prev == tuple((4,1)))
  assert(cell.direction == Direction.DIAGONAL)
  