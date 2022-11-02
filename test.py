from GeneSequencing import *

def test_should_solve_mismatchedLen():
  solver = GeneSequencing()
  answer = solver.align("hi", "hello", False, 100)
  
  cell = solver.costDict[tuple((0,0))]
  assert(cell.costVal == 0)
  assert(cell.prev == None)
  assert(cell.direction == None)
  
  cell = solver.costDict[tuple((0,1))]
  assert(cell.costVal == 5)
  assert(cell.prev == tuple((0,0)))
  assert(cell.direction == Direction.LEFT)
  
  cell = solver.costDict[tuple((0,2))]
  assert(cell.costVal == 10)
  assert(cell.prev == tuple((0,1)))
  assert(cell.direction == Direction.LEFT)

  cell = solver.costDict[tuple((1,0))]
  assert(cell.costVal == 5)
  assert(cell.prev == tuple((0,0)))
  assert(cell.direction == Direction.TOP)

  cell = solver.costDict[tuple((1,1))]
  assert(cell.costVal == -3)
  assert(cell.prev == tuple((0,0)))
  assert(cell.direction == Direction.DIAGONAL)

  cell = solver.costDict[tuple((1,2))]
  assert(cell.costVal == 2)
  assert(cell.prev == tuple((1,1)))
  assert(cell.direction == Direction.LEFT)

  #Final Cell
  cell = solver.costDict[tuple((5,2))]
  assert(cell.costVal == 13)
  assert(cell.prev == tuple((4,2)))
  assert(cell.direction == Direction.TOP)

  assert(answer['align_cost'] == 13)
  assert(answer['seqi_first100'] == 'hi---')
  assert(answer['seqj_first100'] == 'hello')
  