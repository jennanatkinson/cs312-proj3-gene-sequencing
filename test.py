from GeneSequencing import *

def test_should_solve_simpleMatch():
  solver = GeneSequencing()
  answer = solver.align("h", "h", False, 100)
  assert(answer['align_cost'] == MATCH)
  assert(answer['seqi_first100'] == 'h')
  assert(answer['seqj_first100'] == 'h')

def test_should_solve_simpleSub():
  solver = GeneSequencing()
  answer = solver.align("h", "a", False, 100)
  assert(answer['align_cost'] == SUB)
  assert(answer['seqi_first100'] == 'h')
  assert(answer['seqj_first100'] == 'a')

def test_should_solve_simpleIndel():
  solver = GeneSequencing()
  answer = solver.align("h", "", False, 100)
  assert(answer['align_cost'] == INDEL)
  assert(answer['seqi_first100'] == 'h')
  assert(answer['seqj_first100'] == '-')

def test_should_solve_sameSeq():
  solver = GeneSequencing()
  seq = "match"
  answer = solver.align(seq, seq, False, 100)
  assert(answer['align_cost'] == MATCH*len(seq))
  assert(answer['seqi_first100'] == seq)
  assert(answer['seqj_first100'] == seq)

def test_should_solve_closeAlignedSeq():
  solver = GeneSequencing()
  answer = solver.align("polynomial", "exponential", False, 100)
  assert(answer['align_cost'] == -1)

def test_should_solve_largerSeq1Len():
  solver = GeneSequencing()
  answer = solver.align("hello", "hi", False, 100)
  assert(answer['align_cost'] == 13)
  assert(answer['seqi_first100'] == 'hello')
  assert(answer['seqj_first100'] == 'hi---')

def test_should_solve_largerSeq2Len():
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
  