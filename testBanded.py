
from GeneSequencing import GeneSequencing


def main():
  seq1 = "XXXXXXX"
  seq2 = "XX"
  numRows = len(seq2) + 1
  numColumns = len(seq1) + 1
  BAND = 3
  table_data = [[]]
  # Put the first row aka listing the seq1
  for i in range(-1, numColumns):
    if i == -1:
      table_data[0].append(" ")
    elif i == 0:
      table_data[0].append("_")
    else:
      table_data[0].append(seq1[i-1])
  
  for i in range(0, numRows):
    table_data.append([])
    last = len(table_data)-1
    for j in range(-1, numColumns):
      if j == -1 and i == 0:
        table_data[last].append("_")
      elif j == -1:
        table_data[last].append(seq2[i-1])
      elif j >= i-BAND and j <= i+BAND:
        table_data[last].append("*")
      else:
        table_data[last].append(" ")

  formatString = "{: >5} "
  for row in table_data:
    for item in row:
      print(formatString.format(item), end="")
    print()
  print()

def check_align():
  solver = GeneSequencing()
  answer = solver.align("XXXXXXX", "XXXXXXXX", True, 100)

#main()
check_align()
