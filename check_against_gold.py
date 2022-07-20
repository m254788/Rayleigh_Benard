import numpy as np
import sys
gs_name = sys.argv[1]
check_name = sys.argv[2]

# ly51_T10000:
# ly = 51
# aspect_ratio = 2
# Ra = 100000
# Thot = 1.
# Tcold = 0
# maxT = 10000

#GOLDSTANDARD is a txt file
GOLD_STANDARD = np.loadtxt(gs_name,dtype=np.float64)
#checkboard is a binary file
check_board = np.fromfile(check_name,dtype=np.float64)


dif = check_board-GOLD_STANDARD
print(np.linalg.norm(dif,ord=1)/np.linalg.norm(GOLD_STANDARD,ord=1))	
