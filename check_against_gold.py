import numpy as np
import sys
gs_name = sys.argv[1]
check_name = sys.argv[2]

GOLD_STANDARD = np.loadtxt(gs_name,dtype=np.float64)
check_board = np.loadtxt(check_name,dtype=np.float64)


dif = check_board-GOLD_STANDARD
print(np.linalg.norm(dif,ord=1)/np.sum(GOLD_STANDARD))	
