import magritte.core as magritte

import sys

print(magritte.pcmp_comm_rank(), magritte.pcmp_length(int(sys.argv[1])))