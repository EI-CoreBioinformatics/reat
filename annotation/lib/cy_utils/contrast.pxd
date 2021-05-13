
import numpy as np
cimport numpy as np

cpdef np.ndarray[np.float_t, ndim=6] array_compare(np.ndarray[np.int_t] ref, np.ndarray[np.int_t] pred, np.float_t identity)