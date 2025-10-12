import numpy as np

from python import numpy as pnp
import numpy.pybridge

distances: np.ndarray[int,2] = pnp.loadtxt("tests/sequence/data/distances.txt", dtype=pnp.int64)
# Now you can use distances
tree = upgma(distances)