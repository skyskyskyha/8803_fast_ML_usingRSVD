from scipy.sparse import csr_matrix
import numpy as np
import time
from PIL import Image
from fSVT import fastSVT_U  # Replace with actual import

# Load and preprocess the image
pic = Image.open('new4.jpg')
pic = np.array(pic) / 255

# Load 'pic_Omega' (you need to convert this data to Python format or load it if it's a .mat file)
# Assuming 'Ome' is defined in pic_Omega
# Omega = full([Ome; Ome; Ome])  # Convert this line to Python after loading pic_Omega

# Convert pic to sparse matrix
M = csr_matrix(np.vstack((pic[:, :, 0], pic[:, :, 1], pic[:, :, 2])))

# Commented out since parfor is not applicable in Python directly
# for i in range(1):
#     pass

# Measure time
t = time.time()
print(M)
# Call the custom function 'fastSVT_U' (you need to define or translate this function)
[X, iters, k] = fastSVT_U(M, 0.1, 0, 1, 50, 10)

# Calculate elapsed time
t_SVT = time.time() - t

# Calculate error
err = np.sum(np.abs(np.vstack((pic[:, :, 0], pic[:, :, 1], pic[:, :, 2])) - X)) / (2048 * 2048 * 3) * 255
