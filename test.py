import h5py
import numpy as np
import matplotlib.pyplot as plt

file_name = f'./ics_000.hdf5'
print(file_name)
f = h5py.File(file_name, 'r')
all_coordinates = f['PartType0']['Coordinates'][:]

print(np.max(all_coordinates))
print(np.min(all_coordinates))



plt.hist(all_coordinates[:,0], bins=30, edgecolor='black')
plt.xlabel('x')
plt.savefig('histogram_x11.png')
plt.close()

plt.hist(all_coordinates[:,1], bins=30, edgecolor='black')
plt.xlabel('y')
plt.savefig('histogram_y11.png')
plt.close()

plt.hist(all_coordinates[:,2], bins=30, edgecolor='black')
plt.xlabel('z')
plt.savefig('histogram_z11.png')
plt.close()

plt.scatter(all_coordinates[0:10000,0],all_coordinates[0:10000,1], s=0.1)
plt.xlabel('x')
plt.xlabel('y')
plt.savefig('density.png')
plt.close()

