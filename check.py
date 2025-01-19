import h5py
import numpy as np
import matplotlib.pyplot as plt

file_name = f'ics.0.hdf5'
f = h5py.File(file_name, 'r')
all_coordinates = f['PartType0']['Coordinates'][:]

for ii in range(15):
    file_name = f'ics.{ii+1}.hdf5'
    f = h5py.File(file_name, 'r')
    coordinates = f['PartType0']['Coordinates'][:]

    all_coordinates=np.concatenate((all_coordinates, coordinates))


plt.hist(all_coordinates[:,0], bins=30, edgecolor='black')
plt.xlabel('x')
plt.savefig('histogram_x.png')
plt.close()

plt.hist(all_coordinates[:,1], bins=30, edgecolor='black')
plt.xlabel('y')
plt.savefig('histogram_y.png')
plt.close()

plt.hist(all_coordinates[:,2], bins=30, edgecolor='black')
plt.xlabel('x')
plt.savefig('histogram_z.png')
plt.close()

plt.scatter(all_coordinates[0:10000,0],all_coordinates[0:10000,1], s=0.1)
plt.xlabel('x')
plt.xlabel('y')
plt.savefig('density.png')
plt.close()

