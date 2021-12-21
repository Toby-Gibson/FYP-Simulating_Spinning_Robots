# Uses MatPlotLib to visualise the data
# The visualisation is outputted as frames (we are yet to animate these frames into video format
# Similarly we need to add arrows to the particles and add their rotation
# This code uses the text file data outputted using Version_4_5_textfile_data_output.py

from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
import numpy as np
import math

#Extracting the number of particles used in the simulation

file_X_data = open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_X_pos.txt')
next(file_X_data)
with file_X_data as f:
    number_of_particles = int(f.readline())
    circle_radius = int(f.readline())

f.close()

print(circle_radius)

#READING in the X data values
particle_data = []

file_X_data = open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_X_pos.txt')
next(file_X_data)
next(file_X_data)
next(file_X_data)

with file_X_data as f:
    for line in f.readlines():
        particle_data.append(float(line))

particle_X_data = np.array(particle_data)

f.close()


#READING in the Y data values
particle_data = []

file_Y_data = open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_Y_pos.txt')
next(file_Y_data)
next(file_Y_data)
next(file_Y_data)

with file_Y_data as f:
    for line in f.readlines():
        particle_data.append(float(line))

particle_Y_data = np.array(particle_data)

f.close()

# Extracting key values  from the data 
Num_of_Data_Points = int(len(particle_X_data))
Time_Steps = int(float(Num_of_Data_Points) / float(number_of_particles))

# Plotting the graphs for each timestep
x = np.zeros(number_of_particles)
y = np.zeros(number_of_particles)

for j in range(0,Time_Steps):
    for i in range(0,number_of_particles):
        x[i] = particle_X_data[number_of_particles*j+i]
        y[i] = particle_Y_data[number_of_particles*j+i]
    

    
    figure, axes = plt.subplots(figsize=[2.5,2.5], dpi = 100)
    plt.xlim([-circle_radius, circle_radius])
    plt.ylim([-circle_radius, circle_radius])
    draw_circle = plt.Circle((0, 0), float(circle_radius),fill=False)
    plt.scatter(x, y, s = 314.15)
    
    #axes.set_aspect(1)
    axes.add_artist(draw_circle)
    plt.show()
