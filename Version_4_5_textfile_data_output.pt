# VERSION 4 

# Outputs particle positions and change in initial angle to a text file
# If particles overlap by more than an allowed amount then they are moved away from one another.
# Can control features of the simulation such as outputting data in the CONTROL PARAMETERS section

from math import sqrt
from math import pi
import turtle
import random
import numpy as np
import time
import csv

wn = turtle.Screen() #sets details of our screen
wn.bgcolor("black")
wn.setup(600,600)
wn.tracer(0)

#system parameters
circle_radius = 123
number_of_particles = 18
particle_diameter = 50
moment_of_inertia = (particle_diameter**2)/8
rotation_directions = [-1, 1]
rotational_acceleration = 2.0
translational_acceleration = 0.1
translational_friction = 0.0001
rotational_friction = (particle_diameter**2)*(translational_friction)/3
time_step = 0.00005
lj_min = 0.0
lj_diameter = 0
run_number = 0
allowed_overlap = 0.0001

# CONTROL PARAMETERS
OVERLAP_TO_TOUCHING = 1
write_to_file_X_Y = 1
write_to_file_overlap = 0

name = "Sim Params:    particle_num = {}     acc_omega = {}    acc_linear = {}    friction_linear = {}".format(number_of_particles, rotational_acceleration,translational_acceleration,translational_friction)

wn.title(name)

from turtle import * #draws our circle region

color('green','green')
begin_fill()
turtle.penup()
turtle.speed(0)
turtle.goto(0,-circle_radius)
circle(circle_radius)
turtle.hideturtle()

end_fill()

balls = []
arrows = []
degrees = np.zeros(number_of_particles)

if write_to_file_overlap == 1:
    with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_OVERLAP.txt','w') as f:
        f.write(f"\n")
    

#PARTICLE INITILIASATIONS
for i in range(number_of_particles): #makes our particles
    balls.append(turtle.Turtle())
    arrows.append(turtle.Turtle())

for ball in balls:

    ball.shape('circle')
    ball.turtlesize(particle_diameter/20,particle_diameter/20)
    ball.color('red')
    ball.penup() #stop line drawing of ball trajectory
    ball.dx = 0
    ball.dy = 0
    ball.cw_acw = random.choice(rotation_directions)
    ball.dtheta = 0
    #puts the balls in the largest square in a circle (not in the whole circle)
    ball.goto((random.random()-0.5)*circle_radius*sqrt(2), (random.random()-0.5)*circle_radius*sqrt(2))

for arrow in arrows:
    arrow.shape('classic')
    arrow.turtlesize(particle_diameter/20,particle_diameter/20)
    arrow.color('black')
    arrow.penup()

#POTENTIAL SETUP
def lennard_jones(position_i,position_j):
    r_ij = position_i - position_j
    mod_r_ij = sqrt(r_ij.dot(r_ij))
    V_LJ = ((particle_diameter + lj_diameter)/mod_r_ij)**12 - ((particle_diameter + lj_diameter)/mod_r_ij)**6
    return V_LJ, r_ij

# checks the lenard jones potentials of all the particles and returns the particle with the highest potential
def initial_check():
    lj = lj_min
    for i in range(0,number_of_particles):  
        for j in range(i+1,number_of_particles):
            r_i = np.array([balls[i].xcor(), balls[i].ycor()])
            r_j = np.array([balls[j].xcor(), balls[j].ycor()])
            if lennard_jones(r_i,r_j)[0] > lj:
                lj = lennard_jones(r_i,r_j)[0]
    for i in range(0,number_of_particles):  
        if ((balls[i].xcor())**2 + (balls[i].ycor())**2 > (circle_radius - particle_diameter/2)**2):
            if lj < 1:
                lj = 1.0          
    return lj

def initialisation_potential():
    vector_LJ = np.array([0.0,0.0])
    #First we reset the velocities so that they STOP moving once they have cleared a ball
    for i in range(0,number_of_particles):  
        balls[i].dx = 0.0
        balls[i].dy = 0.0

    #now we make the velocity of the ball magnitude 1 in the direction away from the balls it overlaps      
    for i in range(0,number_of_particles):  
        for j in range(i+1,number_of_particles):
            #now we set up the two position arrays for ball i and ball
            r_i = np.array([balls[i].xcor(), balls[i].ycor()])
            r_j = np.array([balls[j].xcor(), balls[j].ycor()])
            if lennard_jones(r_i,r_j)[0] > lj_min: #checks Lennard-Jones potential
                        vector_LJ = lennard_jones(r_j,r_i)[0]*lennard_jones(r_j,r_i)[1]                
                        if vector_LJ[0] != 0 and vector_LJ[1] != 0:          
                            vector_LJ = vector_LJ/sqrt(vector_LJ.dot(vector_LJ))
                            balls[j].dx += vector_LJ[0]
                            balls[j].dy += vector_LJ[1]
                            balls[i].dx += -vector_LJ[0]
                            balls[i].dy += -vector_LJ[1]
        # This implements the Lenard Jones potential as if there was a infinite circle of
        # particles with their centres one radius away from the boundary
        if ((balls[i].xcor())**2 + (balls[i].ycor())**2 > (circle_radius - particle_diameter/2)**2):
            r_i = np.array([balls[i].xcor(), balls[i].ycor()])
            r_j = r_i + (particle_diameter/2)*r_i/(r_i.dot(r_i))
            if lennard_jones(r_i,r_j)[0] > lj_min:
                vector_LJ = lennard_jones(r_j,r_i)[0]*lennard_jones(r_j,r_i)[1]                
                if vector_LJ[0] != 0 and vector_LJ[1] != 0:          
                    vector_LJ = vector_LJ/sqrt(vector_LJ.dot(vector_LJ))
                    balls[i].dx += -vector_LJ[0]
                    balls[i].dy += -vector_LJ[1]

    # Makes the magnitude of all vectors one
    for i in range(0,number_of_particles):  
        if (balls[i].dx != 0 and balls[i].dy != 0):    
            dx = balls[i].dx
            dy = balls[i].dy
            new_dx = dx / (sqrt(dx**2 + dy**2))
            new_dy = dy / (sqrt(dx**2 + dy**2))
            balls[i].dx = new_dx
            balls[i].dy = new_dy

    for i in range(0, number_of_particles): 
        balls[i].setx(balls[i].xcor() + balls[i].dx)
        balls[i].sety(balls[i].ycor() + balls[i].dy) 
    
    wn.update()


while initial_check() > 0.0:
    initialisation_potential()
    wn.update

for ball in balls:
    ball.dx = 0
    ball.dy = 0

#SYSTEM DYNAMICS
def move_particles():
    for i in range(0,number_of_particles):

        balls[i].setx(balls[i].xcor() + time_step*balls[i].dx)
        balls[i].sety(balls[i].ycor() + time_step*balls[i].dy)
        balls[i].lt(time_step*balls[i].dtheta)

        
        arrows[i].setx(balls[i].xcor() + time_step*balls[i].dx)
        arrows[i].sety(balls[i].ycor() + time_step*balls[i].dy)
        arrows[i].lt(time_step*balls[i].dtheta)

        # Calculates the change in angle around own axis for each particle
        if write_to_file_X_Y == 1:
            degrees[i] += time_step*balls[i].dtheta
            degrees[i] = degrees[i] % 360

        balls[i].dx +=  translational_acceleration*(random.randint(-1,1)) - translational_friction*balls[i].dx
        balls[i].dy +=  translational_acceleration*(random.randint(-1,1)) - translational_friction*balls[i].dy
        #translational dynamics defined by random walks but with momentum carried over, translational friction also applied
        
        balls[i].dtheta += rotational_acceleration*balls[i].cw_acw - rotational_friction*balls[i].dtheta
        #rotational speed is increased up to a max, i.e each ball has a driving torque and rotational friction

    

#Checks the collision with each of the balls and all the balls not already checked in the first for loop
def collision_check_particles():
    for i in range(0,number_of_particles):  
        for j in range(i+1,number_of_particles):

            if balls[i].distance(balls[j]) <= particle_diameter: #checks for collision between particles
                r_i = np.array([balls[i].xcor(), balls[i].ycor(), 0])
                r_j = np.array([balls[j].xcor(), balls[j].ycor(), 0])


                # FOR PARTICLE i
                r_ij = r_i - r_j
                mod_r_ij = sqrt(r_ij.dot(r_ij))
                if mod_r_ij < particle_diameter:
                    overlap = particle_diameter - mod_r_ij
                    if OVERLAP_TO_TOUCHING == 1:
                        if overlap > allowed_overlap :
                            #print('particles need to repel')
                            balls[i].setx(balls[i].xcor() + 0.5*overlap*r_ij[0]/abs(r_ij[0]))
                            balls[i].sety(balls[i].ycor() + 0.5*overlap*r_ij[1]/abs(r_ij[1]))
                    if write_to_file_overlap == 1:
                        with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_OVERLAP.txt','a') as f:
                            f.write(f"{overlap}\n")

                kappa = 4*moment_of_inertia/(mod_r_ij**2)

                v_i = np.array([balls[i].dx, balls[i].dy, 0])
                v_j = np.array([balls[j].dx, balls[j].dy, 0])
                omega_i = np.array([0, 0, (2*pi/360)*balls[i].dtheta])
                omega_j = np.array([0, 0, (2*pi/360)*balls[j].dtheta])
                v_ij = v_i - v_j - 0.5*np.cross(omega_i + omega_j, r_ij)

                #velocity parallel/perpendicular to the relative position of their centre of masses
                v_ij_parr = r_ij*(v_ij.dot(r_ij))/((mod_r_ij)**2)
                v_ij_perp = v_ij - r_ij*(v_ij.dot(r_ij))/((mod_r_ij)**2)

                #change in linear velocity (from paper)
                delta_v_i = -(v_ij_parr + (kappa/(1+kappa))*v_ij_perp)
                #change in rotational velocity (from paper)
                delta_omega_i = -np.cross(r_ij, delta_v_i)/(2*moment_of_inertia)

                #Adding the change of both velocities to both the balls involved in the collision
                balls[i].dx += delta_v_i[0]
                balls[i].dy += delta_v_i[1]

                balls[i].dtheta += (360/(2*pi))*delta_omega_i[2]


                # NOW FOR PARTICLE j
                r_ji = r_j - r_i
                mod_r_ji = sqrt(r_ji.dot(r_ji))

                if OVERLAP_TO_TOUCHING == 1:
                    if overlap > allowed_overlap :
                            #print('particles need to repel')
                            balls[j].setx(balls[j].xcor() + 0.5*overlap*r_ji[0]/abs(r_ji[0]))
                            balls[j].sety(balls[j].ycor() + 0.5*overlap*r_ji[1]/abs(r_ji[1]))

                kappa = 4*moment_of_inertia/(mod_r_ji**2)

                v_ji = v_j - v_i - 0.5*np.cross(omega_j + omega_i, r_ji)

                #velocity parallel/perpendicular to the relative position of their centre of masses
                v_ji_parr = r_ji*(v_ji.dot(r_ji))/((mod_r_ji)**2)
                v_ji_perp = v_ji - r_ji*(v_ji.dot(r_ji))/((mod_r_ji)**2)

                #change in linear velocity (from paper)
                delta_v_j = -(v_ji_parr + (kappa/(1+kappa))*v_ji_perp)
                #change in rotational velocity (from paper)
                delta_omega_j = -np.cross(r_ji, delta_v_j)/(2*moment_of_inertia)

                balls[j].dx += delta_v_j[0]
                balls[j].dy += delta_v_j[1]

                balls[j].dtheta += (360/(2*pi))*delta_omega_j[2]

#elastic, smooth boundary, need to figure out how to change this
def collision_check_boundary():
    for ball in balls:
        if (ball.xcor())**2 + (ball.ycor())**2 >= (circle_radius - particle_diameter/2)**2:
            ball.dx *= -1
            ball.dy *= -1
        if OVERLAP_TO_TOUCHING == 1:
            if (ball.xcor())**2 + (ball.ycor())**2 >= (circle_radius - particle_diameter/2)**2 + allowed_overlap:
                boundary_overlap = sqrt((ball.xcor())**2 + (ball.ycor())**2 - (circle_radius - particle_diameter/2)**2)
                ball.setx(ball.xcor() - 0.5*boundary_overlap*ball.xcor()/abs(ball.xcor()))
                ball.sety(ball.ycor() - 0.5*boundary_overlap*ball.ycor()/abs(ball.ycor()))

# Creates the text files that the data will be written to
if write_to_file_X_Y == 1:
    with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_X_pos.txt','w') as f:
        f.write(f"{name}\n{number_of_particles}\n{circle_radius}\n")
    with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_Y_pos.txt','w') as f:
        f.write(f"{name}\n{number_of_particles}\n{circle_radius}\n")
    with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_DEGREES_pos.txt','w') as f:
        f.write(f"{name}\n{number_of_particles}\n{circle_radius}\n")

#MAIN LOOP
while True:
    move_particles()
    collision_check_particles()
    collision_check_boundary()
    
    #if run_number % (1/time_step) == 0:
    wn.update() #updates the window

# Writes to the text files
    if run_number % (0.1/time_step) == 0:
        if write_to_file_X_Y == 1:
            with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_X_pos.txt','a') as f:
                for i in range(0,number_of_particles):
                    f.write(f"{balls[i].xcor()}\n")
            with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_Y_pos.txt','a') as f:
                for i in range(0,number_of_particles):
                    f.write(f"{balls[i].ycor()}\n")
            with open('X:/Physics/StudentProjects/Final Year/2021-22/robots-md/22112021 FYP---Simulating-Spinning-Robots-master/particles_DEGREES_pos.txt','a') as f:
                for i in range(0,number_of_particles):
                    f.write(f"{degrees[i]}\n")

    run_number += 1
    
wn.mainloop() #keeps window open
