import numpy as np
import scipy as sci
from scipy import linalg
import matplotlib.pyplot as plt
import statistics
import timeit

mean_X_velocity = 2000
standard_deviation = 50
mean_time = 0.0025
n = 300
Tracking_station_positions =(30,35,40,45)
station_resolution=0.01

#generates a velocity from a normal distribution
def velocity(mean_X_velocity, standard_deviation):
    
    #z-component velocity generated from normal distribution
    velocity_z = np.random.normal(mean_X_velocity, standard_deviation)
    #velocity vector
    velocity_vector = np.array([0,0, velocity_z])
    
    return velocity_z,velocity_vector
#generates a decay displacement from a normal distribution and a decay time from a exponential distribution
def position(mean_X_velocity, standard_deviation, mean_time):
    
    velocity_z, velocity_vector = velocity(mean_X_velocity, standard_deviation)
    #decay time from exponential dsitribution using the scale parameter which is the inverse of the rate parameter
    decay_time = np.random.exponential((mean_time))
    #decay vector based on displacement of the particle, calcaulted from the decay_time and the speed
    displacement = decay_time*velocity_z
    decay_vertex = [0,0,displacement]
    
    return decay_time, displacement, decay_vertex
#the beam is not relativsitc but the daughter particle is meaning there will be a lorenz contraction acting on the direction of travel
def Lorenz_gamma(z):
    #the contraction co-efficient
    gamma = np.sqrt(1- ((10**7)**2 / ((3*10**8)**2)))
    
    return z*gamma

#################################################################################################################################################################################################
                                                               # PART A - GENERATING DISTRIBUTIONS
#################################################################################################################################################################################################

def Test(mean_X_velocity,standard_deviation,mean_time):
    
    #calling decay vertex and time from position function
    decay_time, displacement, decay_vertex = position(mean_X_velocity, standard_deviation, mean_time)
    velocity_z, velocity_vector = velocity(mean_X_velocity, standard_deviation)
    
    print('velocity_z =', velocity_z)
    print('velocity_vector =' , velocity_vector)
    print('decay_time =', decay_time)
    print('decay_vertex =', decay_vertex)
    print('displacement =', decay_vertex)
    
    return

#print(Test(mean_X_velocity, standard_deviation, mean_time))

#This function creates a distribution list for decay time, decay vertex
def distribution_list(mean_X_velocity,standard_deviation,mean_time,n):
    
    #creating empty lists for decay time, vertex and velocity of particle
    velocity_list=[]
    decay_time_list=[]
    displacement_list=[]
    
    
    for i in range(n):
        #calling decay time, vertex and velocity of the particle
        decay_time, displacement, decay_vertex = position(mean_X_velocity, standard_deviation, mean_time)
        velocity_z, velocity_vector = velocity(mean_X_velocity, standard_deviation)
        
        #appending different values for each variable to a list
        velocity_list.append(velocity_z)
        decay_time_list.append(decay_time)
        displacement_list.append(displacement)
        
    x_velocity = velocity_list
    plt.hist(x_velocity, bins = 50)
    plt.title('Velocity distribution')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Number of particles')
    plt.show()
        
    return velocity_list,decay_time_list,displacement_list

#This function plots a histogram for the generated velocity distribution
def distribution_test_velocity(mean_X_velocity,standard_deviation,mean_time,n):
    #calling the velocity list from the distribution list function
    velocity_z,decay_time,displacement = distribution_list(mean_X_velocity, standard_deviation, mean_time, n)
    
    x_velocity = velocity_z
    plt.hist(x_velocity, bins = 50)
    plt.title('Velocity distribution')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Number of particles')
    plt.show()

    return
#This function plots a histogram for the generated decay_time distribution
def distribution_test_decay_time(mean_X_velocity,standard_deviation,mean_time,n):
    #calling the decay time list from the distribution list function
    velocity_z,decay_time,displacement = distribution_list(mean_X_velocity, standard_deviation, mean_time, n)

    x_decay_time = decay_time
    plt.hist(x_decay_time, bins = 50)
    plt.title('Decay time distribution')
    plt.xlabel('Decay time (/s)')
    plt.ylabel('Number of particles')
    plt.show()

    
    return
#This function plots a histogram for the generated decay vertex distribution
def distribution_test_decay_vertex(mean_X_velocity,standard_deviation,mean_time,n):
    #calling the decay vertex list from the distribution list function
    velocity_z,decay_time,displacement = distribution_list(mean_X_velocity, standard_deviation, mean_time, n)

    x_displacement = displacement
    plt.hist(x_displacement, bins = 50)
    plt.title('Displacement distribution')
    plt.xlabel('Displacement (m)')
    plt.ylabel('Number of Particles')
    plt.show()
    
    return

distribution_test_velocity
distribution_list(mean_X_velocity, standard_deviation, mean_time, n)
################################################################################################################################################################################################
                                                            # PART B GENERATING DECAY ANGLES
#################################################################################################################################################################################################
###########################################################Analytical method#################################################################################################################
#This function is the Analytical method for producing points over the solid angle
def Analytical_method():
    #assuming unit sphere rho is set to 1
    r = 1
    #phi is already uniformaly distributed between 2pi and pi
    phi = 2 * np.pi * np.random.random()
    #using the du,dv method theta is equal to inverse cosine of u
    theta = np.arccos(1-(2*np.random.random()))
    
    # convert to cartesian coordinates
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    return np.array([x, y, z])
#This function plots the theta points generated in the Analytical method to show that it follows a sin distribution
def Theta_histogram():
    
    def Theta_function():
        theta = np.arccos(1-(2*np.random.random()))
        return theta
    
    theta_list =  [Theta_function() for i in range(50000)]
    
    hist1, bins1, patches1 = plt.hist(theta_list, bins=50, density=True, label="Analytic Method")
    bin_centres = (bins1[1:] + bins1[:-1])/2
    plt.plot(bin_centres, np.sin(bin_centres)/2, label=r'$sin(\theta)$')
    
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$P(\theta)$')
    plt.legend()
    plt.show()
    
    return
#This function plots the theta points generated in the Analytical method to show that it follows a sin distribution
def Phi_histogram():
    
    def Phi_function():
        phi = 2 * np.pi * np.random.random()
        return phi

    phi_list =  [Phi_function() for i in range(50000)]
    
    hist1, bins1, patches1 = plt.hist(phi_list, bins=50, density=True, label="Analytic Method")
    bin_centres = (bins1[1:] + bins1[:-1])/2
    plt.plot(bin_centres, 2* np.pi * bin_centres/2, label=r'$sin(\theta)$')
    
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$P(\theta)$')
    plt.legend()
    plt.show()
    
    return
#Plotting the points outputted by the analytical method to show that it creates a spherical shell of co-ordinates
def Analytical_method_plot():
    
    co_ords =np.array([Analytical_method() for i in range(1000)])
    x,y,z = co_ords[:,0], co_ords[:,1], co_ords[:,2]
    
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
    angle = 45
    ax.set_zlabel('Z-Axis', rotation=angle)
    ax.set_ylabel('Y-Axis')
    ax.set_xlabel('X-Axis')
    
    ax.scatter3D(x, y, z, color = "green")
    
    plt.show()
    
    return
###########################################################Analytical method####################################################################################################################
###########################################################Accept Reject method#################################################################################################################
#This function is the accept/reject method for calculating a uniform distribution of points
def accept_reject():
    #generating a random number of points -1 and 2
    x = 2 * np.random.random() - 1
    y = 2 * np.random.random() - 1
    z = 2 * np.random.random() - 1
    
    #the bounds closely approaching 0 to try and mimic that of a infinitesimaly thin shell since it cannot compute for one
    while np.sqrt(x**2+ y**2+ z**2) >1.0000000000000000005 and np.sqrt(x**2+ y**2+ z**2) < 0.9999999999999995:
        # only accepting points that lie on the surface of a sphere wuth radius 1
        x = 2 * np.random.random() - 1
        y = 2 * np.random.random() - 1
        z = 2 * np.random.random() - 1
        
    return np.array([x, y,z])
# creating a list to plot the distribution of accept reject points
def accept_reject_list():
    
    co_ords = np.array([accept_reject() for i in range(1)])
    x,y,z = (co_ords[:,0], co_ords[:,1], co_ords[:,2])
    
    return np.array([x,y,z])
#plotting the accept/reject list of points
def accept_reject_plot():
    
    co_ords = np.array([accept_reject() for i in range(1000)])
    x,y,z = (co_ords[:,0], co_ords[:,1], co_ords[:,2])
    
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection='3d')
    angle = 45
    ax.set_zlabel('Z-Axis', rotation=angle)
    ax.set_ylabel('Y-Axis')
    ax.set_xlabel('X-Axis')
    
    ax.scatter3D(x,y,z, c='g')
    
    plt.show()
    
    return

###########################################################Accept Reject method#################################################################################################################
############################################################## Co-ordinate #####################################################################################################################
#selecting co-ordinates from a list of selected points made by the analytical method
def Co_ord_selection():
    
   decay_time, displacement, decay_vertex = position(mean_X_velocity, standard_deviation, mean_time)
   co_ords =np.array([Analytical_method() for i in range(1000)])
   x_vals,y_vals,z_vals = co_ords[:,0], co_ords[:,1], co_ords[:,2]
   
   #since the daughter particle is moving at relativistic speeds the length is contracted
   x = Lorenz_gamma(np.random.choice(x_vals))
   y = Lorenz_gamma(np.random.choice(y_vals))
   #from lab frame it has moved an additional length equal to what the beam particle has moved before decaying
   z = Lorenz_gamma(np.random.choice(z_vals))
   
   return np.array([x,y,z])

#############################################################Printing##########################################################################################################################
start_analytical_time = timeit.timeit()
Analytical_method_plot()
end_analytical_time = timeit.timeit()
#print(end_analytical_time - start_analytical_time)

start_accept_reject_time = timeit.timeit()
accept_reject_plot()
end_accept_reject_time = timeit.timeit()
#print(end_accept_reject_time - start_accept_reject_time)

#print (Co_ord_selection())
#############################################################Printing##########################################################################################################################
#################################################################################################################################################################################################
                                                            # PART C PROPAGATION PARAMTERES
#################################################################################################################################################################################################
#Function to find the co-efficients for the striaght line equation of the particle
def Tracking_parameters():
    #taking the decay_vertex parameter by calling the position function
    decay_time, displacement, dv = position(mean_X_velocity, standard_deviation, mean_time)
    #position vector outputs the decay vertex as a vector so have to unpack it by calling the z element
    decay_vertex = dv[2]
    # choosing random x and y co-ordinates from the co-ord selection function, to find the angle the particle is projected at
    x,y,z = Co_ord_selection()
    #generating the straight line eqaution co-efficients using the co-ordinates selected
    #the co-ordinates are from the lab frame
    mx = (x/(z+(decay_vertex)))
    my = (y/(z+(decay_vertex)))
    cx = -(decay_vertex) * mx
    cy = -(decay_vertex) * my
    #using d to eliminate particles goin backwards, if z is negative
    if (z) >= displacement:
        d = 1
    else:
        d = -1   
        
    return mx,my,cx,cy,d
#Creating the straight line equation from the co-efficients
def Tracking_line_equations(z):
    
    mx,my,cx,cy,d = Tracking_parameters()
    
    #using the paramteres to determine the striaght line equations to calcaulte future x and y co-ordinates
    y = my*z + cy
    x = mx*z + cx
    
    return x,y,z

#Using this function to limit the area of detection by adding in detector size values
def Tracking_station_size(x_range,y_range,i):
    x,y,z  = Tracking_line_equations(i)
    #x_range and y_range is equal to the bounds of the tracking detector area
    while abs(x) > x_range or abs(y) > y_range:
        x,y,z  = Tracking_line_equations(i)
    
    return x,y,z

def Tracking_station_co_ords(Tracking_station_positions):
    #creating an empty array with the correct dimensions
    station_co_ords = np.empty((0,3))
    #calculating x and y co-ordinates for every tracking station position
    for i in range(len(Tracking_station_positions)):
        #if a limit was to be put on tracking station size then that function would be used instead of Tracking line equations 
        x, y, z = Tracking_line_equations(Tracking_station_positions[i])
        row = np.array([x,y,z])
        #outputting the tracking station co-ordinates as an array
        station_co_ords = np.append(station_co_ords ,[row],axis= 0)
    
    return station_co_ords

print(Tracking_station_co_ords(Tracking_station_positions))
print(Co_ord_selection())
#################################################################################################################################################################################################
                                                            # PART D HIT SMEARING
#################################################################################################################################################################################################

#Smear function that uses the x,y co-ordinates as the centre mean and station resolution as the standard deviation to output a smeared point
def Hit_smearing(i):
    
    smeared_point = np.random.normal(i, station_resolution)
    
    return smeared_point
     
def Smear_list():
    # creating an empty smear list with all the right dimensions
    smear_list = np.empty((0,3))
    #for every tracking station there will be a different smear on the x,y co ordinates. This function appends each stations new co-ordinates together
    for i in range(len(Tracking_station_positions)):
        x ,y, z = Tracking_station_co_ords(Tracking_station_positions)[i,0], Tracking_station_co_ords(Tracking_station_positions)[i,1], Tracking_station_co_ords(Tracking_station_positions)[i,2]
        row = np.array([Hit_smearing(x),Hit_smearing(y),z])
        smear_list = np.append(smear_list ,[row],axis= 0)

    return smear_list
#################################################################################################################################################################################################
                                                            # PART E TRACK RECONSTRUCTION
#################################################################################################################################################################################################

def Track_Matrix_reconstruction():
    #taking the smeared co-ordinate values for all the tracking stations
    x_vals, y_vals, z_vals = Smear_list()[:,0], Smear_list()[:,1], Smear_list()[:,2]
    #creating empty j and Xhit matricies of the right size
    j = np.empty((0,4))
    Xhit = np.empty((0,1))
    
    # 𝑗𝑡𝑟𝑘 matrix of the tracking reconstruction equation. Appending the rows to an empty array per number of stations
    for i in range(len(z_vals)):
        j = np.append(j, [[z_vals[i],0,1,0]], axis = 0)
        j = np.append(j, [[0,z_vals[i],0,1]], axis = 0)
    # 𝑥ℎ𝑖𝑡𝑠 matrix of the tracking reconstruction equation.  Appending the rows to an empty array per number of stations
    for i in range(len(x_vals)):
        Xhit = np.append(Xhit,[[x_vals[i]]], axis = 0)
        Xhit = np.append(Xhit,[[y_vals[i]]], axis = 0)
        
    #print(j,Xhit)
    #Using the (Moore-Penrose) pseudo inverse matrix for a non square matrix
    M = np.matmul(sci.linalg.pinv(j),Xhit)
    #Using the least squares method
    #M = sci.linalg.lstsq(j,Xhit)
    #Calculating the variables for the striaght line equation using the inverse Matrix
    mx, my, cx, cy = M[0],M[1],M[2],M[3]
    #print(mx,my,cx,cy)
    #print(M)
    
    return mx, my, cx, cy

start_matrix = timeit.timeit()
Track_Matrix_reconstruction()
end_matrix = timeit.timeit()
print(end_matrix-start_matrix)
#################################################################################################################################################################################################
                                                            # PART F VERTEX RECONSTRUCTION
#################################################################################################################################################################################################

#This function calcualtes the closest point to the decay vertex that it can find using the co -efficients found from the
#least squares approximation in the function above
def Vertex_reconstruction():
    
    mx,my,cx,cy = Track_Matrix_reconstruction()
    
    Z = (mx*cx + my*cy)/(mx**2 + my**2)
    
    return abs(Z)/10
#calculating the lifetime of the mother particle
def Lifetime(decay_vertex):
    
    lifetime = decay_vertex / mean_X_velocity
    
    return lifetime

#This function simulates multiple particles hit on the detector
def Finding_point():
    #calculating the vertex from multiple particles
    vertex_list =np.array([ Vertex_reconstruction() for i in range(100)])

    #creating a histogram of all the decay verticies
    n, bins, patches = plt.hist(vertex_list, bins = 50)
    plt.title('Number of verticies')
    plt.xlabel('decay vertex')
    plt.ylabel('Number of values')
    plt.plot()
    plt.show()
    # finding the most occupied bin from the distribution (which should be normal)
    mode_index = n.argmax()
    most_populated_decay_vertex = str(bins[mode_index])
    print('the most frequent bin:' , most_populated_decay_vertex )
    #calcualting the lifetime for decay vertex at this bin
    lifetime = str(Lifetime(bins[mode_index]))
    print('this is the particle lifetime' + str(Lifetime(bins[mode_index])))
    #calling the decay_vertex variable
    
    return

#print(Vertex_reconstruction())
start_finding = timeit.timeit()
Finding_point()
end_finding= timeit.timeit()
print(end_finding - start_finding)
################################################################################################################################################################################################
    