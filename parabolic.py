import numpy as np 
import matplotlib.pyplot as plt

print("Enter initial velocity in m/s")
v = float(input())
g = 9.8 #Gravity acceleration in m/s2
theta_d = np.arange(0,100,10) #List of the angle
theta_r = np.radians(theta_d) #Convert degrees to radians
t = np.arange(0,10,0.1) #Define the time interval setting as a continous parameter
#Calculate trajectory for every angle
for u in theta_r: #Trajectory for each angle
    x = [] #Empty list of x and y coordinates
    y = []
    for k in t: #Create a loop for calculating the coordinates
        x.append(v*k*np.cos(u))
        y.append(v*k*np.sin(u)-0.5*g*k*k)
    #Plot the results
    plt.plot(x, y)
    plt.xlabel('Distance (m)')
    plt.ylabel('Height (m)')
    plt.title('Projectile motion')
    plt.legend(theta_d)
plt.show()
