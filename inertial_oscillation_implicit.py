# Numerical solution of the inertial oscillation equation using
# forward-backward time-stepping with Coriolis parameter f
# dudt = f*v
# dvdt = -f*u

import numpy as np    # External library for numerical calculations
import matplotlib.pyplot as plt    # Plotting library

# Put everything inside a main function to avoid global variables
def main():
    # setup parameters
    f = 1e-4                # Coriolis parameter
    nt = 100               # Number of time steps
    dt = 5000             # Time step in seconds
    
    # Initial conditions (in meters or m/s)
    x0 = 0.
    y0 = 0
    u0 = 10.
    v0 = 0.
    # Initialise velocity from initial conditions

    #u = u0
    #v = v0
    # Store all locations for plotting and store initial locations
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    u = np.zeros(nt+1)
    v = np.zeros(nt+1)
    x[0] = x0
    y[0] = y0
    u[0] = u0
    v[0] = v0
    # Loop over all time-steps
    for n in range(nt):
        u[n+1] = (u[n] + f*v[n]*dt)/(1 + (f*dt)**2)
        v[n+1] = (v[n] - f*u[n]*dt)/(1 + (f*dt)**2)
        x[n+1] = x[n] + dt*u[n+1]
        y[n+1] = y[n] + dt*v[n+1]
       
    # Analytic solution for the location as a function of time
    times = np.linspace(0,nt*dt, nt+1)
    xa = x0 + 1/f*(u0*np.sin(f*times) - v0*np.cos(f*times) + v0)
    ya = y0 + 1/f*(u0*np.cos(f*times) + v0*np.sin(f*times) - u0)
    ua = u0*np.cos(f*times) + v0*np.sin(f*times)
    va = -u0*np.sin(f*times) + v0*np.cos(f*times)
    
    # Plot the solution in comparison to the analytic solution
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.plot(xa, ya, '-k+', label='analytical')
    plt.plot(x, y, '-bo', label='implicit')
    plt.legend(loc='best',fancybox=True)
    plt.xlabel('x',fontsize=15)
    plt.ylabel('y',fontsize=15)
    #plt.axhline(0, linestyle=':', color='black')
    #ax1.show()
    
    plt.subplot(2,1,2)
    plt.plot(ua, va, '-k+', label='analytical')
    plt.plot(u, v, '-bo', label='implicit')
    plt.legend(loc='best')
    plt.xlabel('u',fontsize=15)
    plt.ylabel('v',fontsize=15)
    #plt.axhline(0, linestyle=':', color='black')


    # plot solutions of x, y, u and v compared to the analytical solution
    plt.figure(2,figsize=(12,8))
    plt.subplot(2,2,1)
    plt.plot(xa,'-k+',label='analytical')
    plt.plot(x,'-bo',label='forward-backward')
    plt.ylabel('x',fontsize=15)

    plt.subplot(2,2,2)
    plt.plot(ya,'-k+',label='analytical')
    plt.plot(y,'-bo',label='forward-backward')
    plt.ylabel('y',fontsize=15)
    #plt.legend(loc='best')
    
    plt.subplot(2,2,3)
    plt.plot(ua,'-k+',label='analytical')
    plt.plot(u,'-bo',label='forward-backward')
    plt.ylabel('u',fontsize=15)
    plt.xlabel('Number of time steps',fontsize=15)
    plt.subplot(2,2,4)
    plt.plot(va,'-k+',label='analytical')
    plt.plot(v,'-bo',label='forward-backward')
    plt.ylabel('v',fontsize=15)
    plt.xlabel('Number of time steps',fontsize=15)
    plt.show()

# Execute the code
main()


