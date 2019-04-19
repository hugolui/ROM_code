import numpy as np
from func import func
 
def RK45(x0,parameters,h,n_var,nt):
    
    X_s = np.zeros((n_var,nt), dtype = "float32") # Solution
    X_s[:,0] = x0 # Initial solution
    
    # Runge-Kutta coefficients
    A2 = -6234157559845/12983515589748.
    A3 = -6194124222391/4410992767914.
    A4 = -31623096876824/15682348800105.
    A5 = -12251185447671/11596622555746.
    B1 = 494393426753/4806282396855.
    B2 = 4047970641027/5463924506627.
    B3 = 9795748752853/13190207949281.
    B4 = 4009051133189/8539092990294.
    B5 = 1348533437543/7166442652324.
    
    for n in range(nt-1):  

        # Step 1
        dx_s1 = h*func(X_s[:,[n]],parameters)
        x_s1 =  X_s[:,[n]] + B1*dx_s1
        
        # Step 2
        dx_s2 = A2*dx_s1 + h*func(x_s1,parameters)
        x_s2 = x_s1 + B2*dx_s2
        
        # Step 3
        dx_s3 = A3*dx_s2 + h*func(x_s2,parameters)
        x_s3 = x_s2 + B3*dx_s3
        
        # Step 4
        dx_s4 = A4*dx_s3 + h*func(x_s3,parameters)
        x_s4 = x_s3 + B4*dx_s4
        
        # Step 5
        dx_s5 = A5*dx_s4 + h*func(x_s4,parameters)
        X_s[:,n+1] = (x_s4 + B5*dx_s5).T
    
    
    return X_s