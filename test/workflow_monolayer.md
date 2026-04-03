Workflow:

Use the 1D model to identify the repulsion force, k, required such that 
90%  relaxation occurs at t=1.

For us, this turns out to be k = 88.54.

Now we know the repulsion force, and the time, we slow down the area 
growth rate, growth_rate.

This gives growth_rate = 0.2

Now we know the area growth rate, and the force, we can proceed.

The way we do this is with the following algorithm

for (time = 0 : dt: final_time)

    Compute f_i, a_i
    cell_inhibition_i = Is_Cell_Inhibited()

    if(~cell_inhibition_i)
        Grow_And_Divide_Cell()
    
    Update_Cell_Positions()

def Grow_And_Divide_Cell()
    Phase_Time








