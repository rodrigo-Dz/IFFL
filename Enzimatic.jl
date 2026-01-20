# IFFL sequestration based model

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions
    using OrderedCollections  
      
    # ODE system with feedback
    function IFFL(du, u, p, t)
        U1, U2, X, C = u
        m1, m2, mx, g, d, r, f, A = p
    
        if f != 0
            Y = A*(sin(2*pi*t*f) + 1)
        else
            Y = A
        end
    
        du[1] = m1*Y - d*U1 - g*U1*U2 
        du[2] = m2*X - d*U2 - g*U1*U2 + r*C
        du[3] = mx*Y - d*X
        du[4] = g*U1*U2 - r*C - d*C
    end

    # Parameters 
    p = OrderedDict(
    :m1 => 1.0,
    :m2 => 1.0,
    :mx => 1.0, 
    :g => 100.0, 
    :d => 1.0,
    :r => 100,
    :f => 0.5,
    :A => 5.0
    )
 

    u0 = [0.0, 0.0, 0.0, 0.0]     # Initial conditions for the ODE system
end