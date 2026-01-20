# IFFL sequestration based model

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions
    using OrderedCollections  
      
    # ODE system with feedback
    function IFFL(du, u, p, t)
        U, Up, X = u
        k1, k2, mU, mX, d, f, A = p
    
        if f != 0
            Y = A*(sin(2*pi*t*f) + 1)
        else
            Y = A
        end
    
        du[1] = mU - d*U  -  k1*Y*U + k2*X*Up
        du[2] = - d*Up + k1*Y*U - k2*X*Up
        du[3] = mX*Y - d*X
    end

    # Parameters 
    p = OrderedDict(
        :k1 => 1, 
        :k2 => 100,
        :mU => 1.0,
        :mX => 1.0,
        :d => 1.0,
        :f => 0.5,
        :A => 5.0
    )
 

    u0 = [1.0, 1.0, 1.0]     # Initial conditions for the ODE system
end