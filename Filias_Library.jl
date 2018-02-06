function squareroot(x)
	z = x # Initial starting point for Newton’s method
	while abs(z*z - x) > 1e-13
	z = z - (z*z-x)/(2z)
	end
	return z
end
JuMP.register(m, :squareroot, 1, squareroot, autodiff=true)
