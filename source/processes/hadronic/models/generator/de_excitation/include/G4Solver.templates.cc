template <class Function>
G4bool G4Solver<Function>::Bisection(Function & theFunction)
{
	// Check the interval before start
	if (a > b || abs(a-b) <= tolerance) {
		G4cerr << "G4Solver::Bisection: The interval must be properly set." << G4endl;
		return false;
	}
	G4double fa = theFunction(a);
	G4double fb = theFunction(b);
	if (fa*fb > 0.0) {
		G4cerr << "G4Solver::Bisection: The interval must include a root." << G4endl;
		return false;
	}
	
	G4double eps=tolerance*(b-a);
	
	
	// Finding the root
	for (G4int i = 0; i < MaxIter; i++) {
		G4double c = (a+b)/2.0;
		if ((b-a) < eps) {
			root = c;
			return true;
		}
		G4double fc = theFunction(c);
		if (fc == 0.0) {
			root = c;
			return true;
		}
		if (fa*fc < 0.0) {
			a=c;
			fa=fc;
		} else {
			b=c;
			fb=fc;
		}
	}
	G4cerr << "G4Solver::Bisection: Excedded maximum number of iterations whithout convegence." << G4endl;
	return false;
}
	

template <class Function>
G4bool G4Solver<Function>::RegulaFalsi(Function & theFunction)
{
	// Check the interval before start
	if (a > b || abs(a-b) <= tolerance) {
		G4cerr << "G4Solver::RegulaFalsi: The interval must be properly set." << G4endl;
		return false;
	}
	G4double fa = theFunction(a);
	G4double fb = theFunction(b);
	if (fa*fb > 0.0) {
		G4cerr << "G4Solver::RegulaFalsi: The interval must include a root." << G4endl;
		return false;
	}
	
	G4double eps=tolerance*(b-a);
	
	
	// Finding the root
	for (G4int i = 0; i < MaxIter; i++) {
		G4double c = (a*fb-b*fa)/(fb-fa);
		G4double delta = G4std::min(abs(c-a),abs(b-c));
		if (delta < eps) {
			root = c;
			return true;
		}
		G4double fc = theFunction(c);
		if (fc == 0.0) {
			root = c;
			return true;
		}
		if (fa*fc < 0.0) {
			b=c;
			fb=fc;
		} else {
			a=c;
			fa=fc;
		}
	}
	G4cerr << "G4Solver::Bisection: Excedded maximum number of iterations whithout convegence." << G4endl;
	return false;

}	

template <class Function>
G4bool G4Solver<Function>::Brent(Function & theFunction)
{
	
	const G4double precision = 3.0e-8;
	
	// Check the interval before start
	if (a > b || abs(a-b) <= tolerance) {
		G4cerr << "G4Solver::Brent: The interval must be properly set." << G4endl;
		return false;
	}
	G4double fa = theFunction(a);
	G4double fb = theFunction(b);
	if (fa*fb > 0.0) {
		G4cerr << "G4Solver::Brent: The interval must include a root." << G4endl;
		return false;
	}
	
	G4double c = b;
	G4double fc = fb;
	G4double d = 0.0;
	G4double fd = 0.0;
	G4double e = 0.0;
	G4double fe = 0.0;
	
	
	
	for (G4int i=0; i < MaxIter; i++) {
		// Rename a,b,c and adjust bounding interval d
		if (fb*fc > 0.0) {
			c = a;
			fc = fa;
			d = b - a;
			e = d;
		}
		if (abs(fc) < abs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		G4double Tol1 = 2.0*precision*abs(b) + 0.5*tolerance;
		G4double xm = 0.5*(c-b);
		if (abs(xm) <= Tol1 || fb == 0.0) {
			root = b;
			return true;
		}
		// Inverse quadratic interpolation
		if (abs(e) >= Tol1 && abs(fa) > abs(fb)) {
			G4double s = fb/fa;
			G4double p = 0.0;
			G4double q = 0.0;
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0 - s;
			} else {
				q = fa/fc;
				G4double r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			// Check bounds
			if (p > 0.0) q = -q;
			p = abs(p);
			G4double min1 = 3.0*xm*q-abs(Tol1*q);
			G4double min2 = abs(e*q);
			if (2.0*p < G4std::min(min1,min2)) {
				// Interpolation
				e = d;
				d = p/q;
			} else {
				// Bisection
				d = xm;
				e = d;
			}
		} else {
			// Bounds decreasing too slowly, use bisection
			d = xm;
			e = d;
		}
		// Move last guess to a 
		a = b;
		fa = fb;
		if (abs(d) > Tol1) b += d;
		else {
			if (xm >= 0.0) b += abs(Tol1);
			else b -= abs(Tol1);
		}
		fb = theFunction(b);
	}
	G4cerr << "G4Solver::Brent: Number of iterations exceeded." << G4endl;
	return false;
}

template <class Function> 	
void G4Solver<Function>::SetIntervalLimits(const G4double Limit1, const G4double Limit2)
{
	if (abs(Limit1-Limit2) <= tolerance) {
		G4cerr << "G4Solver::SetIntervalLimits: Interval must be wider than tolerance." << G4endl;
		return;
	}
	if (Limit1 < Limit2) {
		a = Limit1;
		b = Limit2;
	} else {
		a = Limit2;
		b = Limit1;
	}
	return;
}