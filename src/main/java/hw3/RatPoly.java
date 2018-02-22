package hw3;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

/**
 * <b>RatPoly</b> represents an immutable single-variate polynomial expression.
 * RatPolys are sums of terms with rational coefficients and non-negative exponents.
 * <p>
 *
 * Examples of RatPolys include "0", "x-10", and "x^3-2*x^2+5/3*x+3", and "NaN".
 */
// See RatNum's documentation for a definition of "immutable".
public final class RatPoly {

    /** Holds all the RatNum coefficients in this RatPoly */
    private final RatNum[] coeffs;
    /** Holds the degree of this RatPoly */
    private final int degree;

    // Abstraction Function:
    // RatPoly, p, represents the polynomial equal to the sum of the terms:
    // sum (0 <= i < length(p)): p.coeffs[i]*x^i
    // If there are no coefficients, p represents the "0" polynomial.
    //
    // Representation Invariant for every RatPoly p:
    // coeffs != null &&
    // foreach i, 0<=i<coeffs.length: coeffs[i] != null &&
    // if p is the "0" polynomial
    //      coeffs is an empty array and degree = 0
    // else
    //      degree = coeffs.length - 1 &&
    //      coeffs[degree] != 0
    //
    // In other words:
    // * The coeffs field always points to some usable object array.
    // * No coefficient is null.
    // * The degree field is the highest power and should be one less 
    // * than the size of the coeffs array. 
    // * The coefficient of the highest-power term must be non-zero.
    
    
    /** A constant holding a Not-a-Number (NaN) value of type RatPoly */
    public static final RatPoly NaN = new RatPoly(new RatNum[] { RatNum.NaN });

    /** A constant holding a zero value of type RatPoly */
    public static final RatPoly ZERO = new RatPoly();

    /**
     * @effects Constructs a new Poly with value "0".
     */
    public RatPoly() {
        coeffs = new RatNum[0];
        degree = 0;        
        checkRep();
    }

    
    /**
     * @param c The constant in the term which the new RatPoly equals.
     * @param e The exponent in the term which the new RatPoly equals.
     * @requires e >= 0
     * @effects Constructs a new RatPoly equal to "c*x^e". If c is zero, constructs
     *          a "0" polynomial.
     */
    public RatPoly(int c, int e) {
        // TODO: Fill in this method, then remove the RuntimeException
    	if (c == 0) // constructs a "0" polynomial as indicated
    	{
    		coeffs = new RatNum[0];
    		degree = 0;
    	}
    	else // Otherwise creates an array of size e+1 full of zeroes except for c (as a RatNum) at the end
    	{
    		coeffs = new RatNum[e+1];
    		for(int i = 0; i < e+1; i++)
    		{
    			if (i == e)
    			{
    				RatNum coe = new RatNum(c);
    				coeffs[i] = coe;
    			}
    			else
    			{
    				RatNum coe = RatNum.ZERO;
    				coeffs[i] = coe;
    			}
    		}
    		degree = e; // Sets degree equal to the highest degree, or e.
    	}
    	checkRep();  // Ensures structure representation is being followed.
    }

    
    /**
     * @param coeffs An array of coefficients to be contained in the new RatPoly.
     * @requires 'coeffs' is non-empty and it satisfies clauses given in rep. invariant
     * @effects Constructs a new Poly using 'coeffs' as part of the representation.
     *          The method does not make a copy of 'coeffs'.
     */
    public RatPoly(RatNum[] coeffs) {
        this.coeffs = coeffs;
        this.degree = coeffs.length - 1;
        // The spec tells us that we don't need to make a copy of 'coeffs' 
        // (argument satisfies the clauses of the rep. invariant
        checkRep();
    }

    /**
     * Returns the degree of this RatPoly.
     *
     * @requires !this.isNaN()
     * @return the largest exponent with a non-zero coefficient, or 0 if this is
     *         "0".
     */
    public int degree() {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep(); // Checks representation on entry and makes sure it is not null, if so returns 0
    	if (!this.isNaN())
    	{
    		int coe = coeffs.length; // Returns length-1 or 0, whichever is bigger.
    		if (coe == 0)
    		{
    			checkRep();
    			return 0;
    		}
    		return coeffs.length-1; 
    	}
    	else
    	{
    		checkRep();
    		return 0;
    	}
    }

    /**
     * Gets the coefficient of the term of power 'pow'
     *
     * @param pow The power for which to find the corresponding coefficient.
     * @requires !this.isNaN()
     * @return the RatNum that is the coefficient of the term of power 'pow'.
     *         "0" if this is "0" || pow < 0 || pow >= coeffs.size 
     */
    public RatNum getCoeff(int pow) {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
    	if (!this.isNaN()) // checks if pow is a valid power then returns the coefficient in that location
    	{
    		if(pow < 0 || pow >= coeffs.length)
    		{
    			checkRep();
    			return RatNum.ZERO;
    		}
    		checkRep();
    		return coeffs[pow];
    	}
    	else
    	{
    		checkRep();
    		return RatNum.ZERO;
    	}
    }

    /**
     * Returns true if this RatPoly is not-a-number.
     *
     * @return true if and only if this has some coefficient = "NaN".
     */
    public boolean isNaN() {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
    	for(int i = 0; i < coeffs.length; i++) // checks if any individual element isNaN using the RaTNum
    	{ // check available, if so returns true.
    		if(coeffs[i].isNaN())
    		{
    			checkRep();
    			return coeffs[i].isNaN();
    		}
    	}
    	checkRep();
    	return false;
       // throw new RuntimeException("RatPoly.isNaN() is not yet implemented");    	
    }

        
    /**
     * Scales coefficients within 'arr' by 'scalar' (helper procedure).
     *
     * @param arr The RatNums to be scaled.
     * @param scalar the value by which to scale coefficients in arr.
     * @requires arr, scalar != null
     * @modifies arr
     * @effects Forall i s.t. 0 <= i < arr.length, arr_post[i] = arr_pre[i]*scalar
     *          
     */
    private static void scaleCoeff(RatNum[] arr, RatNum scalar) {
        // TODO: Fill in this method, then remove the RuntimeException
    	for(int i = 0; i < arr.length; i++) // Changes the arrary to have every non NaN element multiplied
    	{ // by the scalar. Because it is a helper function no need to check Rep
    		if(!(arr[i].isNaN()))
    		{
    			arr[i] = scalar.mul(arr[i]);
    		}
    	}
    }


    /**
     * Return the additive inverse of this RatPoly.
     *
     * @return a RatPoly equal to "0 - this"; if this.isNaN(), returns some r
     *         such that r.isNaN()
     */
    public RatPoly negate() {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
    	if(this.isNaN()) // Returns this (this is NaN) if it is NaN, zero if it is zero
    	{
    		return this;
    	}
    	if(this.equals(ZERO))
    	{
    		return ZERO;
    	}
    	else // Otherwise it takes a clone (new copy) then finds the negation with sub, checks rep, returns
    	{
			RatNum zero = new RatNum(0);
			RatNum[] newcoe = this.coeffs.clone();
			RatPoly coeffspost = new RatPoly(newcoe);
    		for(int i = 0; i < coeffs.length; i++)
    		{
    			coeffspost.coeffs[i] = zero.sub(coeffspost.coeffs[i]);
    		}
    		coeffspost.checkRep();
    		return coeffspost;
    	}
    }

    /**
     * Addition operation.
     *
     * @param p The other value to be added.
     * @requires p != null
     * @return a RatPoly, r, such that r = "this + p"; if this.isNaN() or
     *         p.isNaN(), returns some r such that r.isNaN()
     */
    public RatPoly add(RatPoly p) {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
    	if(this.isNaN()) // If either is NaN, returns that one, if either is zero, returns the other
    	{
    		return this;
    	}
    	if(p.isNaN())
    	{
    		return p;
    	}
    	if(this.equals(ZERO))
    	{
    		return p;
    	}
    	if(p.equals(ZERO))
    	{
    		return this;
    	}
    	else // Takes whichever polynomial has a higher degree and adds the other one to it, so
    	{ // there is no need to add more to the end. Then takes a new copy of the longer one and 
    		if(p.degree() > this.degree()) // adds everything to it, if things get reduced to zero it checks
    		{ // and changes the return as needed.
    			int num_of_z = 0;
    			RatNum[] newcoe = p.coeffs.clone();
    			RatPoly r = new RatPoly(newcoe);
    			for(int i = 0; i < this.degree()+1; i++)
    			{
    				r.coeffs[i] = coeffs[i].add(p.getCoeff(i));
    				if (r.coeffs[i].equals(RatNum.ZERO))
    				{
    					num_of_z+= 1;
    				}
    			}
    			if(num_of_z == p.degree()+1)
    			{
    				r = new RatPoly();
    				return r;
    			}
    	    	r.checkRep();
    			return r;
    		}
    		else // The same but using the other polynomial as a base
    		{
    			int num_of_z = 0;
    			RatNum[] newcoe = this.coeffs.clone();
    			RatPoly r = new RatPoly(newcoe);
    			for(int i = 0; i < p.degree()+1; i++)
    			{
    				r.coeffs[i] = coeffs[i].add(p.getCoeff(i));
    				if (r.coeffs[i].equals(RatNum.ZERO))
    				{
    					num_of_z+= 1;
    				}
    			}
    			if(num_of_z == p.degree()+1)
    			{
    				r = new RatPoly();
    				return r;
    			}
    	    	r.checkRep();
    			return r;
    		}
    	}
    }
    
    /**
     * Subtraction operation.
     *
     * @param p The value to be subtracted.
     * @requires p != null
     * @return a RatPoly, r, such that r = "this - p"; if this.isNaN() or
     *         p.isNaN(), returns some r such that r.isNaN()
     */
    public RatPoly sub(RatPoly p) { // The same as add but subtracting instead of adding, RatNum's .sub
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
		RatNum zero = new RatNum(0); // The only real big difference is the active removal of trailing 
		RatPoly zeroo = new RatPoly(); // zeros necessary to ensure that the div. operation works. 
    	if(this.isNaN()) 
    	{
    		return this;
    	}
    	if(p.isNaN())
    	{
    		return p;
    	}
    	if(this.equals(zeroo)) // Returns the negation if the other side if zero, vica versa
    	{
    		return p.negate();
    	}
    	if(p.equals(zeroo))
    	{
    		return this;
    	}
    	if(this.equals(p)) // Returns the 0 polynomial if this and p are the same polynomial
    	{
    		return new RatPoly();
    	}
    	else
    	{
    		if(p.degree() > this.degree())
    		{
    			int num_of_z = 0;
    			RatNum[] newcoe = p.coeffs.clone();
    			RatPoly r = new RatPoly(newcoe);
    			for(int i = 0; i < this.degree()+1; i++)
    			{
    				r.coeffs[i] = coeffs[i].sub(p.getCoeff(i));
    				if (r.coeffs[i].equals(zero))
    				{
    					num_of_z+= 1;
    				}
    			}
    			if(num_of_z == p.degree()+1)
    			{
    				r = new RatPoly();
    				r.checkRep();
    				return r;
    			}
    			
    			int trailing_z = 0;
    			for(int i = r.coeffs.length-1; i >= 0; i--)
    			{
    				if(r.coeffs[i].equals(RatNum.ZERO))
    				{
    					trailing_z += 1;
    				}
    				else
    				{
    					break;
    				}
    			}
    			if(trailing_z != 0)
    			{
    				RatNum[] fixed = new RatNum[r.coeffs.length-trailing_z];
    				for(int i = 0; i < p.degree()+1-trailing_z;i++)
        			{
        				fixed[i] = r.coeffs[i];
        			}
    				r = new RatPoly(fixed);
    			}
    			r.checkRep();
    			return r;
    		}
    		else
    		{
    			int num_of_z = 0;
    			RatNum[] newcoe = this.coeffs.clone();
    			RatPoly r = new RatPoly(newcoe);
    			for(int i = 0; i < p.degree()+1; i++)
    			{
    				if(i < this.degree() + 1)
    				{
    					r.coeffs[i] = coeffs[i].sub(p.getCoeff(i));
    					if (r.coeffs[i].equals(zero))
        				{
        					num_of_z+= 1;
        				}
    				}
    			}
    			if(num_of_z == p.degree()+1)
    			{
    				r = new RatPoly();
    				r.checkRep();
    				return r;
    			}
    			int trailing_z = 0; // Checks for and removes trailing zeros
    			for(int i = r.coeffs.length-1; i >= 0; i--)
    			{
    				if(r.coeffs[i].equals(RatNum.ZERO))
    				{
    					trailing_z += 1;
    				}
    				else
    				{
    					break;
    				}
    			}
    			if(trailing_z != 0)
    			{
    				RatNum[] fixed = new RatNum[r.coeffs.length-trailing_z];
    				for(int i = 0; i < p.degree()+1-trailing_z;i++)
        			{
        				fixed[i] = r.coeffs[i];
        			}
    				r = new RatPoly(fixed);
    			}
    			r.checkRep();
    			return r;
    		}
    	}
    }

    /**
     * Multiplication operation.
     *
     * @param p The other value to be multiplied.
     * @requires p != null
     * @return a RatPoly, r, such that r = "this * p"; if this.isNaN() or
     *         p.isNaN(), returns some r such that r.isNaN()
     */
    public RatPoly mul(RatPoly p) {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
    	if(this.isNaN()) // Checks rep, then if anything is NaN, then if anything is zero and returns that
    	{
    		return this;
    	}
    	if(p.isNaN())
    	{
    		return p;
    	}
    	if(this.equals(ZERO) || p.equals(ZERO))
    	{
    		RatPoly r = new RatPoly();
    		r.checkRep();
    		return r;
    	}
		RatNum[] newcoeffs = new RatNum[this.degree()+p.degree()+1]; // Otherwise it creates a RatNum array
		for(int i = 0; i < newcoeffs.length; i++) // of appropriate length, also a temp array clone of p.
		{
			newcoeffs[i] = RatNum.ZERO;
		}
		RatNum[] temp;
    	for(int i = 0; i < this.coeffs.length; i++) // uses the temp array to be scalar'd by every element
    	{ // in this.coeff, then adding the array to the base sum, getting the multiplication done.
			temp = p.coeffs.clone();
    		RatPoly.scaleCoeff(temp, this.coeffs[i]);
    		for(int j = 0; j < temp.length; j++)
    		{
    			newcoeffs[j + i] = newcoeffs[j + i].add(temp[j]); 
    		}
    	}
    	RatPoly r = new RatPoly(newcoeffs); // After its been summed up, creates a polynomial and checks rep
    	r.checkRep();
    	return r;
    }

	public RatPoly div(RatPoly p) // Division function
	{
		// TODO Auto-generated method stub
		checkRep();
		if(p.equals(ZERO)) // Numerous checks involving NaNs or zeros
		{
			return RatPoly.NaN;
		}
		if(this.isNaN() || p.isNaN())
		{
			return RatPoly.NaN;
		}
		if(this.equals(ZERO))
		{
			return RatPoly.ZERO;
		}
		int deg = this.degree()-p.degree(); // Finds the initial degree, if one is initially bigger return 0
		if(deg < 0) // if it is the denominator.
		{
			return RatPoly.ZERO;
		}
		RatNum[] q = new RatNum[deg+1]; // creates an array for the quotient of a specific size based on deg
		for(int i = 0; i < q.length; i++) // fills it with zeros
		{
			q[i] = RatNum.ZERO;
		}
		RatNum[] temp_pc= p.coeffs.clone(); // I create a bunch of helper Polynomials and arrays
		RatPoly temp_p = new RatPoly(temp_pc);
		RatNum[] temp_thisc = this.coeffs.clone();
		RatPoly temp_this = new RatPoly(temp_thisc);
		RatNum[] temp;
		for(int i = 0; i < this.degree()+1; i++) //iterates through every degree of object being divided,
		{ // which would be max number of times
			temp_pc = p.coeffs.clone();
			temp_p = new RatPoly(temp_pc);
			if(deg < 0) // if deg is less than 0, it halts
			{
				break;
			}
			if(deg == 0) // if the degree is zero is adds the last element in then halts
			{
				q[0] = temp_this.coeffs[temp_this.degree()].div(temp_p.coeffs[temp_p.degree()]);
				deg = -2; // redundant
				break;
			}  // finds the coefficient of the highest powers
			q[deg] = temp_this.coeffs[temp_this.degree()].div(temp_p.coeffs[temp_p.degree()]);
			RatPoly.scaleCoeff(temp_pc, q[deg]); // scales the coefficients to be the same
			temp = new RatNum[temp_pc.length+deg]; // creates a new array and scales up p so it can be
			for(int j = 0; j < temp.length; j++) // subtracted (essentially multiplying x^deg
			{
				temp[j] = RatNum.ZERO;
			}
			for(int j = 0; j < temp_pc.length; j++)
			{
				temp[j+deg] = temp_pc[j];
			}
			temp_p = new RatPoly(temp); 
			if(temp_this.sub(temp_p).equals(ZERO)) // if they are equal returns q as is instead of subbing
			{
				return new RatPoly(q);
			}
			else // Otherwise it subtracts, resets p to a nonscaled version, and finds the new degree,
			{ // leaving temp_this altered so the division can make progress
				temp_this = new RatPoly(temp_this.sub(temp_p).coeffs);
				temp_p = new RatPoly(p.coeffs.clone());
				deg = temp_this.degree() - temp_p.degree();
			}
		}
		checkRep(); // checks rep and returns the resulting quotient polynomial
		return new RatPoly(q);
	}
    /**
     * Returns the value of this RatPoly, evaluated at d. Evaluate using Horner's
     * rule as you did in Homework 2.
     *
     * @param d The value at which to evaluate this polynomial.
     * @return the value of this polynomial when evaluated at 'd'. For example,
     *         "x+2" evaluated at 3 is 5, and "x^2-x+1" evaluated at 3 is 7. if
     *         (this.isNaN() == true), return RatNum.NaN
     */
    public double eval(double d) { // Checks rep, then returns NaN if NaN or coeffs[0] if length 1.
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
    	if(this.isNaN())
    	{
    		return Double.NaN;
    	}
    	if(this.coeffs.length == 1)
    	{
    		return coeffs[0].doubleValue();
    	}
    	double sum = 0;
    	for(int i = 1; i < coeffs.length+1; i++) // Otherwise it sums using Horner's rule then returns sum.
    	{
    		sum = sum*d;
    		sum += coeffs[coeffs.length-i].doubleValue();
    	}
    	return sum;
    }  
    
    
    /**
     * Return the derivative of this RatPoly.
     *
     * @return a RatPoly, q, such that q = dy/dx, where this == y. In other
     *         words, q is the derivative of this. If this.isNaN(), then return
     *         some q such that q.isNaN()
     *
     * <p>
     * The derivative of a polynomial is the sum of the derivative of each term.
     */
    public RatPoly differentiate() {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep();
    	if(this.isNaN()) // Checks rep, returns NaN if NaN, and returns 0 polynomial if degree and poly is 0
    	{
    		return this;
    	}
    	if(this.degree() == 0 || this == ZERO)
    	{
    		return new RatPoly();
    	}
    	else // Otherwise I create an array, alter is using the coeffs and the power rule of differentiation
    	{
    		RatNum[] newcoeffs = new RatNum[this.degree()];
        	for(int i = 1; i < this.coeffs.length; i++)
        	{
        		RatNum temp = new RatNum(i);
        		newcoeffs[i-1] = this.coeffs[i].mul(temp);
        	}
        	RatPoly q = new RatPoly(newcoeffs); // Then with the completed array I check rep and return
        	q.checkRep();
        	return q;
    	}
    }

    /**
     * Returns the antiderivative of this RatPoly.
     *
     * @param integrationConstant The constant of integration to use when
     *  computating the antiderivative.
     * @requires integrationConstant != null
     * @return a RatPoly, q, such that dq/dx = this and the constant of
     *         integration is "integrationConstant" In other words, q is the
     *         antiderivative of this. If this.isNaN() or
     *         integrationConstant.isNaN(), then return some q such that
     *         q.isNaN()
     *
     * <p>
     * The antiderivative of a polynomial is the sum of the antiderivative of
     * each term plus some constant.
     */
    public RatPoly antiDifferentiate(RatNum integrationConstant) {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep(); // The same idea as differentiation but inverted and with more base cases for an 
    	if(this.isNaN()) // empty or nearly empty polynomial, the loop is nearly the same as before.
    	{
    		return this;
    	}
    	if(this == ZERO && integrationConstant == RatNum.ZERO)
    	{
    		return new RatPoly();
    	}
    	if(this == ZERO && integrationConstant != RatNum.ZERO)
    	{
    		return new RatPoly(integrationConstant.intValue(),0);
    	}
    	else
    	{
    		RatNum[] newcoeffs = new RatNum[this.degree()+2];
			newcoeffs[0] = integrationConstant;
        	for(int i = 0; i < this.coeffs.length; i++)
        	{
        		RatNum temp = new RatNum(i+1);
        		newcoeffs[i+1] = this.coeffs[i].div(temp);
        	}
        	RatPoly q = new RatPoly(newcoeffs); // Then I still check rep and return the completed array
        	q.checkRep();
        	return q;
    	}
    }

    /**
     * Returns the integral of this RatPoly, integrated from lowerBound to
     * upperBound.
     *
     * <p>
     * The Fundamental Theorem of Calculus states that the definite integral of
     * f(x) with bounds a to b is F(b) - F(a) where dF/dx = f(x) NOTE: Remember
     * that the lowerBound can be higher than the upperBound.
     *
     * @param lowerBound The lower bound of integration.
     * @param upperBound The upper bound of integration.
     * @return a double that is the definite integral of this with bounds of
     *         integration between lowerBound and upperBound. If this.isNaN(),
     *         or either lowerBound or upperBound is Double.NaN, return
     *         Double.NaN.
     */
    public double integrate(double lowerBound, double upperBound) {
        // TODO: Fill in this method, then remove the RuntimeException
    	checkRep(); // I check rep then return NaN if NaN, and call differentiate and eval functions to 
    	if(this.isNaN() || Double.isNaN(lowerBound) || Double.isNaN(upperBound)) // determine the integral. 
    	{
    		return Double.NaN;
    	}
    	RatNum temp = new RatNum(0);
    	RatPoly i = this.antiDifferentiate(temp);
    	checkRep(); // Check rep just in case.
    	return i.eval(upperBound)-i.eval(lowerBound);
    }

    
    /**
     * Returns a string representation of this RatPoly.
     *
     * @return A String representation of the expression represented by this,
     *         with the terms sorted in order of degree from highest to lowest.
     *         <p>
     *         There is no whitespace in the returned string.
     *         <p>
     *         If the polynomial is itself zero, the returned string will just
     *         be "0".
     *         <p>
     *         If this.isNaN(), then the returned string will be just "NaN"
     *         <p>
     *         The string for a non-zero, non-NaN poly is in the form
     *         "(-)T(+|-)T(+|-)...", where "(-)" refers to a possible minus
     *         sign, if needed, and "(+|-)" refer to either a plus or minus
     *         sign, as needed. For each term, T takes the form "C*x^E" or "C*x"
     *         where C > 0, UNLESS: (1) the exponent E is zero, in which case T
     *         takes the form "C", or (2) the coefficient C is one, in which
     *         case T takes the form "x^E" or "x". In cases were both (1) and
     *         (2) apply, (1) is used.
     *         <p>
     *         Valid example outputs include "x^17-3/2*x^2+1", "-x+1", "-1/2",
     *         and "0".
     *         <p>
     */
    @Override
    public String toString() {
        if (coeffs.length == 0) {
        	return "0";
        }
    	if (isNaN()) {
            return "NaN";
        }
        
        StringBuilder output = new StringBuilder();
        boolean isFirst = true;
        for (int i=coeffs.length-1; i>=0; i--) {
        	// We print nothing for terms with 0 coefficients
        	if (coeffs[i].compareTo(RatNum.ZERO) == 0) continue;
            String term = formatTerm(coeffs[i],i);
        	if (isFirst) {
                isFirst = false;
                output.append(term);
            } 
            else {
            	if (coeffs[i].isNegative()) {
            		output.append(term);
                } else {
                    output.append("+" + term);
                }
            }
        }
        
        return output.toString();
    }
    /**
     * Helper function. Formats the term with coefficient c and exponent 
     * e according to the toString rules.
     * 
     * @param c is the coefficient of the term, e is the exponent
     * @requires c non-null and e >= 0
     * @returns a new string representing the formatted term
     * 
     */
    private static String formatTerm(RatNum c, int e) {
    	StringBuilder output = new StringBuilder();
    	if (e == 0) {
    		// if exponent is 0, add the string representation of coefficient
    		output.append(c.toString());
    	}
    	else if (c.compareTo(new RatNum(1)) == 0) {
    		// if e != 0 and coefficient is 1, skip coefficient
    		output.append("x");
    	}
    	else if (c.compareTo(new RatNum(-1)) == 0) {
    		// if e != 0 and coefficient is -1, skip coefficient
    		output.append("-x");
    	}
    	else {
    		// if e != 0 and |coefficient| != 1, add as expected
    		output.append(c.toString());
    		output.append("*x");
    	}
    	if (e>1) output.append("^"+e);
    	
    	return output.toString();
    }
    
    /**
     * Builds a new RatPoly, given a descriptive String.
     *
     * @param polyStr A string of the format described in the @requires clause.
     * @requires 'polyStr' is an instance of a string with no spaces that
     *           expresses a poly in the form defined in the toString() method.
     *           <p>
     *
     * Valid inputs include "0", "x-10", and "x^3-2*x^2+5/3*x+3", and "NaN".
     *
     * @return a RatPoly p such that p.toString() = polyStr
     */
    public static RatPoly valueOf(String polyStr) {

    	final class RatTerm {
    	   RatNum coefficient;
    	   int exponent;
     	}
    	if (polyStr.equals("0")) return RatPoly.ZERO;
    	if (polyStr.equals("NaN")) return RatPoly.NaN;
    	
        List<RatTerm> parsedTerms = new ArrayList<RatTerm>();

        // First we decompose the polyStr into its component terms;
        // third arg orders "+" and "-" to be returned as tokens.
        StringTokenizer termStrings = new StringTokenizer(polyStr, "+-", true);
        int degree = -1;
        boolean nextTermIsNegative = false;
        while (termStrings.hasMoreTokens()) {
            String termToken = termStrings.nextToken();
            
            if (termToken.equals("-")) {
                nextTermIsNegative = true;
            } else if (termToken.equals("+")) {
                nextTermIsNegative = false;
            } else {
                // Not "+" or "-"; must be a term
            	RatTerm term = new RatTerm();
                // If termToken has "x", decompose into coeff and exponent
            	// otherwise, treat as just a coefficient
                if (termToken.contains("x")) { 
                    int xIndex = termToken.indexOf("x");
                    String c = termToken.substring(0,xIndex); // the coefficient
                    if (c.equals("")) term.coefficient = new RatNum(1); 
                    else term.coefficient = RatNum.valueOf(c.substring(0,c.length()-1));
                    String e = termToken.substring(xIndex+1); // the exponent
                    if (e.equals("")) term.exponent = 1;
                    else term.exponent = Integer.parseInt(e.substring(1));
                }
                else { // Token is the 0-power term
                	term.coefficient = RatNum.valueOf(termToken);
                	term.exponent = 0;
                }
                // Skip the terms with 0-coefficients
                if (term.coefficient.compareTo(RatNum.ZERO) == 0) continue;
                // Record the degree of the polynomial. Test succeeds only for first (highest-power) term.
                if (term.exponent > degree) degree = term.exponent;
                // at this point, coeff and expt are initialized.
                // Need to fix coeff if it was preceeded by a '-'
                if (nextTermIsNegative) {
                	term.coefficient = term.coefficient.negate();
                }

                // accumulate terms of polynomial in 'parsedTerms'
                parsedTerms.add(term);
            }
        }
        // If degree = -1, then all terms had 0 coefficients. Return the 0 poly
        if (degree == -1) return RatPoly.ZERO;
        // Otherwise, construct the array of coefficients from the parsedTerms list
        RatNum[] coefficients = new RatNum[degree+1];
        // initializes the coefficients to 0
        for (int i=0; i<coefficients.length; i++) {
        	coefficients[i] = RatNum.ZERO;
        }
        for (int i=0; i<parsedTerms.size(); i++) {
        	RatTerm term = parsedTerms.get(i);
        	coefficients[term.exponent] = term.coefficient;
        }
        
        return new RatPoly(coefficients);
    }

    /**
     * Standard hashCode function.
     *
     * @return an int that all objects equal to this will also return.
     */
    @Override
    public int hashCode() {
        // all instances that are NaN must return the same hashcode;
        if (this.isNaN()) {
            return 0;
        }
        return Arrays.hashCode(coeffs);
    }

    /**
     * Standard equality operation.
     *
     * @param obj The object to be compared for equality.
     * @return true if and only if 'obj' is an instance of a RatPoly and 'this'
     *         and 'obj' represent the same rational polynomial. Note that all
     *         NaN RatPolys are equal.
     */
    @Override
    public boolean equals(/*@Nullable*/ Object obj) {
        if (obj instanceof RatPoly) {
            RatPoly rp = (RatPoly) obj;

            // special case: check if both are NaN
            if (this.isNaN() && rp.isNaN()) {
                return true;
            } else {
                return Arrays.equals(coeffs,rp.coeffs);
            }
        } else {
            return false;
        }
    }

    /**
     * Checks that the representation invariant holds (if any).
     */
    // Throws a RuntimeException if the rep invariant is violated.
    private void checkRep() throws RuntimeException {
        if (coeffs == null) {
            throw new RuntimeException("coeffs == null");
        }
        if (coeffs.length == 0) {
        	if (degree != 0) {
        		throw new RuntimeException("Degree of the 0 polynomial is not 0");
        	}
        }
        else {
        	if (coeffs.length-1 != degree) {
        		throw new RuntimeException("degree != coeffs.length-1");    
        	}
 		    for (int i = 0; i < coeffs.length; i++) {
		        if (coeffs[i] == null) {
		        	throw new RuntimeException("coefficient "+i+" is null");
		        }
 		    }
		    if (coeffs[degree].compareTo(RatNum.ZERO) == 0) {
		    	throw new RuntimeException("coeffs[degree] is 0");
		    }
        }
    }

}
