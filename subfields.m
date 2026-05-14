//Quadratic extensions of Q(t) contained in Q(t)(f^(4)(x)-t) - magma code

P<t> := FunctionField(Rationals());
R<x> := PolynomialRing(P);

// f = g_4 - t*h_4
// g_4, h_4 - numerator and denominator of 4-th iterate of the function 2/(x-1)^2

f := 2*x^16 - 32*x^15 + 272*x^14 - 1568*x^13 + 6712*x^12 - 22304*x^11 + 58736*x^10 - 123424*x^9 + 205964*x^8 - 268128*x^7 + 263024*x^6 - 182112*x^5 + 78392*x^4 - 16736*x^3 + 1808*x^2 - 96*x + 2 - t*(x^16 - 16*x^15 + 72*x^14 + 112*x^13 - 1892*x^12 + 5232*x^11 + 1272*x^10 - 39824*x^9 + 107334*x^8 - 153520*x^7 + 138744*x^6 - 85296*x^5 + 35996*x^4 - 9392*x^3 + 1352*x^2 + 80*x + 1);

G, r := GaloisGroup(f);
M := MaximalSubgroups(G);
n := #M;
S<t> := PolynomialRing(Integers());
for i in [1..n] do
	// p - defining polynomial of a subfield corresponding to M[i] subgroup
	// a - M[i] invariant

	p, a := GaloisSubgroup(f,M[i]`subgroup);
	
	D := Discriminant(p);
	
	// For quadratic polynomial p we have K({x|p(x)=0}) = K (sqrt(Discr(p)))
	// We're going to determine the square-free part of Discr(p)
	
	l := S!1;
	if LeadingCoefficient(D) lt 0 then
		l := S!(-1);
	end if;
	F := Factorization(S!D);
	for k in F do
		if k[2] mod 2 eq 1 then
			l := l*k[1];
		end if;
	end for;
	print "Subfield", i;
	printf "Q(t)(sqrt(%o))\n",l;
	// R!p, a;
end for;
