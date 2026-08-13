// Generators of G - magma code

Geom := function(n)
    if n le 12 then
        k := 12;
        S := SymmetricGroup(2^n);
        m := 1;

        // Set here S!(1,2) if automorphism acts non-trivially on first level and set S!(1) if it acts trivially on the first level.

        a1 := S!(1,2); a2 := S!(1,2); a3 := S!(1); 

        // Define s to be the element that switches the first half and the second half of the tree on level n that is switching {1,2,...,2^{n-1}} with {2^{n-1}+1,...,2^n} in that order. 
        // s helps us later to define the new level with the help of knowing how the automorphism work on one higher level. 

        m := 2;
        while m in [2..n] do
            s := S!(1); L := [];
            for i in [1..2^(m-1)] do
                L[i] := S!(i, i + 2^(m-1));
            end for;
            for i in [1..2^(m-1)] do
                s := s * L[i];
            end for;

	    // save temporarily what the a_i do at one level before

            b1 := a1; b2 := a2; b3 := a3;

	    // now if a_i= (g,h)tau put here g*s*h*s*tau; and sigma = s and have s*s=id. 
	    // E.g. a_2 = ((a_3)^{-1}, (a_2)^{-1})sigma, then change a_i from one level higher to b_i and use formula to get
	    // a_2 = (b_3)^{-1} * s * (b_2)^{-1} * s * s   =  (b_3)^{-1} * s * (b_2)^{-1}.

	    a3 := b2 * s * b3 * s;
            a2 := Inverse(b3)*s*Inverse(b2);
            a1 := s;
            m := m + 1;
        end while;
        return [a1, a2, a3];
    end if;
end function;
