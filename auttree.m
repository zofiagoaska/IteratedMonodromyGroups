// Automorphism of a binary tree - magma code

AutTree := function(n)
    
// Return the automorphism group of a full binary tree of height n
// acting naturally on the 2^n leaves (numbered 1..2^n)

    if n lt 1 then
        error "n must be >= 1";
    end if;

    S := SymmetricGroup(2^n);
    gens := [];

    // Each level m = 1..n corresponds to subtree swaps at that level

    for m in [1..n] do
        for i in [1..2^(n-m)] do

            // Subtree swap covering a block of size 2^(m-1)

            start := (i-1)*2^m + 1;
            sigma := Identity(S);
            for j in [0..2^(m-1)-1] do
                sigma *:= S!(start + j, start + j + 2^(m-1));
            end for;
            Append(~gens, sigma);
        end for;
    end for;

    return sub< S | gens >;
end function;
