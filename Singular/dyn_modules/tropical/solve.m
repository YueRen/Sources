declare attributes RngInt: debug;

ZZ := Integers();
ZZ`debug := [* *];

intrinsic get_points(F, p:All := true, bound := Infinity()) -> Any
  {}
  //"Computing points with precision ", p;
  //"F contains ", #F;
  Ctx_input := Parent(F[1]);
  Ct := CoefficientRing(Ctx_input);
  R := PolynomialRing(Ct);

  h := hom< Ctx_input -> R | [ R.1 : i in [1..Rank(Ctx_input)]]>;

  g := F[1];
  
  //if #F eq 1 then
    ////"F contains h(g)", h(g);
    ////"precision", p div 2;
    //ZZ := Integers();
    //Append(~ZZ`debug, < h(g), p >);
  //end if;

  roots := PuiseuxExpansion(h(g), p - 3);

  if not All then
    roots := [ [roots[1]] ];
  end if;

  if #F eq 1 then
    return [ [r] : r in roots | Abs(Valuation(r)) le bound ];
  end if;

  Fdash := [];

  Ctx_new := PolynomialRing(Ct, Rank(Ctx_input) - 1);
  
  VF := [];
  for r in roots do
    replace_x1 := hom< Ctx_input -> Ctx_new | [ Ctx_new.i : i in [1..Rank(Ctx_new)]] cat [ r ]>;
    Fdash := [ replace_x1(f) : f in F[2..#F]];
    //"Fdash", Fdash;
    succeeded := false;
    while not succeeded do
      if p lt 0 then
        error("adsaD");
      end if;
      try
        VFdash := get_points(Fdash, p - 1: All:=All, bound := bound);
        succeeded := true;
      catch e
        p := p - 1;
      end try;
    end while;
        
    VF := VF cat [ oldpoints cat [ r ] : oldpoints in VFdash];
  end for;

  return VF;
end intrinsic;

intrinsic get_points(F:All:=true, bound := Infinity()) -> Any
  {}
  p := 4*Rank(Parent(F[1]));
  succeeded := false;

  while not succeeded do
    try
      P := get_points(F, p:All:=All, bound := bound);
      if #P eq 0 then
        succeeded := true;
      end if;
      for i in [1..#P] do
        point := P[i];
        for coord in point do
          if RelativePrecision(coord) eq 0 then
            break i;
          end if;
        end for;
        if i eq #P then
          succeeded := true;
        end if;
      end for;
    catch e
      ii := 0;
      //"Increasing precision to ", 2*p;
    end try;
    p := 2*p;
  end while;

  return P;
end intrinsic;

intrinsic TropicalVariety(F: WithMultiplicity := true, All:= true, bound := Infinity()) -> Any
  {}
  V := get_points(F: All:=All, bound := bound);
  tropF := [ [ Valuation(x) : x in p ] : p in V];
  if WithMultiplicity then
    return Sort(tropF);
  else 
    return Sort(SetToSequence(Set(tropF)));
  end if;
end intrinsic;

intrinsic TropicalVarietyOfList(LF: WithMultiplicity := true, All:= true, bound := Infinity()) -> Any
  {}
  res := [];
  for F in LF do
    res := res cat TropicalVariety(F: WithMultiplicity := WithMultiplicity, All := All, bound := bound);
  end for;
  return res;
end intrinsic;

intrinsic TropicalVarietiesListOfLists(LLF: WithMultiplicity := true, All:= true, bound := Infinity()) -> Any
  {}
  res := [];
  for i in [1..#LLF] do
    LF := LLF[i];
    Append(~res, TropicalVarietyOfList(LF: WithMultiplicity := WithMultiplicity, All := All, bound := bound));
  end for;
  n := #res div 2;
  rus := [];
  for i in [1..#res] do
    t := [];
    S := res[i];
    for j in [1..#S] do
      if i le #res/2 then
        Append(~t, Insert(S[j], i, 1));
      else
        Append(~t, Insert(S[j], i - n, -1));
      end if;
    end for;
    Append(~rus, t);
  end for;
  rus := &cat(rus);
  rus := Sort(SetToSequence(Set(rus)));
  return rus;
end intrinsic;

intrinsic PrintForSingularList(F, SS) -> Any
  {}
  toprint := "list TT;\n";
  for i in [1..#SS] do
    S := SS[i];
    if #S eq 0 then
      continue;
    end if;
    
    toprint := toprint * "matrix T[" * Sprint(#S) * "][" * Sprint(#S[1]) * "] = \n";

    for P in S do
      for p in P do
        toprint := toprint * Sprint(p) * ", ";
      end for;
      toprint := toprint * "\n";
    end for;

    toprint := Substring(toprint, 1, #toprint - 3) * ";\n";

    toprint := toprint * "TT[" * Sprint(i) * "] = T;\nkill T;\n";
  end for;

  PrintFile(F, toprint:Overwrite:=true);

  return true;
end intrinsic;

intrinsic PrintForSingular(F, S) -> Any
  {}
  pp := [ i : i in [1..#S] | #S[i] ne 0];

  if #pp eq 0 then
    toprint := "matrix TT;";
  else
    toprint := "matrix TT[" * Sprint(#pp) * "][" * Sprint(#S[pp[1]]) * "] = \n";

    for P in S[pp] do
      for p in P do
        toprint := toprint * Sprint(p) * ", ";
      end for;
      toprint := toprint * "\n";
    end for;

    toprint := Substring(toprint, 1, #toprint - 3) * ";";
  end if;

  PrintFile(F, toprint:Overwrite:=true);
  return true;
end intrinsic;
