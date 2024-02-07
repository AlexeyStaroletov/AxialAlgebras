##### Calculations in Octonion and Albert algebras
##### for details see https://arxiv.org/abs/2305.10958

bas:=
[
[1,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0],
[0,0,0,1,0,0,0,0],
[0,0,0,0,1,0,0,0],
[0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,1]
];   ##### 8-dimensional space

ii:=bas[1];  ### identity element
i0:=bas[2];  ### i0^2=-1
i1:=bas[3];  ### i1^2=-1
i2:=bas[4];  ### i2^2=-1
i3:=bas[5];  ### i3^2=-1 
i4:=bas[6];  ### i4^2=-1
i5:=bas[7];  ### i5^2=-1
i6:=bas[8];  ### i6^2=-1

mas:=[ [1,2,3,4,5,6,7,8],
       [2,-1,5,8,-3,7,-6,-4],
       [3,-5,-1,6,2,-4,8,-7],
       [4,-8,-6,-1,7,3,-5,2],
       [5,3,-2,-7,-1,8,4,-6],
       [6,-7,4,-3,-8,-1,2,5],
       [7,6,-8,5,-4,-2,-1,3],
       [8,4,7,-2,6,-5,-3,-1]
];    ### table of products

Prod:=function(u,v)    ##### this function finds the product of two octonions in the basis ii,i0,...,i6
local ans, n, sgn, i, j;

ans:=[0,0,0,0,0,0,0,0];

for i in [1..8] do
for j in [1..8] do
  sgn:=1;
  n := mas[i][j];
  if (n < 0) then
    n := -n;
    sgn := -1;
  fi;
  ans:= ans + sgn * u[i] * v[j] * bas[n];
od;
od;
return ans;
end;

Conj:=function(u)  ##### this function returns the conjugate of an octonion
local ans, i;
  ans:=[0,0,0,0,0,0,0,0];
  ans[1]:=u[1];
  for i in [2..8] do
    ans[i] := -u[i];
  od;
  return ans;
end;


ProdAlbert := function(u, v)  ##### this function finds the product of two matrices in the Albert algebra constructed by the octonion algebra
                          ##### we follow Section 4.8 of the book R.A. Wilson, The finite simple groups. Graduate Texts in Mathematics, 251. Springer-Verlag London, DOI:10.1007/978-1-84800-988-2

local d1, d2, e1, e2, f1, f2, D1, D2, E1, E2, F1, F2, tmp1, tmp2, tmp3;

d1:=u[1]; e1:=u[2]; f1:= u[3];
D1:=u[4]; E1:=u[5]; F1:= u[6];
d2:=v[1]; e2:=v[2]; f2:= v[3];
D2:=v[4]; E2:=v[5]; F2:= v[6];

tmp1:= Prod(F1, Conj(F2)) + Prod(F2, Conj(F1)) + Prod(Conj(E1), E2) + Prod(Conj(E2), E1);
tmp2:= Prod(Conj(F1), F2) + Prod(Conj(F2), F1) + Prod(D1, Conj(D2)) + Prod(D2, Conj(D1));
tmp3:= Prod(Conj(D1), D2) + Prod(Conj(D2), D1) + Prod(E1, Conj(E2)) + Prod(E2, Conj(E1));

return [ d1*d2 + 1/2*tmp1[1], e1*e2 + 1/2*tmp2[1], f1*f2 + 1/2*tmp3[1], 
1/2*(Prod(Conj(F1), Conj(E2)) + Prod(Conj(F2), Conj(E1)) + (e1+f1)*D2 + (e2+f2)*D1),
1/2*(Prod(Conj(D1), Conj(F2)) + Prod(Conj(D2), Conj(F1)) + (d1+f1)*E2 + (d2+f2)*E1),
1/2*(Prod(Conj(E1), Conj(D2)) + Prod(Conj(E2), Conj(D1)) + (d1+e1)*F2 + (d2+e2)*F1)
];

end;

#### define 4 idempotents that generate the Albert algebra

idem1:=1/2*[1,1,0, 0*i1, 0*i1, i0];
idem2:=1/2*[1,0,1, 0*i1, i1, 0*i1];
idem3:=1/2*[0,1,1, i2, 0*i1, 0*i1];
idem4:=[1/9,4/9,4/9, 4/9*i4, 2/9*i3, 2/9*i6];

ProdAlbert(idem1, idem1) = idem1;
ProdAlbert(idem2, idem2) = idem2;
ProdAlbert(idem3, idem3) = idem3;
ProdAlbert(idem4, idem4) = idem4;

Convert:=function(u)   ##### this function transforms a matrix from the Albert algebra into a vector of length 27.
  local ans, i, j;
  ans:=[];
  Add(ans, u[1]); Add(ans, u[2]); Add(ans, u[3]);
  for i in [4..6] do
    for j in [1..8] do
      Add(ans, u[i][j]);
    od;
  od;
  return ans;
end;

flag:=true;
L:=[Convert(idem1), Convert(idem2), Convert(idem3), Convert(idem4)];
res:=[idem1, idem2,idem3,idem4];
while flag = true do
  flag:=false;
  V:=VectorSpace(Rationals, L);
  for i in [1..Size(L)] do
  if flag = true then
    break;
  fi;
  for j in [1..Size(L)] do
    x:=ProdAlbert(res[i], res[j]);
    if not Convert(x) in V then
      flag:=true;
      Add(L, Convert(x));
      Add(res, x);
      break;
    fi;
  od;od;
od;

##### calculate 27 elements of the Albert algebra

a:=idem1;
b:=idem2;
c:=idem3;
d:=idem4;
ab:=ProdAlbert(a,b);
ac:=ProdAlbert(a,c);
ad:=ProdAlbert(a,d);
bc:=ProdAlbert(b,c);
bd:=ProdAlbert(b,d);
cd:=ProdAlbert(c,d);
abc:=ProdAlbert(a,bc);
bac:=ProdAlbert(b,ac);
cab:=ProdAlbert(c,ab);
abd:=ProdAlbert(a,bd);
acd:=ProdAlbert(a,cd);
bad:=ProdAlbert(b,ad);
bcd:=ProdAlbert(b,cd);
cad:=ProdAlbert(c,ad);
cbd:=ProdAlbert(c,bd);
ab_cd:=ProdAlbert(ab,cd);
ac_bd:=ProdAlbert(ac,bd);
d_abc:=ProdAlbert(d,abc);
d_bac:=ProdAlbert(d,bac);
a_bcd:=ProdAlbert(a,bcd);
ab_cad:=ProdAlbert(ab,cad);
ab_cbd:=ProdAlbert(ab,cbd);
ac_bcd:=ProdAlbert(ac,bcd);

##### write 27 elements into one matrix

L:=[Convert(a), Convert(b), Convert(c), Convert(d),
    Convert(ab), Convert(ac), Convert(ad), Convert(bc), Convert(bd), Convert(cd),
    Convert(abc), Convert(abd), Convert(acd), Convert(bac), Convert(bad), Convert(bcd), Convert(cab), Convert(cad), Convert(cbd),
    Convert(ab_cd), Convert(ac_bd), Convert(d_abc), Convert(d_bac), Convert(a_bcd), Convert(ab_cad), Convert(ab_cbd), Convert(ac_bcd)];

    
###### verify the equality for the determinant of L

Determinant(L) = 1/(2^78*3^36);
