##### This file contains commands for checking the linearized Jordan identity in some Matsuo algebras.
##### Calculations are used in the paper https://arxiv.org/abs/2305.10958

##### First we have a list of groups to check, choose any of them and continue.

##### 3^7:Sym(8)=Wr(3,8) 

N:=Group((1,2,3));
H:=SymmetricGroup(8);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,9), (1,2));
G:=NormalClosure(G, [a]);

##### 3^6:Sym(7)=Wr(3,7)

N:=Group((1,2,3));
H:=SymmetricGroup(7);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,8), (1,2));
G:=NormalClosure(G, [a]);

##### 3^5:Sym(6)=Wr(3,6)

N:=Group((1,2,3));
H:=SymmetricGroup(6);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,7), (1,2));
G:=NormalClosure(G, [a]);

##### 3^4:Sym(5)=Wr(3,5)

N:=Group((1,2,3));
H:=SymmetricGroup(5);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,6), (1,2));
G:=NormalClosure(G, [a]);

##### 3^3:Sym(4)=Wr(3,4)

N:=Group((1,2,3));
H:=SymmetricGroup(4);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,5), (1,2));
G:=NormalClosure(G, [a]);

##### 2^7:Sym(8)=Wr(2,8)

N:=Group((1,2));
H:=SymmetricGroup(8);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,9), (1,2));
G:=NormalClosure(G, [a]);

##### 2^6:Sym(7)=Wr(2,7)

N:=Group((1,2));
H:=SymmetricGroup(7);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,8), (1,2));
G:=NormalClosure(G, [a]);

##### 2^5:Sym(6)=Wr(2,6)

N:=Group((1,2));
H:=SymmetricGroup(6);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,7), (1,2));
G:=NormalClosure(G, [a]);

##### 2^4:Sym(5)=Wr(2,5)

N:=Group((1,2));
H:=SymmetricGroup(5);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,6), (1,2));
G:=NormalClosure(G, [a]);

##### 2^3:Sym(4)=Wr(2,4)

N:=Group((1,2));
H:=SymmetricGroup(4);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,5), (1,2));
G:=NormalClosure(G, [a]);

##### Alt(4)^4:Sym(4)=Wr(Alt(4),4)

N:=AlternatingGroup(4);
H:=SymmetricGroup(4);
G:=WreathProduct(N,H);
a:=Image(Embedding(G,5), (1,2));
G:=NormalClosure(G, [a]);

##### O^+(8,2)

G:=Image(IsomorphismPermGroup(GO(1,8,2)));
cc:=ConjugacyClasses(G);;
ind:=-1;
for i in [1..Size(cc)] do
  Representative(cc[i]);
  if Order(Representative(cc[i])) = 2 and Size(cc[i]) = 120 then
    ind := i;
    break;
  fi;  
od;
a:=Representative(cc[ind]);

##### O^-(6,2)

G:=Image(IsomorphismPermGroup(GO(-1,6,2)));
cc:=ConjugacyClasses(G);;
ind:=-1;
for i in [1..Size(cc)] do
  Representative(cc[i]);
  if Order(Representative(cc[i])) = 2 and Size(cc[i]) = 36 then
    ind := i;
    break;
  fi;  
od;
a:=Representative(cc[ind]);

##### Sp(6,2)
G:=Sp(6,2);
cc:=ConjugacyClasses(G);;
ind:=-1;
for i in [1..Size(cc)] do
  Representative(cc[i]);
  if Order(Representative(cc[i])) = 2 and Size(cc[i]) = 63 then
    ind := i;
    break;
  fi;  
od;
a:=Representative(cc[ind]);

##### ^+Omega^-6(3)

G:=Image(IsomorphismPermGroup(GO(-1,6,3)));
cc:=ConjugacyClasses(G);;
ind:=-1;
for i in [1..Size(cc)] do
  Representative(cc[i]);
  if Order(Representative(cc[i])) = 2 and Size(cc[i]) = 126 then
    ind := i;
    break;
  fi;  
od;
a:=Representative(cc[ind]);
G:=NormalClosure(G, [a]);

##### SU(4,2)
G:=SU(4,2);
cc:=ConjugacyClasses(G);;
ind:=-1;
for i in [1..Size(cc)] do
  Representative(cc[i]);
  if Order(Representative(cc[i])) = 2 and Size(cc[i]) = 45 then
    ind := i;
    break;
  fi;  
od;
a:=Representative(cc[ind]);

##### SU(5,2)
G:=Image(IsomorphismPermGroup(SU(5,2)));
cc:=ConjugacyClasses(G);;
ind:=-1;
for i in [1..Size(cc)] do
  Representative(cc[i]);
  if Order(Representative(cc[i])) = 2 and Size(cc[i]) = 165 then
    ind := i;
    break;
  fi;  
od;
a:=Representative(cc[ind]);

##### 4^3:SU(3,2)' (we use representation of this group from J.I. Hall, L.H. Soicher, Presentations of some 3-transposition groups, Comm. Algebra, 23:7 (1995), 2517–2559, DOI:10.1080/00927879508825358)

F:=FreeGroup(4);
a:=F.1;; b:=F.2;; c:=F.3;; d:=F.4;;
rel:=[a^2, b^2, c^2, d^2, (a*b)^3, (a*c)^3, (a*d)^2, (b*c)^3, (b*d)^3, (c*d)^3, (b^a*c)^3, (b^d*c)^3, (a^b*c^d)^3] ;
G:=F/rel;
H:=Group(G.1,G.4, (G.1*G.2*G.3)^2, (G.1*G.2*G.3^G.4)^2 );
n:=Index(G,H);
Add(listOfInd, n);
iso:=FactorCosetAction(G,H);
a:=Image(iso,G.1);;
b:=Image(iso,G.2);;
c:=Image(iso,G.3);;
d:=Image(iso,G.4);;
G:=Group(a,b,c,d);
cc:=ConjugacyClasses(G);;
ind:=-1;
for i in [1..Size(cc)] do
  Representative(cc[i]);
  if Order(Representative(cc[i])) = 2 and Size(cc[i]) = 36 then
    ind := i;
    break;
  fi;  
od;
a:=Representative(cc[ind]);

##### POmega(8,2):Sym(3) (we use representation of this group from J.I. Hall, L.H. Soicher, Presentations of some 3-transposition groups, Comm. Algebra, 23:7 (1995), 2517–2559, DOI:10.1080/00927879508825358)

F:=FreeGroup(5);
a:=F.1;; b:=F.2;; c:=F.3;; d:=F.4;; e:=F.5;;
rel:=[a^2, b^2, c^2, d^2, e^2, (a*b)^3, (a*c)^3, (a*d)^2, (a*e)^2, (b*c)^3, (b*d)^3, (b*e)^2, (c*d)^3, (c*e)^3, (d*e)^2, 
(b^a*c)^3, (b^d*c)^3, (a^b*c^d)^3, (d^((a*b*c)^2)*e)^3] ;
G:=F/rel;
H:=Group(G.5, G.1, G.2, G.4, (G.5*G.3*G.2*G.4*G.3)^2, (G.5*G.3*G.2*G.1*G.3)^2);
n:=Index(G,H);
iso:=FactorCosetAction(G,H);
a:=Image(iso,G.1);;
b:=Image(iso,G.2);;
c:=Image(iso,G.3);;
d:=Image(iso,G.4);;
e:=Image(iso,G.5);;
G:=Group(a,b,c,d,e);

##### Common code for all groups

##### the following function verifies that a list of involutions is a set of 3-transpositions

check := function(class)
 local i,a;
 a:=Representative(class);
 for i in class do
   if Order(a*i) > 3 or not Order(i) = 2 then
     Print("False\n");
     return;
   fi; 
 od;
 Print("True\n");
 return;
end;

#####

F:=Rationals;       ### we work over the field of rationals
eta:=1/2;           ### parameter in Matuso algebra
inv:=Orbit(G,a);    ### conjugacy class of 3-transpositions
check(inv);         ### verify the 3-transposition property for inv
V:=F^(Size(inv));
L:=[];
for i in [1..Size(inv)] do          ### standard basis on involutions from inv for the Matsuo algebra
  q:=List([1..Size(inv)], k->0);
  q[i]:=1;
  Add(L, q);
od;

Gram:=[];                           ### build the Gram matrix 
for i in [1..Size(inv)] do
  Add(Gram, []);
  for j in [1..Size(inv)] do
    if i = j then
      Add(Gram[i], 1);
    fi;
    if (Order(inv[i]*inv[j]) = 3) then
      Add(Gram[i], 1/4);
    fi;
    if (Order(inv[i]*inv[j]) = 2) then
      Add(Gram[i], 0);
    fi;
  od;  
od;

bas:=NullspaceMat(Gram);               ### find a basis of the radical of the Frobenius form
if Size(bas) = 0 then
  W:=F^0;
else
  W:=VectorSpace(F, bas);
fi;

num := [];

for i in [1..Size(inv)] do                  ## num contains for each pair of non-commuting involutions inv[i] and inv[j] the number of the involuion inv[i]*inv[j]*inv[i] in inv
    Add(num, List([1..Size(inv)], k->0));   ## we will use it to calculate the product in the Matsuo algebra
    for j in [1..Size(inv)] do
        Add(num[i], 0);
        if Order(inv[i] * inv[j]) = 3 then
            tmp := inv[i] * inv[j] * inv[i];
            for k in [1..Size(inv)] do
            if (inv[k] = tmp) then
                num[i][j] := k;
                break;
            fi;
            od;
        fi;
    od;
od;

prod := function (a,b)                  ## function for the product in Matsuo algebra: a and b are vectors in the basis of involutions (inv)

    local i,j,n, ord, q;
    
    n := Size(a);
    q := List([1..n], i->0);
  
    for i in [1..n] do
    for j in [1..n] do
        if not(a[i] = 0) and not(b[j] = 0) then
            ord := Order(inv[i] * inv[j]);
            if ord = 1 then
                q[i] := q[i] + a[i] * b[j];
            fi;
            if ord = 3 then
                q[i] := q[i] + a[i] * b[j] * eta/2;
                q[j] := q[j] + a[i] * b[j] * eta/2;
                q[num[i][j]] := q[num[i][j]] - a[i] * b[j] * eta/2;
            fi;
            
        fi;
        
    od;
    od;
       
    return q;
        
end;

assos := function(a,b,c)            ### this function calculates the associator for a, b, c
  return prod(a, prod(b,c)) - prod(c, prod(a,b));
end;

res:=true;

x:=L[1];                      ##### now we verify the linearized Jordan identity on elements from inv. Since inv is a conjugacy class, we can fix the first element x.
for j in [1..Size(inv)] do
for k in [1..Size(inv)] do
for t in [1..Size(inv)] do
 y:=L[j];
 z:=L[k];
 w:=L[t];
 if not ( assos(prod(x,z), y, w) + assos(prod(z,w), y, x) + assos(prod(w,x), y, z) in W) then
   res:=false;
 fi;
od;od;od;

res;  ##### this variable tells us the result of verification.

