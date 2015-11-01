################################################################################
# spectra.gap
# Briana Foster-Greenwood
# Last update: 28 Oct 2015

################################################################################
# This code is written for use in GAP 3 with the CHEVIE package.  At the
# time of writing, CHEVIE was not available for GAP 4.  The CHEVIE 
# package contains functionality for working with complex reflection groups.
#
# For information on GAP: http://www.gap-system.org/gap.html
#                         http://webusers.imj-prg.fr/~jean.michel/gap3/

################################################################################
# This code accompanies the preprint:
# B. Foster-Greenwood and C. Kriloff. Spectra of Cayley graphs of complex 
# reflection groups. http://arxiv.org/abs/1502.07392.  Submitted, 2015.
#-------------------------------------------------------------------------------
# Background (references given in article): Given a class function f on a 
# finite group G, define a matrix with |G| rows and columns indexed by the 
# elements of G (in some ordering).  Let the (g,h)-entry be the value 
# f(gh^{-1}).  The eigenvalues of this matrix can be calculated via a character
# theoretic formula:
#         |G|/X(1) * <f,X> with multiplicity X(1)^2
# where X ranges over the irreducible characters of G.  If f is integer-valued
# and constant on rational conjugacy classes, then the eigenvalues are
# guaranteed to be integers.  (Group elements g and h are rationally conjugate
# if the cyclic subgroups <g> and <h> are conjugate.  Equivalently, h must
# be conjugate to a generator g^d of <g>, where d is necessarily relatively
# prime to the order of g.)
#-------------------------------------------------------------------------------
# In the article, we are interested in the case where G is a complex 
# reflection group and f is a class function that determines a matrix
# related to a Cayley graph for G.  If T is the set of all reflections in G, 
# then the left Cayley graph of G with respect to T has a vertex for each 
# element of the group and an edge joining h to g if gh^{-1} is a reflection.  
# (In this case, hg^{-1} is also a reflection, so there is an edge from g to h
# as well.)
#
# The following class functions lead to various matrices relating to the Cayley
# graph of G with respect to T:
# (a) REFLECTION INDICATOR FUNCTION
#     g -> 1 if g is a reflection and g -> 0 otherwise
#     gives rise to ADJACENCY MATRIX of the Cayley graph 
# (b) (ABSOLUTE) REFLECTION LENGTH FUNCTION
#     g -> min number of factors to express g as a product of reflections
#     gives rise to DISTANCE MATRIX of the Cayley graph 
# (c) CODIMENSION FUNCTION
#     g -> codimension of the fixed point space of g
#     gives rise to CODIMENSION MATRIX (with respect to group representation)
# 
# The purpose of this code is to calculate eigenvalues of the adjacency,
# distance, and codimension matrices via the character formula and to
# determine if absolute reflection length is constant on rational conjugacy
# classes.  The files 
#    spectraprompt.gap
#    constantprompt.gap
# call the relevant functions from this file spectra.gap. 
#-------------------------------------------------------------------------------

################################################################################
################################################################################
############                    CLASS FUNCTIONS                     ############
################################################################################ 
################################################################################

# This section is for computing the reflection indicator function,
# absolute reflection length function, and codimension function.

################################################################################
#                                 GET REF REP
# Input:  GROUP     - irreducible complex reflection group
#         CHARTABLE - character table of GROUP
#
# Output: position of a reflection representation in the list 
#         CHARTABLE.irreducibles of irreducible characters of the reflection
#         group
#-------------------------------------------------------------------------------
GetRefRep:=function(group,charTable)

   local i, refRepCharValueList;
   
   # use built-in function to get list of character values for reflection rep
   refRepCharValueList:=ReflectionCharacter(group);
		 
   # use inner-products of characters to find position in character table	 
   for i in [1..Length(charTable.irreducibles)] do
   
      if ScalarProduct(ClassFunction(charTable,charTable.irreducibles[i]),
	               ClassFunction(charTable,refRepCharValueList))=1 then
	  
	     return i;

      fi;
	  
   od;

end;

################################################################################
#                                 IRR  CODIM
# Input:  CHARTABLE - the character table of a group G
#         CHARNMBR - position of an irreducible character in 
#                    charTable.irreducibles
#         CLASSNMBR - position of a conjugacy class in charTable.classes
#
# Output:  codimension of subspace fixed by an elt in the conjugacy class
#          given by CLASSNMBR in representation with specified character
#
# WARNING:  THIS FUNCTION ASSUMES THE IDENTITY IS IN THE FIRST COLUMN
#           OF THE CHARTABLE
#-------------------------------------------------------------------------------
IrrCodim:=function(charTable,charNmbr,classNmbr)

   local evals;

   evals:=Eigenvalues(charTable,charTable.irreducibles[charNmbr],classNmbr);
   
   return Eigenvalues(charTable,charTable.irreducibles[charNmbr],1)[1]
             -evals[Length(evals)];

end;

################################################################################
#                              CODIM CLASS FUNCTION
# Input:  CHARTABLE - character table for a group
#         CHARNMBR - position of an irreducible character in
#                    charTable.irreducibles
# 
# Output:  a class function with g -> codim(V^g), where V^g=subspace fixed by g
#-------------------------------------------------------------------------------
CodimClassFunction:=function(charTable,charNmbr)

   local c;

   return ClassFunction(charTable,List([1..Length(charTable.irreducibles)],
             c->IrrCodim(charTable,charNmbr,c)));

end;

################################################################################
#                               REFLECTION CLASSES
# Input:   CHARTABLE - character table for a reflection group
#          CHARNMBR - position of a reflection representation in the
#                     the list CHARTABLE.irreducibles
#
# Output:  list of positions in CHARTABLE.classes corresponding to conjugacy
#          classes of reflections
#-------------------------------------------------------------------------------
RefClasses:=function(charTable,charNmbr)

   local classNmbr;

   # reflections are the elements with a codimension one fixed point space
   return Filtered([1..Length(charTable.irreducibles)],
                    classNmbr->IrrCodim(charTable,charNmbr,classNmbr)=1);

end;

################################################################################
#                          REFLECTION CLASS FUNCTION
# Input: CHARTABLE - character table of a reflection group
#        CHARNMBR - position of a reflection representation in 
#                   CHARTABLE.irreducibles
#
# Output: a class function with g->1 if g is a reflection
#                               g->0 if g is not a reflection
#-------------------------------------------------------------------------------
RefClassFunction:=function(charTable,charNmbr)

   local charfunctionrefs, refs, r, i;

   charfunctionrefs:=List([1..Length(charTable.irreducibles)],i->0);
   
   refs:=RefClasses(charTable,charNmbr);
   
   for r in refs do
   
      charfunctionrefs[r]:=1;
	  
   od;
   
   return ClassFunction(charTable,charfunctionrefs);

end;

################################################################################
#                              LENGTH PARTITION
# Input:  CHARTABLE - character table of a group G
#         GENCLASSNUMS - conjugacy class numbers for a set of generators of
#                        the group G  (ASSUMING 1 is not in GENCLASSNUMS)
#
# Output: a partition of the conjugacy classes of group G according to length 
#         with respect to the conjugacy classes given in GENCLASSNUMS
#
#         Definition of length(g):
#         length(g)=min number of elts needed to express g as a product of
#                   elts from the conjugacy classes in GENCLASSNUMS
#         This length function is constant on conjugacy classes.
#
#         The partition is returned as a list.  Each element of the list is
#         itself a list:  the first entry is a length (possibly "infinity"); 
#         the second entry is a list of class numbers of elts with that 
#         length.
#-------------------------------------------------------------------------------
LengthPartition:=function(charTable, genClassNums)

   local P, remainingClasses, i, classi, classj, classk, 
         previouslength, currentlength;

   # P will be the partition:  
   #     [0,[1]] says the identity has length zero
   #     [1, genClassNums] says the generating classes have length 1
   # IT IS ASSUMED THAT GENCLASSNUMS WAS PASSED IN NOT CONTAINING THE IDENTITY
   P:=[ [0,[1]], [1,genClassNums] ];

   # list all the class numbers [1,2,3,...]
   remainingClasses:=Set(List([1..Length(charTable.irreducibles)],i->i));
   
   # remove 1 and all the genClassNums from remainingClasses 
   #   (since their length is already known)
   # note: SubtractSet is destructive
   SubtractSet(remainingClasses,P[1][2]);
   SubtractSet(remainingClasses,P[2][2]);
   
   previouslength:=1;
   currentlength:=2;   
   
   while Length(P[previouslength+1][2])>0 and Size(remainingClasses)>0 do
   
      # There were some classes with previouslength, so now we need to look for
	  # classes of currentlength, which will be stored in slot 
	  # currentlength+1 of P
   
      Add(P,[currentlength,[]]);
   
      for classi in genClassNums do
	  
	  	 for classj in P[previouslength+1][2] do
		  
		    for classk in remainingClasses do
                			   
               if (not classk in P[currentlength+1][2]) and 
		         ClassMultCoeffCharTable(charTable,classi,classj,classk)>0 then
				  
                 # can multiply an element from generating class with
                 # an element with previouslength and get an element in classk,
                 # so classk has currentlength				 
				 AddSet(P[currentlength+1][2],classk);
					 
			   fi;
			     
			od;
			  
		 od;
		  
      od;
	   
	  SubtractSet(remainingClasses,P[currentlength+1][2]);
	   
	  previouslength:=currentlength;
	  currentlength:=currentlength+1;
	   
   od;
   
   if Length(P[previouslength+1][2])=0 and Size(remainingClasses)>0 then
   
	  P[previouslength+1][1]:="infinity";
	  AddSet(P[previouslength+1][2],remainingClasses);
	  
   fi;

   return P;

end;

################################################################################
#                             LENGTH CLASS FUNCTION
# Input:  CHARTABLE - character table of a group
#         GENCLASSNUMS - positions of conjugacy classes whose elements
#                        generate the group (the identity should NOT be in this
#                        list)
#
# Output: a class function with g->length(g)
#         [length(g)=min number of elts needed to express g as a product of
#                    elts from the conjugacy classes in GENCLASSNUMS          ]
#         if the classes specified in GENCLASSNUMS do not generate the group,
#         an error message is printed and nothing is returned
#-------------------------------------------------------------------------------
LengthClassFunction:=function(charTable,genClassNums)

   local lengthlist, i, P, Pd, c;
   
   # make a list to store lengths of elts in each conj class
   lengthlist:=List([1..Length(charTable.irreducibles)],i->-1);
   
   P:=LengthPartition(charTable,genClassNums);
   
   # each Pd in the length partition P is a list [d,[c_1,...,c_r]]
   #   -the first slot gives a length d
   #   -the second slot is the list of positions of the conjugacy classes
   #    having length d
   # in case a class is not generated, it will be in a slot ["infinity",[...]]
   for Pd in P do
   
      for c in Pd[2] do
	  
	     lengthlist[c]:=Pd[1];
		 
      od;
   
   od;
   
   # < evaluates to true when objects of different types are compared
   if Maximum(lengthlist) < "infinity" then
   
      return ClassFunction(charTable,lengthlist);
   
   fi;
   
   Print("Some classes not generated.  Could not make class function.\n");

end;

################################################################################
################################################################################
############               RATIONAL CLASS FUNCTIONS                 ############
################################################################################ 
################################################################################

# The functions in this section are for determining if a class function
# is actually constant on rational conjugacy classes.

################################################################################
#                             CLASS NUM OF POWER
# Input: CHARTABLE - the character table of a finite group
#        CLASSNUM - the position of a class in the list of conjugacy classes
#                   (as ordered in CHARTABLE)
#        POWER - a positive integer
#
# Output: Let g be a representative of the conjugacy class whose position in
#         the conjugacy class list is CLASSNUM.  This function returns the 
#         position of the conjugacy class of the element g^POWER
#
# WARNING: The p-th powermap is a list that stores in position i the 
#          position of the class containing the p-th powers of the elements 
#          in the i-th class.  
#
#          If a power map is not bound to the character table, then we
#          call PowerMap, which returns a "parametrized map", i.e., a list of 
#          possibilities for the powermap (as far as can be deduced without 
#          the group available).  If PowerMap returns more than one possibility, 
#          we default to the first option in the list and print a warning.
#          
#          In practice:  For the exceptional reflection groups, not all 
#          powermaps were bound to the character table, but when PowerMap was
#          called, the ambiguous situtation did not arise.    
#-------------------------------------------------------------------------------
ClassNumOfPower:=function(charTable,classNum,power)

  local newClassNum, primeFactorizationOfPower, i, p, powerMapCandidates;

  # reduce POWER modulo the order of the elements in conjugacy class CLASSNUM
  power:=RemInt(power,charTable.orders[classNum]);

  if power = 1 then       # g^POWER = g
 
    return classNum;      # return conjugacy class of g

  fi;

  if power = 0 then       # g^POWER = 1,

    return Position(charTable.orders,1);   # return conjugacy class of identity

  fi;

  # initialization
  newClassNum:=classNum;

  # list the prime factors of power, with multiplicities
  primeFactorizationOfPower:=FactorsInt(power);

  for i in [1..Length(primeFactorizationOfPower)] do
  
     # Get the i-th prime factor of POWER
     p:=primeFactorizationOfPower[i];

     # if p-th powermap is bound to character table, use it
     if IsBound(charTable.powermap[p]) then
        
        newClassNum:=charTable.powermap[p][newClassNum];
        
     else  # try to compute the powermap on the fly
        
        # Print("*");    # Print to track how many times PowerMap is called
        
        # Get candidates for p-th powermap
        powerMapCandidates:=Powermap(charTable,p);

        # Output a warning if the powermap is not definitively known
        if Length(powerMapCandidates) > 1 then

           Print("WARNING: Powermap ",i," not determined, using \n",
                 powerMapCandidates[1],"\n");

        fi;

        newClassNum:=powerMapCandidates[1][newClassNum];
        
     fi;

  od;

  return newClassNum;

end;

################################################################################
#                           RATIONAL CONJUGACY CLASS
# Input:  CHARTABLE - character table of a finite group
#         CLASSNUM - position of conjugacy class (as ordered in CHARTABLE)
# 
# Output: An element h is rationally conjugate to g if h is conjugate to
#         g^d for some power d relatively prime to the order of g
#         (in other words <h> and <g> are conjugate subgroups).  Rational
#         conjugacy classes form a coarser partition of the group than
#         ordinary conjugacy classes.  This function outputs a list of 
#         positions of conjugacy classes (as ordered in CHARTABLE).  The union 
#         of the corresponding conjugacy classes is the rational conjugacy class
#         of an element g in the conjugacy class with position CLASSNUM.
#-------------------------------------------------------------------------------
RationalConjClass:=function(charTable,classNum)

  local orderg, relprimeorders, d, i, ratclass;

  # Find order of the elements in CLASSNUM
  orderg:=charTable.orders[classNum];
  
  # Print("Order: ", orderg, "\n");

  # Make a list of powers relatively prime to order of g
  relprimeorders:=Filtered([1..orderg],d->Gcd(d,orderg)=1);
  
  # Print("Relatively Prime Orders: ",relprimeorders,"\n");

  # For each power relatively prime to order(g), find the conjugacy class of g^d
  ratclass:=List([1..Length(relprimeorders)],
         i->ClassNumOfPower(charTable,classNum,relprimeorders[i]));  

  return Set(ratclass);

end;

################################################################################
#                          ALL RAT CLASSES
# Input:    CHARTABLE - character table for a finite group G
#
# Output:   A partition of the conjugacy classes of G into
#           rational conjugacy classes
# 
# WARNING:  Creating the character table and determining rational
#           conjugacy classes should be done in two commands:
#
#           gap> t:=CharTable(ComplexReflectionGroup(G));
#           gap> RationalConjugacyClasses(t);
#
#           For some reason, doing everything at once via
#           gap> RationalConjugacyClasses(CharTable(ComplexReflectionGroup(G));
#           produces strange results.
# 
#--------------------------------------------------------------------
RationalConjugacyClasses:=function(charTable)

  local numClasses, ratclasses, i, nextrationalclass, remaining;

  # initialize list to hold the partition of the conjugacy classes into
  # rational conjugacy classes
  ratclasses:=[];

  # number of conjugacy classes is same as number of irreducible characters
  numClasses:=Length(charTable.irreducibles);

  # this is the list of conjugacy classes that are not in a block
  # of the partition yet
  remaining:=Set(List([1..Length(charTable.irreducibles)],i->i));

  while Length(remaining)>0 do

    # Get the rational conjugacy class of the first class not yet in partition
    # and add it to the list of rational conjugacy classes
    nextrationalclass:=RationalConjClass(charTable,remaining[1]);
    Add(ratclasses,nextrationalclass);
    
    # Update remaining so it only contains conjugacy classes whose rational
    # conjugacy class is not yet included in the ratclasses list.
    SubtractSet(remaining,nextrationalclass);  # destructive!

  od;

  return ratclasses;

end;

################################################################################
#                            VALUES ON CLASS UNION
# Input: CHARTABLE - character table of a finite group
#        CLASSNUMS - a list of positions of conjugacy classes
#        CLASSFN   - a class function
# 
# Output: a list of the values of CLASSFN on the conjugacy classes specified
#         by CLASSNUMS
#
# Note: This function is short but naming it makes some other code easier
#       to read.
#-------------------------------------------------------------------------------
ValuesOnClassUnion:=function(charTable,classNums,classFn)

  local i;

  return List([1..Length(classNums)], i->classFn.values[classNums[i]]);

end;

################################################################################
#                          VALUES ON CLASS UNIONS
# Input: CHARTABLE   - character table of a finite group
#        CLASSUNIONS - a list of lists of conjugacy class positions
#        CLASSFN     - a class function
# 
# Output: a list of lists: if [k_1,...,k_p] is the i-th list in CLASSUNIONS,
#         then [CLASSFN(k_1),...,CLASSFN(k_p)] is the i-th list in the output.
#-------------------------------------------------------------------------------
ValuesOnClassUnions:=function(charTable,classUnions,classFn)

  local i;

  return List([1..Length(classUnions)], 
              i->ValuesOnClassUnion(charTable,classUnions[i],classFn));

end;

################################################################################
#                       IS CONSTANT ON CLASS UNION
# Input: CHARTABLE - character table of a finite group
#        CLASSNUMS - a list of positions of conjugacy classes 
#                    (as ordered in CHARTABLE)
#        CLASSFN   - a class function
#
# Output: True  - if the function CLASSFN is constant on the union of the 
#                 conjugacy classes specified by CLASSNUMS
#         False - otherwise
#--------------------------------------------------------------------
IsConstantOnClassUnion:=function(charTable,classNums,classFn)

  local valuesOnClassNums, i;

  valuesOnClassNums:=ValuesOnClassUnion(charTable,classNums,classFn);

  return Maximum(valuesOnClassNums)=Minimum(valuesOnClassNums);

end;

################################################################################
#                       IS CONSTANT ON CLASS UNIONS
# Input: CHARTABLE   - character table of a finite group
#        CLASSUNIONS - a list of lists of positions of conjugacy classes 
#        CLASSFN     - a class function
#
# Output: True - if the function CLASSFN is constant on each union of 
#                conjugacy classes in the list CLASSUNIONS
#         False - otherwise
#
# Example: If CLASSUNIONS is RationalConjugacyClasses(CHARTABLE), then this 
#          function determines if CLASSFN is constant on rational conjugacy 
#          classes.
#-------------------------------------------------------------------------------
IsConstantOnClassUnions:=function(charTable,classUnions,classFn)

  local i;

  for i in [1..Length(classUnions)] do
  
    if not IsConstantOnClassUnion(charTable,classUnions[i],classFn) then
    
       return false;
       
    fi;
    
  od;
  
  return true;

end;

################################################################################
################################################################################
############                         SPECTRUM                       ############
################################################################################ 
################################################################################

# The functions in this section are for using the character formula to 
# calculate the eigenvalues of the adjacency, distance, and codimension 
# matrices.

################################################################################
#                                  SPECTRUM
# Input:  CHARTABLE - character table of a finite group
#         CLASSFN   - a class function
#
# Output:  list of pairs [lambda,m], where lambda is an eigenvalue of the 
#          matrix and m is the multiplicity of the eigenvalue  
#          
# WARNING:  THIS FUNCTION ASSUMES THE IDENTITY IS IN THE FIRST COLUMN
#           OF THE CHARTABLE
#-------------------------------------------------------------------------------
Spectrum:=function(charTable,classFn)
 
   local charNmbr, numIrrChars, i, j, 
         spectrum, multiplicities, distinctspectra, combined;
 
   # Get number of irreducible characters (equals number of conjugacy classes)
   numIrrChars:=Length(charTable.irreducibles);

   # Each irreducible character X determines eigenvalue |G|/X(1)*<X,classFn>
   # of the matrix M_classFn
   spectrum:=List([1..numIrrChars],
     i->charTable.order/charTable.irreducibles[i][1]*
     ScalarProduct(classFn,ClassFunction(charTable,charTable.irreducibles[i])));

   # multiplicity of eigenvalue corresponding to character X is X(1)^2
   multiplicities:=List([1..numIrrChars],i->charTable.irreducibles[i][1]^2);

   # get the list of distinct eigenvalues
   distinctspectra:=Reversed(Set(spectrum));

   # list for storing final output
   combined:=[];

   for i in distinctspectra do

      # Multiple characters may determine the same eigenvalue, so
      # Filtered creates a list of all character numbers determining eval i
      # Sublist pulls out the multiplicities X(1)^2 for those characters
      # Sum adds up the character-wise multiplicities to get the total 
      #   multiplicity m
      # Add puts pair [i,m] into list of data: i=eigenvalue, m=multiplicity     
      Add(combined,[i,Sum(Sublist(multiplicities,
                                Filtered([1..numIrrChars],j->spectrum[j]=i)))]);     

   od;

   return combined;
	  
end;

################################################################################
#                              TEX FORM SPECTRUM
# Input:   SPECTRUM - a list of pairs [lambda,m]
#
# Output:  prints list of lambda^{m} separated by commas and with multiplicities 
#          m enclosed in ^{} so the output can be copy-pasted into a 
#          LaTeX document
#-------------------------------------------------------------------------------
TeXForm:=function(spectrum)

   local i;

   Print(spectrum[1][1],"^{",spectrum[1][2],"}");

   for i in [2..Length(spectrum)] do

       Print(", ",spectrum[i][1],"^{",spectrum[i][2],"}");

   od;

   Print("\n");

end;

################################################################################
#EOF
