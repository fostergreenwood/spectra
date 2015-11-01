################################################################################
# constantprompt.gap
# Briana Foster-Greenwood
# Updated: 28 Oct 2015

################################################################################
# See notes at beginning of spectra.gap for more information.

################################################################################
#                             HOW TO USE THIS FILE
# 1. Start a GAP 3 session and load the package CHEVIE (some distributions
#    automatically load the package, so you may not have to do anything).
#    gap> LoadPackage("chevie");
# 2. Modify the definition of "helperfile" (see next section) to be the 
#    path of the file spectra.gap (which you must have downloaded and saved).
# 3. Read in this file.  The input string is the path for the file.
#    gap> Read("/home/username/mygapfiles/constantprompt.gap");
# 4. A prompt will appear telling you what to type to determine if
#    reflection length is constant on rational conjugacy classes.
# 5. You can use LogTo if you want to write the results to a file:
#    gap> LogTo("/home/username/mygapfiles/output.txt");
#    And when finished:
#    gap> LogTo();

################################################################################
#                            DEFINE HELPER FILE
helperfile:="/home/username/mygapfiles/spectra.gap";
Read(helperfile);

################################################################################
#                             STARTUP MESSAGE
Print("\n",
       "To check whether absolute reflection length is\n",
       "constant on rational conjugacy classes for group G_k:\n",
       "  Type: IsLengthRationalClassFunction(k)\n",
       "To check for all groups G4-G37:\n",
       "  Type: IsLengthRationalClassFunctionALL()\n",	
       "Prepend the word Trace to the function name for detailed output.\n\n");

################################################################################
#                      IS LENGTH RATIONAL CLASS FUNCTION
# Input:  STNMBR - a number between 4 and 37 representing a complex
#                  reflection group (as numbered in the Shephard-Todd table)
#
# Output: true  - if absolute reflection length is constant on rational
#                 conjugacy classes
#         false - otherwise
#
# WARNING:  THIS FUNCTION ASSUMES THE IDENTITY IS IN THE FIRST COLUMN
#           OF THE CHARTABLE
#-------------------------------------------------------------------------------
IsLengthRationalClassFunction:=function(STNmbr)

   local refgroupG, t, charNmbr, L;
   
   # Create the exceptional reflection group
   refgroupG:=ComplexReflectionGroup(STNmbr);
   
   # Get the character table of the group
   t:=CharTable(refgroupG);
   
   # Determine position of reflection character in character table
   charNmbr:=GetRefRep(refgroupG,t);
   
   # Compute absolute reflection length
   L:=LengthClassFunction(t,RefClasses(t,charNmbr));
	  
   # Determine if length is constant on rational conjugacy classes
   return IsConstantOnClassUnions(t,RationalConjugacyClasses(t),L);
	  
end;

################################################################################
#                    IS LENGTH RATIONAL CLASS FUNCTION ALL
# This function runs IsLengthRationalClassFunction for all exceptional 
# reflection groups G4-G37.
#-------------------------------------------------------------------------------
IsLengthRationalClassFunctionALL:=function()
   
   local STNmbr;

   for STNmbr in [4..37] do

     Print("G",STNmbr,": ",IsLengthRationalClassFunction(STNmbr),"\n");

   od;

end;

################################################################################
#                    TRACE IS LENGTH RATIONAL CLASS FUNCTION
# Input:  STNMBR - a number between 4 and 37 representing a complex
#                  reflection group in the Shephard-Todd table
#
# Output: This function determines if reflection length is constant on
#         rational conjugacy classes for the exceptional reflection group
#         numbered STNMBR.  The function does not return anything but prints
#         what step it is performing, prints the rational conjugacy classes
#         (as a list of lists of positions of conjugacy classes), prints the
#         values of the reflection length function (same format as list of 
#         rational conjugacy lasses, but conjugacy class positions are replaced
#         with reflection length on that class), and finally, prints whether 
#         or not length is constant on rational conjugacy classes.
#
# WARNING:  THIS FUNCTION ASSUMES THE IDENTITY IS IN THE FIRST COLUMN
#           OF THE CHARTABLE
#-------------------------------------------------------------------------------
TraceIsLengthRationalClassFunction:=function(STNmbr)

   local refgroupG, t, charNmbr, L, ratclasses;
   
   # Create the exceptional reflection group
   refgroupG:=ComplexReflectionGroup(STNmbr);

   # Get the character table of the group
   Print("-----G",STNmbr,"-----\n","Getting character table...\n");
   t:=CharTable(refgroupG);
   
   # Determine position of reflection character in character table
   Print("Getting reflection rep...\n");
   charNmbr:=GetRefRep(refgroupG,t);   
   Print("RefRep: ",charNmbr,"\n");
   
   # Compute absolute reflection length
   Print("Computing reflection length...\n");
   L:=LengthClassFunction(t,RefClasses(t,charNmbr));
	  
   # Get rational conjugacy classes   
   Print("Finding rational conjugacy classes...\n");   
   ratclasses:=RationalConjugacyClasses(t);   
   Print("Rational Conjugacy Classes\n",ratclasses,"\n",
         "Values on Rat. Conj. Classes\n",ValuesOnClassUnions(t,ratclasses,L));

   # Determine if length is constant on rational conjugacy classes
   Print("\n Is Constant on Rat. Conj. Classes? ",
         IsConstantOnClassUnions(t,ratclasses,L),"\n\n");

   return;
	  
end;

################################################################################
#                  TRACE IS LENGTH RATIONAL CLASS FUNCTION ALL
# This function runs TraceIsLengthRationalClassFunction for all exceptional
# reflection groups G4-G37.
#-------------------------------------------------------------------------------
TraceIsLengthRationalClassFunctionALL:=function()
   
   local STNmbr;

   for STNmbr in [4..37] do

     TraceIsLengthRationalClassFunction(STNmbr);

   od;

end;

################################################################################
#EOF
