################################################################################
# spectraprompt.gap
# Briana Foster-Greenwood
# last update: 27 Oct 2015

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
#    gap> Read("/home/username/mygapfiles/spectraprompt.gap");
# 4. A prompt will appear telling you what to type to compute the spectrum of
#    the distance matrix, adjacency matrix, or codimension matrix.
# 5. You can use LogTo if you want to write the results to a file:
#    gap> LogTo("/home/username/mygapfiles/output.txt");
#    And when finished:
#    gap> LogTo();

################################################################################
#                             DEFINE HELPER FILE
helperfile:="/home/username/mygapfiles/spectra.gap";

################################################################################
#                      REQUIRE PACKAGES AND HELPER FILE
RequirePackage("chevie");
Read(helperfile);

################################################################################
#                              STARTUP MESSAGE
Print("\n",
       "To compute spectra for group G_k (for k between 4 and 37):\n",
       "   MatrixSpectrum(k,typestring)\n",
       "To compute spectra for G4-G37:\n",
       "   MatrixSpectrumALL(typestring)\n",
       "To compute spectra for group G(r,p,n):\n",
       "   MatrixSpectrum(r,p,n,typestring)\n",
       "The typestring can be \"distance\", \"codimension\", or \"adjacency\".\n",
       "An eigenvalue lambda occurring with multiplicity m is displayed ",
       "as lambda^{m}.\n\n");	

################################################################################
#                               MATRIX SPECTRUM 
# Input:   r,p,n - parameters for monomial reflection group G(r,p,n)
#          TYPE  - string "adjacency", "codimension", or "distance"
#          FORMAT - string "pairs" or "exponents" specifying the display format
#          k,type,format   OR r,p,n,type,format
#
# Output:  if FORMAT = "pairs":
#             a list of pairs [e,m], where e is an eigenvalue of the 
#             distance/codim/adjacency matrix and m is the multiplicity 
#             of the eigenvalue  
#          if FORMAT = "exponents":
#             a list of e^{m} with eigenvalues e in descending order and 
#             multiplicities m enclosed in ^{} so the output can be copy-pasted 
#             into a LaTeX document
#
# WARNING: THIS FUNCTION ASSUMES THE IDENTITY IS IN THE FIRST COLUMN
#          OF THE CHARTABLE
#
# NOTE:    Group is attached to the character table manually.
#-------------------------------------------------------------------------------
MatrixSpectrum:=function(arg)

   # arg indicates the function may take a variable number of arguments
   # within the function, arg behaves like a list of the arguments
   # MatrixSpectrum:=function(STNmbr,type) or
   # MatrixSpectrum:=function(r,p,n,type)

   local type, refgroupG, t, charNmbr, classFn;
   
   # Create the reflection group
   if Length(arg)=2 then
   
      refgroupG:=ComplexReflectionGroup(arg[1]);
      type:=arg[2];
   
   elif Length(arg)=4 then
      
      refgroupG:=ComplexReflectionGroup(arg[1],arg[2],arg[3]);
      type:=arg[4];
      
   else
     
      Print("invalid number of arguments");
      
      return;
      
   fi;

   # Get the character table

   t:=CharTable(refgroupG);
   Print("-----",t.name,"-----\n");

   # Attach the group to the character table
   t.group:=refgroupG;
   
   # Get position in t.irreducibles of character of reflection representation
   charNmbr:=GetRefRep(t.group,t);
   
   # Get the class function that gives rise to the TYPE of matrix
   if type="distance" then 

      classFn:=LengthClassFunction(t,RefClasses(t,charNmbr));

   elif type="codimension" then

      classFn:=CodimClassFunction(t,charNmbr);

   elif type="adjacency" then

      classFn:=RefClassFunction(t,charNmbr);

   else 

      Print("Invalid spectra type.\n");

      return;

   fi;

   # Print the spectrum as a list of eigenvalue^{multiplicity}
   TeXForm(Spectrum(t,classFn));
	  
end;

################################################################################
#                              MATRIX SPECTRUM ALL
# This function runs MatrixSpectrum(type) for all exceptional reflection
# groups G4-G37.
#-------------------------------------------------------------------------------
MatrixSpectrumALL:=function(type)

  local STNmbr;
  
  for STNmbr in [4..37] do
  
    MatrixSpectrum(STNmbr,type);
  
  od;

end;

################################################################################
#EOF
