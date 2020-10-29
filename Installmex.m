%%***********************************************************************
%% compile mex files
%%
%%***********************************************************************

   function Installmex

   warning off
   
   cdirectory = pwd;  
   fprintf(' current directory is:  %s\n',cdirectory);
   
   src =[cdirectory,filesep,'mexfun'];
   eval(['cd ', 'mexfun']);
   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   fprintf (' %s\n',src); 
   
   mex Amap.cpp
   
   mex Atmap.cpp
   
   mex BCD_Zpart.cpp
   
   mex dinfpart.cpp
   
   mex Fmap.cpp
   
   mex simplex_y.cpp
   
   mex  COMPFLAGS="$COMPFLAGS -openmp" Simplex_matcol.cpp