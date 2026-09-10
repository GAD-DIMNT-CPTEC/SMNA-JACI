-------------------------------------------------------------------
 check data: scripts for checking observed data routinely received
-------------------------------------------------------------------

To use the verification scripts, it is necessary to previously compile and install MBUFRTOOLS

 __________________________________________

 1- Compilation and Intallation  MBUFRTOOLS
 ___________________________________________
   1.1 - Compilation in Egeon using gfortran 
  
     ln -s makefile_gfortran  makefile.comp
     
     make 
     
   1.2  - setting environment variables
     Include in .bashrc or in scripts the followed command 
  
     export MBUFR_TABLES=path_to_where_is_bufrtables
     
     example:    
     export MBUFR_TABLES=$HOME/smna/branch/SMNA_v3.0.0.t11889/SMG/utils/mbufrtools/bufrtables
     
   1.3 - check if it is correctly installed
      
      The following command should correctly decode the example file  
      
      ./bin/bufrdump -i ./data_example/BUOY_ARGOS.bufr -o teste.txt
     
 
  ________________________
  
  2 - check_data_scripts
  _______________________
 
       
       Currently, only one data verification script is available in check_data directory
       cd check_data
       
       If necessay edit the dirconf.sh to specify the input directory and than      
       use simple_check.sh to verify the contant of avalible files in an input_dir

       simple_check.sh  to see options or 
       simple_check.sh <date> <sinoptic time> 
       
       Examples 
       
       simple_check.sh yesterday 00 <- Check the contents of of all PREPBUFR file received yesterday relative to 00 UTC
       simple_check.sh yesterday 12 <- Check the contents of of all PREPBUFR file received yesterday relative to 12 UTC
       simple_check.sh yesterday    <- Check the contents of of all PREPBUFR file received yesterday all times
       simple_check.sh now          <- Check the contents of of all PREPBUFR file received today
       simple_check.sh 20240701 00  <- Check the contents of of all PREPBUFR file received on 2024-07-01 at 00 UTC
       
       It is also possible to redirect the output to a file
       
       simple_check.sh 20240701 00 > result_2024070100.txt 
       
       
     
       
       
       
       
        
