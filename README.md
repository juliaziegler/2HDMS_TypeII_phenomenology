Code used to produce the plots in https://arxiv.org/abs/2308.05653

what does the code do: 
A bash script executes python scripts along with SPheno and MicrOMEGAs,
to scan over a defined parameter range in the 2HDMS and produces plots such as in the paper above.

how to use the code:
To run the code the correct SPheno and MicrOMEGAs versions have to be installed with the 2HDMS model.
The versions used can be found in the folder 'codes_and_programs'.

The file paths have to be adjusted in the python the bash scripts.
Then the bash script has to be executed (this will do the rest for you and execute all other needed codes).

You can choose to either scan varying two parameters over some range (define these in 'RunMicromegas_NEW.sh') 
or do a random scan (define this in 'RunMicromegas_NEW_random_scan.sh' and 'Basis_Change_NEW_random_scan.py')
