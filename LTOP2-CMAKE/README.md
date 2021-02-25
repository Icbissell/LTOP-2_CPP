# LTOP2-CMAKE

This is an example of the LTOP-2 program, modelling how precipitation would fall based on the LTOP-2 model. For the gaussian demonstration, a simple gaussian hill is written directly into an array. For the Himalayan demonstration, the script reads from a DEM file (hGrid3.csv), writes this to an array, and then processes. A precipitation and topography plot are generated and saved to the build file. The maximum precipitation in (m/a) and (mm/hr) are also outputted to the user's screen. 

In order to run the complete program, you must have cmake and gnuplot installed, which can be installed on Linux with:
	
	sudo apt-get install cmake 
	sudo apt-get install gnuplot 

To run this script, simply 

1) Save the folder to your computer
2) Navigate to the folder through your command prompt
3) Enter 

	
		cmake CMakeLists.txt 
		make all
		
4) And finally, run the program with 

		./ltop-2
		

If one encounters an error message along the lines of `bootstrap.sh: Permission denied` after `make all` in step 3, run 
		
		find . -name "*.sh" | xargs chmod a+x		
