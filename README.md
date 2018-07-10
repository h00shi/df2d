		This project hosted at git-hub holds the code written as my
		BASC project.


	 The code is documented using doxygen. If you want to generate it
	 again use:

	 $doxygen doxy.config

	 To get to know the code it would be a good idea to read the
	 documentation.  In the doxygen mainpage, I have explained how to
	 compile the source code and use the df2d executable file.

	 The folder df2d/code/examples contains input files that df2d
	 needs to run the three benchmark problems that I have solved in
	 my thesis. The bench3 directory is commented, but the others do
	 not have comments, as all the folders are almost similar. If you
	 want to obtain results you have to run df2d in that folder. For
	 example typing:

	 $ make bench2

	 (if df2d exefile is located in df2d/code/bin) or:

	 $ /path/to/df2d -s $ /path/to/df2d

	 In that folder will start df2d, and the results will be put in
	 result folder. I have included the results for bench3, but to
	 conserve space the results of other bench folders are not
	 included. You can run df2d on them to obtain results.

	 
	 Shayan Hoshyari 
	 June 2015 

