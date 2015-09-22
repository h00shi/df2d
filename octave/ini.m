#############################################################################
##                                    ini.m                                ##
#############################################################################


## file ini.m 
## helps visualize the results of DF_1d in GNU octave. 
## written by Shayan Hoshyari
## January 2015


##global variables
global q x s p pc;

##function cre : creates the variables after a file has been read
function cre
	global q x s p pc;

x = q(:,4);
s = q(:,1);
p = q(:,2);
pc = q(:,3);

endfunction

##function lres: reads a restart file
function lres(n,dir)
	global q x s p pc;

	if (nargin == 0)
		printf("give me input man.\n");
		return;
	elseif (nargin == 1)
	dir = ".";
	endif
	
q = load([dir,"/restart/restart.",mat2str(n)],"-ascii");
cre;

endfunction

##function lini: reads an initialize file
function lini(dir)
	global q x s p pc;

	if (nargin == 0)
		dir = ".";
	endif
	
q = load([dir,"/initial/initial"],"-ascii");
cre;

endfunction

