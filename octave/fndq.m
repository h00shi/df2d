## reads a .flow file so you can plot the data
## how to use:
## 1- call fndq(<address-of-flow-file)
## 2- the data is stored in the variables: t,qwo,qno,qwin,qnin 

##global variables
global pool t qno qwo qnin qwin vw;

##function cre : creates the variables after a file has been read
function cre
  global pool t qno qwo qnin qwin;

t =    pool(:,2);
qwin = pool(:,5);
qnin = pool(:,3) - qwin;
qwo =  pool(:,6);
qno =  pool(:,4) - qwo;
vw =   pool(:,7);

endfunction

##function readq: reads the data and store them where they should
function readq(adr)
  global pool;
pool = load(adr,"-ascii");
cre;
endfunction
