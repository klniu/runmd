# make_pov.py
# Do "run make_pov.py" from within pymol and then execute the script
# with "make_pov('povray.inp')" to create the povray.inp file.
#
# written by Robert Campbell 2003
#
from pymol import cmd
 
def make_pov(file):
	(header,data) = cmd.get_povray()
	povfile=open(file,'w')
	povfile.write(header)
	povfile.write(data)
	povfile.close()
