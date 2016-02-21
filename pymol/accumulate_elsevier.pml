reset
bg_color white
select drg, resn DRG
select oil, resn OIL
select sol, resn SOL
select ball_atoms, name N+O+S+H
select other_atoms, drg and not ball_atoms

#orient
turn x, 90
zoom center, 110

# hide the oil
hide everything, oil
show line, oil
set_color vdgray=[0.99, 0.99, 0.99]
color vdgray, oil

# hide the hydrogen in water
hide everything, name HW1+HW2
select ow, name OW
set sphere_scale, 0.1, ow
show spheres, ow
alter ow, vdw=1.2
rebuild

hide everything, drg

# turn on the stick-ball representation for ball_atoms in drg in PyMOL
set_bond stick_radius, 0.14, drg, drg
set sphere_scale, 0.25, drg
show sticks, drg
show spheres, drg
#set_bond stick_radius, 0.14, ball_atoms, ball_atoms
#set sphere_scale, 0.25, ball_atoms
#show sticks, ball_atoms
#show spheres, ball_atoms
#
## turn on other atoms in drg to lines
#show lines, other_atoms

# change the color of some atoms
color purple, name NA

# output
set antialias, 3
set direct, 0.6
set orthoscopic, on
ray 720, 2000, renderer=1
#run /home/klniu/backup/python/md/pymol/make_pov.py
#png /home/klniu/Downloads/test.png, dpi=600
