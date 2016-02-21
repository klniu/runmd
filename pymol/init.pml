select drg, resn DRG
select oil, resn OIL
select sol, resn SOL

# turn on the representation for drg in PyMOL
#set_bond stick_color, black, drg, drg
set_bond stick_radius, 0.14, drg, drg
set sphere_scale, 0.25, drg
show sticks, drg
show spheres, drg
