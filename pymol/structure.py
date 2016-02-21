from pymol import cmd
cmd.bg_color('white')
cmd.select('drg', 'all')

cmd.set('valence', 1)
cmd.set_bond('stick_radius', '0.14', 'drg', 'drg')
cmd.set('sphere_scale', '0.25', 'drg')
cmd.show('sticks', 'drg')
cmd.show('spheres', 'drg')

# output
cmd.set('antialias', '2')
cmd.set('direct', '0.6')
cmd.set('orthoscopic', 'on')
cmd.set('ray_trace_frames', '1')
#cmd.ray(renderer=-1)
#cmd.png('~/Downloads/charge.png', dpi=600)

