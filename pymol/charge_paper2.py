from pymol import cmd
cmd.bg_color('white')
cmd.select('drg', 'all')

cmd.set('valence', '1')
cmd.set_bond('stick_radius', '0.14', 'drg', 'drg')
cmd.set('sphere_scale', '0.25', 'drg')
cmd.show('sticks', 'drg')
cmd.show('spheres', 'drg')
cmd.set_view ([0.987891138,   -0.139472052,   -0.067891687,
  0.152521998,    0.793636620,    0.588958621,
 -0.028259384,   -0.592185259,    0.805302858,
  0.000017954,    0.000006792,  -52.386489868,
 -1.638074398,   -1.409737468,   -0.143483341,
-34.060420990,  138.833740234,   20.000000000])

charges = ["-0.103" ," 0.115" ," 0.016" ,"-0.082" ," 0.068" ,"-0.022" ,"-0.017" ," 0.017" ,"-0.001" ,"-0.051" ," 0.120" ,"-0.186" ," 0.163" ,"-0.115" ," 0.216" ,"-0.142" ," 0.130" ,"-0.294" ," 0.114" ,"-0.004" ," 0.102" ,"-0.098" ,"-0.004" ," 0.102" ,"-0.294" ," 0.114" ," 0.882" ,"-0.582" ,"-0.582" ,"-0.582"]
cmd.set('label_size', '22')
cmd.set('label_position', (0, 2, 2))
for i, charge in zip(range(1, 31), charges):
    cmd.label('id %s' % i, charge)

# output
cmd.set('antialias', '2')
cmd.set('direct', '0.6')
cmd.set('orthoscopic', 'on')
cmd.set('ray_trace_frames', '1')
#cmd.ray(renderer=-1)
#cmd.png('~/Downloads/charge.png', dpi=600)

