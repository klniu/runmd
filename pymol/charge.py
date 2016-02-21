from pymol import cmd
cmd.bg_color('white')
cmd.select('drg', 'all')

cmd.set('valence', '1')
cmd.set_bond('stick_radius', '0.14', 'drg', 'drg')
cmd.set('sphere_scale', '0.25', 'drg')
cmd.show('sticks', 'drg')
cmd.show('spheres', 'drg')

cmd.set_view([0.984318674,   -0.042038180,   -0.171310961,\
     0.001240913,   -0.969511747,    0.245043248,\
    -0.176387057,   -0.241413504,   -0.954258323,\
    -0.000004895,    0.000000954,  -62.691986084,\
     1.964830756,   -0.868917644,   -0.235743612,\
    37.651191711,   87.732429504,  -20.000000000])

charges = ["0.152", "0.453", "-0.224", "0.343", "0.269", "0.472", "-0.447", "-0.416", "0.161", "0.203", "-0.055", "0.060", "0.008", "-0.007", "-0.015", "0.051", "0.049", "-0.218", "0.137", "-0.240", "0.140", "0.066", "0.072", "-0.062", "-0.006", "0.131", "-0.101", "0.152", "-0.051", "0.410", "-0.623", "0.422", "-0.574", "0.405", "-0.342", "0.168", "-0.031", "0.126", "-0.290", "0.158", "-0.190", "0.950", "-0.522", "-0.572", "-0.572"]
cmd.set('label_size', '22')
#cmd.set('label_position', (2, 0, 0))
for i, charge in zip(range(1, 46), charges):
    cmd.label('id %s' % i, charge)

# output
cmd.set('antialias', '2')
cmd.set('direct', '0.6')
cmd.set('orthoscopic', 'on')
cmd.set('ray_trace_frames', '1')
#cmd.ray(renderer=-1)
#cmd.png('~/Downloads/charge.png', dpi=600)

