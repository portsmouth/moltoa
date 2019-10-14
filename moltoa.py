
########################################################################################################
#
# Parses molecular structure files in the SDF/MOL format, and converts to Arnold ASS files for rendering 
#
########################################################################################################

import argparse
import os, sys, math, itertools, traceback
import mendeleev

def bondtype_to_desc(bond_type):
    if bond_type is 1: return 'single'
    if bond_type is 2: return 'double'
    if bond_type is 3: return 'triple'
    if bond_type is 4: return 'aromatic'
    if bond_type is 5: return 'single or double'
    if bond_type is 6: return 'single or aromatic'
    if bond_type is 7: return 'double or aromatic'
    if bond_type is 8: return 'any'
    if bond_type is 9: return 'coordination'
    if bond_type is 10: return 'hydrogen'
    return 'unknown'

bondtype_to_strength = {
    1: 1.0,  #single
    2: 2.0,  #double
    3: 3.0,  #triple
    4: 1.5,  #aromatic
    5: 1.5,  #single or double
    6: 1.25, #single or aromatic
    7: 1.75, #double or aromatic
    8: 1,    #any
    9: 1,    #coordination
    10: 0.1  #hydrogen
}

BOND_SCALE = 0.2
ATOM_SCALE = 0.22

def element_to_color(el):
    try:
        if el.symbol.strip().lower()=="h": return (0.5, 0.9, 0.5) # hack to make Hydrogens (slightly) green
        if el.symbol.strip().lower()=="n": return (0.5, 0.5, 0.9) # hack to make Nitrogens blue
        if el.symbol.strip().lower()=="c": return (0.8, 0.8, 0.6) # hack to make Carbons (slightly) yellow
        C_hex = el.cpk_color
        h = C_hex.lstrip('#')
        C_rgb = tuple(float(int(h[i:i+2], 16))/255.0 for i in (0, 2, 4))
        return C_rgb
    except:
        print 'element_to_color failed for element %s' %  el.symbol
        traceback.print_exc()
        return [1.0, 0.5, 0.5]

def element_to_radius(el):
    try:
        R_pm = el.vdw_radius # in pm
        if not R_pm: 
            print 'element %s has no vdw_radius' %  el.symbol
            quit()
        R_ang = R_pm / 100.0    # convert to Angstroms
        return R_ang
    except:
        print 'element_to_radius failed for element %s' %  el.symbol
        traceback.print_exc()
        return 1.0

def emit_atom(element_symbol, position, element_data):
    element = element_data[element_symbol]
    C = element['color']
    R = element['radius']
    return '''

sphere
{{
 center 1 1 POINT {px} {py} {pz}
 radius 1 1 FLOAT {radius_atom}
 shader {shader_atom}
 declare atom_strength constant FLOAT
 atom_strength {strength} 
}}
'''.format(px=str(position[0]),
           py=str(position[1]),
           pz=str(position[2]),
           radius_atom=str(ATOM_SCALE*R),
           strength=0.1,
           shader_atom='shader-atom-%s' % element_symbol)


def lerp(t, a, b):
    return t*b + (1.0-t)*a

def emit_bond(element_symbol_pair, position_pair, bond_type, element_data):
    Ea = element_symbol_pair[0]
    Eb = element_symbol_pair[1]
    Pa = position_pair[0]
    Pb = position_pair[1]
    Ph = ((Pa[0] + Pb[0])/2.0, (Pa[1] + Pb[1])/2.0, (Pa[2] + Pb[2])/2.0)
    element_a = element_data[Ea]
    element_b = element_data[Eb]
    Ca = element_a['color']
    Cb = element_b['color']
    Ra = element_a['radius']
    Rb = element_b['radius']
    bond_radius = BOND_SCALE*min(Ra, Rb)
    num_cvs = 128
    ass = '''
curves
{{
    num_points {num_cvs}
    points {num_cvs} 1 POINT '''.format(num_cvs=num_cvs)
    for n in range(0, num_cvs):
        t = float(n)/float(num_cvs-1)
        ass += '%s %s %s ' % (lerp(t, Pa[0], Pb[0]), lerp(t, Pa[1], Pb[1]), lerp(t, Pa[2], Pb[2]))
    ass += '''
    radius %s 1 FLOAT ''' % str(num_cvs-2)
    for n in range(0, num_cvs-2):
        t = float(n)/float(num_cvs-3)
        ass += '%s ' % str((1.0 - 2.5*t*(1.0-t))*bond_radius)
    ass += '''
    basis "b-spline"
    mode "thick"
    min_pixel_width 0
    receive_shadows off
    self_shadows off
    matrix
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1
    shader {shaderab}
    opaque on
    matte off
}}
'''.format(shaderab='shader-bond-%s-%s-%s' % (bond_type, Ea, Eb),
           shadera='shader-bond-%s-%s' % (bond_type, Ea), 
           shaderb='shader-bond-%s-%s' % (bond_type, Eb))
    return ass


def init_ass(molecule, style):

    if style=='toon':
        filter_ass = '''
contour_filter
{
name filter
width 3.0
}'''
    else:
        filter_ass = '''
gaussian_filter
{
name filter
width 2
}'''

    return '''
options
{{
    AA_samples 3
    outputs "RGBA RGBA filter dr"
    xres 1024
    yres 1024
    GI_transmission_depth 6
    GI_diffuse_depth 6
    GI_specular_depth 6
    GI_transmission_samples 1
    background my_sky
}}

{filter}

user_data_float
{{
 name atom_strength_ud
 attribute atom_strength
 default 1.0
}}

driver_png
{{
    name dr
    filename "{molecule}.png"
}}\n
'''.format(molecule=molecule, filter=filter_ass)


def emit_bond_shaders(style,
                      element_symbol_a, element_a,
                      element_symbol_b, element_b):

    Ca = element_a['color']
    Cb = element_b['color']
    shaders = ''
    for bond_type, strength in bondtype_to_strength.iteritems():
        scatter = 2.0
        shaders += '''
    ramp_rgb
    {{
        name ramp-{name} 
        type "v"
        color 2 1 RGB {ca_r} {ca_g} {ca_b} {cb_r} {cb_g} {cb_b}
        position 2 1 FLOAT 0.0 1.0
        interpolation 2 1 INT 2 2 
    }}'''.format(name='shader-bond-%s-%s-%s' % (bond_type, element_symbol_a, element_symbol_b),
                                ca_r=Ca[0], ca_g=Ca[1], ca_b=Ca[2],
                                cb_r=Cb[0], cb_g=Cb[1], cb_b=Cb[2])

        if style=='toon':
            shader = '''   
    toon
    {{
        name {name}
        base_color ramp-{name}
        edge_color 0 0 0
        emission {em_toon}
        emission_color ramp-{name}
        specular 0.2
        specular_color ramp-{name}
        specular_roughness 0.333
        highlight_size 0.5
        IOR 1.4
        edge_opacity 0.75
        edge_width_scale 0.5
        lights "my_quad_light"
        highlight_color ramp-{name}
        highlight_size 1.0
    }}'''
        else:
            shader = '''   
    standard_surface
    {{
        name {name}
        base {strength}
        base_color ramp-{name} 
        emission {em_standard}
        emission_color ramp-{name}
        specular 1.0
        specular_IOR 1.7
        specular_roughness 0.3
        metalness 0.333
    }}
    '''
        shaders += shader.format(name='shader-bond-%s-%s-%s' % (bond_type, element_symbol_a, element_symbol_b),
                                ca_r=Ca[0], ca_g=Ca[1], ca_b=Ca[2],
                                cb_r=Cb[0], cb_g=Cb[1], cb_b=Cb[2],
                                strength=0.333*strength,
                                em_toon=0.75*(strength-1.0),
                                em_standard=0.333*(strength-1.0))
    return shaders

def emit_atom_shader(style, element_symbol, element):
    C = element['color']
    R = element['radius']
    if style=='toon':
        shader = '''
    toon
    {{
        name shader-atom-{element}
        enable 1
        base_color {cx} {cy} {cz}
        edge_color {cx} {cy} {cz}
        specular 0.2
        specular_color {cx} {cy} {cz}
        specular_roughness 0.666
        highlight_size 0.5
        IOR 1.4
        edge_opacity 0.5
        edge_width_scale 0.75
        lights "my_quad_light"
        highlight_color 1 1 1
        highlight_size 1.0
    }}'''
    else:
        shader = '''
    standard_surface
    {{
        name shader-atom-{element}
        base 0.333
        base_color {cx} {cy} {cz}
        specular 1.0
        specular_IOR 1.7
        specular_color {cx} {cy} {cz}
        specular_roughness 0.3
        metalness 0.5
    }}'''
    return shader.format(element=element_symbol, cx=C[0], cy=C[1], cz=C[2])


def drawProgressBar(percent, barLen=40):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.2f%%" % (progress, percent * 100))
    sys.stdout.flush()

####################################################
# parse v3000
####################################################

def v3000_parser(args, data, line,
                 atoms, bonds):

    print 'Parsing V3000 SDF file ...'
    try:
        counts = data[line].split('COUNTS')[1].split()
        num_atoms = int(counts[0])
        num_bonds = int(counts[1])
    except Exception as e:
        print 'Error parsing SDF counts: %s' % str(e)
        traceback.print_exc()
        quit()

    if args.verbosity:
        print 'num atoms: ', num_atoms
        print 'num bonds: ', num_bonds
        
    line += 1

    print 'Parsing atom block ...'
    try:
        while True:
            if 'BEGIN ATOM' in data[line]: break
            line += 1
            continue
        line += 1
        ia = 0
        while True:
            if 'END ATOM' in data[line]: break
            atom = data[line].split('V30')[1].split()
            atom_index   = int(atom[0])
            atom_element = atom[1]
            atom_x       = float(atom[2])
            atom_y       = float(atom[3])
            atom_z       = float(atom[4])
            atoms.append({'index':atom_index,
                          'element':atom_element,
                          'x':atom_x,
                          'y':atom_z,
                          'z':atom_y})
            if args.verbosity:
                print '\natom %d: %s (%s %s %s)' % (atom_index, atom_element, atom_x, atom_y, atom_z)
            line += 1
            drawProgressBar(ia/float(num_atoms))
            ia += 1
    except Exception as e:
        print 'Error parsing SDF atom block, line %d: %s' % (line, str(e))
        traceback.print_exc()
        quit()

    line += 1

    print '\nParsing bond block ...'
    try:
        while True:
            if 'BEGIN BOND' in data[line]: break
            line += 1
            continue
        line += 1
        ib = 0
        while True:
            if 'END BOND' in data[line]: break
            bond = data[line].split('V30')[1].split()
            bond_index = int(bond[0])
            bond_type  = int(bond[1])
            bond_atom1 = int(bond[2])-1 # subtract 1 since indices are 1-based
            bond_atom2 = int(bond[3])-1 # subtract 1 since indices are 1-based
            bonds.append({'index':bond_index,
                          'type':bond_type,
                          'atom1':bond_atom1,
                          'atom2':bond_atom2})
            if args.verbosity:
             print '\nbond %d: %s (%s-%s)' % (bond_index, bondtype_to_desc(bond_type), bond_atom1, bond_atom2)
            line += 1
            drawProgressBar(ib/float(num_bonds))
            ib += 1
    except Exception as e:
        print 'Error parsing SDF bond block, line %d: %s' % (line, str(e))
        traceback.print_exc()
        quit()


####################################################
# parse v2000
####################################################

def v2000_parser(args, data, line,
                 atoms, bonds):

    print 'Parsing V2000 SDF file ...'
    try:
        counts = data[line].split()
        num_atoms = int(counts[0])
        num_bonds = int(counts[1])
    except Exception as e:
        print 'Error parsing SDF counts: %s' % str(e)
        traceback.print_exc()
        quit()

    if args.verbosity:
        print 'num atoms: ', num_atoms
        print 'num bonds: ', num_bonds
        
    line += 1

    print 'Parsing atom block ...'
    try:
        ia = 0
        while ia < num_atoms:
            atom = data[line].split()
            atom_x       = float(atom[0])
            atom_y       = float(atom[1])
            atom_z       = float(atom[2])
            atom_element = atom[3]
            atoms.append({'index':ia,
                          'element':atom_element,
                          'x':atom_x,
                          'y':atom_z,
                          'z':atom_y})
            if args.verbosity:
                print '\natom %d: %s (%s %s %s)' % (ia, atom_element, atom_x, atom_y, atom_z)
            line += 1
            drawProgressBar(ia/float(num_atoms))
            ia += 1
    except Exception as e:
        print 'Error parsing SDF atom block, line %d: %s' % (line, str(e))
        traceback.print_exc()
        quit()

    print '\nParsing bond block ...'
    try:
        ib = 0
        while ib < num_bonds:
            bond = data[line].split()
            bond_atom1 = int(bond[0])-1 # subtract 1 since indices are 1-based
            bond_atom2 = int(bond[1])-1 # subtract 1 since indices are 1-based
            bond_type  = int(bond[2])
            bonds.append({'index':ib,
                          'type':bond_type,
                          'atom1':bond_atom1,
                          'atom2':bond_atom2})
            if args.verbosity:
             print '\nbond %d: %s (%s-%s)' % (ib, bondtype_to_desc(bond_type), bond_atom1, bond_atom2)
            line += 1
            drawProgressBar(ib/float(num_bonds))
            ib += 1
    except Exception as e:
        print 'Error parsing SDF bond block, line %d: %s' % (line, str(e))
        traceback.print_exc()
        quit()


####################################################
# main
####################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert molecular SDF (V2000, V3000) files to Arnold ASS.')
    parser.add_argument('file')
    parser.add_argument('--style', nargs='?', default='standard', choices=['standard', 'toon'])
    parser.add_argument('--verbosity', help="increase output verbosity")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print 'No file %s exists' % args.file
        quit()
    if not args.file.endswith(".sdf"):
        print 'expecting a .sdf file'
        quit()

    style = args.style

    print 'Reading SDF data from %s' % args.file

    ########################################
    # Parse SDF file
    ########################################

    print 'Opening SDF file ...'
    with open(args.file, 'r') as f:
        data = f.readlines()

    atoms = []
    bonds = []
    try:
        line = 0
        while True:
            if 'BEGIN CTAB' in data[line]: 
                v3000_parser(args, data, line+1, atoms, bonds)
                break
            if 'V2000' in data[line]: 
                v2000_parser(args, data, line, atoms, bonds)
                break
            line += 1    
    except Exception as e:
        print 'Error parsing SDF header: %s' % str(e)
        traceback.print_exc()
        quit()

    print '\nCaching element data ...'
    element_data = {}
    try:
        element_set = {}
        for atom in atoms:
            element_symbol = atom['element']
            element_set[element_symbol] = None
        num_elements = len(element_set)
        ie = 0
        for element_symbol, element in element_set.iteritems():
            drawProgressBar(ie/float(num_elements))
            el = mendeleev.element(element_symbol)
            R = element_to_radius(el)
            C = element_to_color(el)
            element_symbol = el.symbol.lower()
            element_data[element_symbol] = {'radius':R, 'color':C}
            ie += 1
    except Exception as e:
        print 'Error caching element data: %s' % str(e)
        traceback.print_exc()
        quit()

    # Export .ass file
    molecule = os.path.splitext(args.file)[0]
    ass = init_ass(molecule, style)

    centroid = []
    Pmax = []
    Pmin = []

    try:
        print '\nExporting atoms to ASS file %s.ass' % molecule
        num_atoms = len(atoms)
        for atom_index in range(0, num_atoms):
            drawProgressBar(atom_index/float(num_atoms))
            atom = atoms[atom_index]
            el = atom['element'].lower()
            P = (atom['x'], atom['y'], atom['z'])
            if len(centroid)==0: 
                centroid = P
            else:
                centroid = [centroid[0]+P[0], centroid[1]+P[1], centroid[2]+P[2]]
            if len(Pmax)==0: 
                Pmax = P
            else:
                Pmax = [max(Pmax[0], P[0]), max(Pmax[1], P[1]), max(Pmax[2], P[2])]
            if len(Pmin)==0: 
                Pmin = P
            else:    
                Pmin = [min(Pmin[0], P[0]), min(Pmin[1], P[1]), min(Pmin[2], P[2])]
            ass += emit_atom(el, P, element_data)
        for i in [0,1,2]: centroid[i] /= num_atoms
    except Exception as e:
        print '\nError exporting atoms: %s' % str(e)
        traceback.print_exc()
        quit()

    try:
        print '\nExporting bonds to ASS file %s.ass' % molecule
        num_bonds = len(bonds)
        for bond_index in range(0, num_bonds):
            drawProgressBar(bond_index/float(num_bonds))
            bond = bonds[bond_index]
            bond_type = bond['type']
            atom1_index = int(bond['atom1'])
            atom2_index = int(bond['atom2'])
            if atom1_index<0 or atom1_index>=len(atoms):
                print 'Out-of-range atom index %d in bond %d: %s' % (atom1_index, bond_index, str(bond))
                quit()
            if atom2_index<0 or atom2_index>=len(atoms):
                print 'Out-of-range atom index %d in bond %d: %s' % (atom2_index, bond_index, str(bond))
                quit()
            atom1 = atoms[atom1_index]
            atom2 = atoms[atom2_index]
            el1 = atom1['element'].lower()
            el2 = atom2['element'].lower()
            P1 = (atom1['x'], atom1['y'], atom1['z'])
            P2 = (atom2['x'], atom2['y'], atom2['z'])
            ass += emit_bond((el1, el2), (P1, P2), bond_type, element_data)
    except Exception as e:
        print '\nError exporting bonds: %s' % str(e)
        traceback.print_exc()
        quit()

    try:
        print '\nGenerating atom shaders ...'
        n = 0
        num_elements = len(element_data)
        for element_symbol, element in element_data.iteritems():
            drawProgressBar(n/float(num_elements))
            ass += emit_atom_shader(style, element_symbol, element)
            n += 1
        print '\nGenerating bond shaders ...'
        n = 0
        for element_symbol, element in element_data.iteritems():
            drawProgressBar(n/float(num_elements*num_elements))
            ass += emit_bond_shaders(style, element_symbol, element, element_symbol, element)
            n += 1
        for element_symbol_a, element_a in element_data.iteritems():
            for element_symbol_b, element_b in element_data.iteritems():
                if element_symbol_a == element_symbol_b: continue
                drawProgressBar(n/float(num_elements*num_elements))
                ass += emit_bond_shaders(style, element_symbol_a, element_a, element_symbol_b, element_b)
                n += 1
    except Exception as e:
        print '\nError generating shaders: %s' % str(e)
        traceback.print_exc()
        quit()

    # Generate camera
    extent = (Pmax[0]-Pmin[0], Pmax[1]-Pmin[1], Pmax[2]-Pmin[2])
    max_extent = max(extent[0], extent[1], extent[2])
    cam_trg = centroid
    cam_pos = (centroid[0], centroid[1] - 0.75*max(extent[0], extent[1], extent[2]), centroid[2])
    cam_dist = 0.0
    for i in (0,1,2): cam_dist += (cam_pos[i]-cam_trg[i])**2.0
    cam_dist = cam_dist**0.5
    if style != 'toon': aperture = 0.005*cam_dist
    else:               aperture = 0.0
    ass += '''
persp_camera
{{
 position {px} {py} {pz} 
 look_at {tx} {ty} {tz}
 up 0 0 1
 fov 60
 focus_distance {fd}
 aperture_size {ap}
}}
'''.format(px=cam_pos[0], py=cam_pos[1], pz=cam_pos[2],
           tx=cam_trg[0], ty=cam_trg[1], tz=cam_trg[2],
           fd=cam_dist, ap=aperture)

    # Basic lighting
    ass += '''
physical_sky
{{
 name my_sky
 enable_sun on
 sun_size 5
 elevation 80
 azimuth 30
 intensity 1
 turbidity 3
 enable_sun false
}}

skydome_light
{{
 name mylight
 resolution 1000
 color my_sky
 samples 2
}}

quad_light
{{
    name my_quad_light
    vertices 4 1 POINT
    {x1} {y1} {z1}
    {x2} {y2} {z2}
    {x3} {y3} {z3}
    {x4} {y4} {z4}
    intensity 10
    samples 3
    normalize off
}}
'''.format( x1=centroid[0]+2.0*extent[0], y1=centroid[1]-max(extent[0], extent[1], extent[2]), z1=centroid[2]-extent[2],
            x2=centroid[0]+2.0*extent[0], y2=centroid[1]-max(extent[0], extent[1], extent[2]), z2=centroid[2]+extent[2],
            x3=centroid[0]+2.0*extent[0], y3=centroid[1]+max(extent[0], extent[1], extent[2]), z3=centroid[2]+extent[2],
            x4=centroid[0]+2.0*extent[0], y4=centroid[1]+max(extent[0], extent[1], extent[2]), z4=centroid[2]-extent[2])

    # Output file
    open("%s.ass"%molecule, 'w').write(ass)

    print '\nDone.' 
