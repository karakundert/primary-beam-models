from numpy import array

pylabcolourstring = [
        "black", "blue", "green", "cyan", "red", "yellow",
        "darkorange","aqua", "blueviolet", "brown", "burlywood",
        "chartreuse", "chocolate", "coral", "cornflowerblue",
        "crimson", "darkblue", "darkcyan", "darkgoldenrod",
        "darkgray", "darkgreen", "darkkhaki", "darkmagenta",
        "darkolivegreen", "darkorange", "darkorchid", "darkred",
        "darksalmon", "darkseagreen","darkslateblue", "darkslategray",
        "darkturquoise", "darkviolet", "deeppink", "deepskyblue",
        "dimgray", "dodgerblue", "firebrick", "forestgreen",
        "fuchsia", "gainsboro", "gold", "goldenrod", "gray",
        "greenyellow", "hotpink", "indianred", "indigo",
        "khaki", "lavender", "lawngreen", "lightblue", "lightcoral",
        "lightgoldenrodyellow", "lightgreen", "lightgrey", "lightpink",
        "lightsalmon", "lightseagreen", "lightskyblue", "lightslategray",
        "lightsteelblue", "lime", "limegreen", "maroon", "mediumaquamarine",
        "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen",
        "mediumslateblue", "mediumspringgreen", "mediumturquoise",
        "mediumvioletred", "midnightblue", "moccasin", "navy", "olive",
        "olivedrab", "orange", "orchid", "orangered", "palegoldenrod",
        "palegreen", "palevioletred", "peru", "pink", "plum",
        "powderblue", "purple", "rosybrown", "royalblue", "saddlebrown",
        "salmon", "sandybrown", "seagreen", "sienna", "silver",
        "skyblue", "slateblue", "slategray", "springgreen", "cornsilk",
        "steelblue", "tan", "teal", "thistle", "tomato", "turquoise", "violet",
        "wheat", "yellowgreen", "aquamarine", "azure", "bisque", "cadetblue"]

cor_vals_single = []
cor_vals_multi  = []
rmsValues_single = []
rmsValues_multi = []

master_dirlist = [ '/lustre/kkundert/Code/results/7m_aper_test/rot_aper/Data001-',
                   '/lustre/kkundert/Code/results/results/handmade_tests/plus/no_perturbation/Data000-',
                   '/lustre/kkundert/Code/results/handmade_tests/plus_PA/no_perturbation/Data000-',
                   '/lustre/kkundert/Code/results/handmade_tests/plus_cross/rot_aper/Data001-',
                   '/lustre/kkundert/Code/results/handmade_tests/plus_cross_PA/rot_aper/Data001-',
                   '/lustre/kkundert/Code/results/handmade_tests/plus_pointing_uncor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/handmade_tests/plus_pointing_cor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/handmade_tests/plus_pointing_cor_PA/pointing/Data000-',
                   '/lustre/kkundert/Code/results/handmade_tests/plus_cross_pointing_cor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus/Data000-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus_PA/no_perturbation/Data000-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus_cross/rot_aper/Data001-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus_cross_PA/rot_aper/Data001-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus_pointing_uncor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus_pointing_cor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus_pointing_cor_PA/pointing/Data000-',
                   '/lustre/kkundert/Code/results/ray_trace_tests/plus_cross_pointing_cor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus/no_perturbation/Data000-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus_PA/no_perturbation/Data000-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus_cross/rot_aper/Data001-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus_cross_PA/rot_aper/Data001-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus_pointing_uncor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus_pointing_cor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus_pointing_cor_PA/pointing/Data000-',
                   '/lustre/kkundert/Code/results/stuart_beam_tests/plus_cross_pointing_cor/pointing/Data000-',
                   '/lustre/kkundert/Code/results/illumination_offset_tests/no_perturbation/Data000-']


# choose the comparable data set
numRegex = r"[+-]?([0-9]+\.)?[0-9]+([eE][+-][0-9]+)?"
requiredWhiteSpaceRegex = r"[ \t]+"
optionalWhiteSpaceRegex = r"[ \t]*"
rmsRegex = r"(?P<rms>%s)" % numRegex
newline = r"\n"
headerRegex = (
        optionalWhiteSpaceRegex
        + "#"
        + optionalWhiteSpaceRegex
        + "Npts"
        + requiredWhiteSpaceRegex
        + "Sum"
        + requiredWhiteSpaceRegex
        + "Mean"
        + requiredWhiteSpaceRegex
        + "Rms"
        + requiredWhiteSpaceRegex
        + "Std dev"
        + requiredWhiteSpaceRegex
        + "Minimum"
        + requiredWhiteSpaceRegex
        + "Maximum"
        + optionalWhiteSpaceRegex
        )
regex = (
        headerRegex
        + newline
        + optionalWhiteSpaceRegex
        + 3 * (numRegex + requiredWhiteSpaceRegex)
        + rmsRegex
        + 3 * (requiredWhiteSpaceRegex + numRegex)
        + optionalWhiteSpaceRegex
        )
pattern = re.compile(regex, re.IGNORECASE)
 
ia.open('/lustre/kkundert/Code/handmade_tests/plus/no_perturbation/Data000-single/simple.flux')
pb_vals = ia.getchunk()
ia.close()

##########################################################

def find_cor_vals(pathname = 'Data000-single/', pb = []):

    pathname_single = pathname+'single/'
    pathname_multi = pathname+'multi/'

    single_coords = [[493,532]]
    multi_coords = [[465,536],[525,476],[573,560]]

    i = 0
    
    while i < 2:
        cor_vals = []

        if i == 0:
            pathname = pathname_single
            coord_list = single_coords
        else:
            pathname = pathname_multi
            coord_list = multi_coords

        ia.open(pathname+'simple.image')
        pix = ia.getchunk()
        ia.close()

        for coords in coord_list:
            cor_vals.append(pix[coords[0],coords[1],0,0] / pb[coords[0],coords[1],0,0])

        print cor_vals

        outfile = pathname+'corrected_sources.txt'

        f = open(outfile,'w')
        f.write(str(cor_vals))
        f.close()

        i = i + 1

##########################################################

def residual_statistics(pathname = 'Data000-single/'):
    
    pathname_single = pathname+'single/'
    pathname_multi = pathname+'multi/'

    i = 0

    while i < 2:

        if i == 0:
            pathname = pathname_single
        else:
            pathname = pathname_multi

        filename = pathname+'simple.residual'

        ia.open(filename)
        ia.statistics(logfile=pathname+'residuals.txt')
        ia.close()

        i = i + 1

##########################################################

def read_cor_file(filename = 'corrected_sources.txt'):

    f = open(filename, 'r')
    cor_vals = eval(f.read())
    return cor_vals

##########################################################

for dirname in master_dirlist:

    pathname_single = dirname+'single/'
    pathname_multi = dirname+'multi/'

    find_cor_vals(pathname = dirname, pb = pb_vals)
    residual_statistics(pathname = dirname)

    n = 0

    while n < 2:

        if n == 0:
            filename = pathname_single+'corrected_sources.txt'
            cor_vals_single.append(read_cor_file(filename = filename))
        else:
            filename = pathname_multi+'corrected_sources.txt'
            cor_vals_multi.append(read_cor_file(filename = filename))

        n = n + 1
    
    # generates plots for all statistical locations
    with open(pathname_single+"/residuals.txt") as f:
        contents = f.read()
        match = pattern.search(contents)
        assert match
        rms_val = float(match.group('rms'))
        rmsValues_single.append(rms_val)

    with open(pathname_multi+"/residuals.txt") as f:
        contents = f.read()
        match = pattern.search(contents)
        assert match
        rms_val = float(match.group('rms'))
        rmsValues_multi.append(rms_val)

i = 0
j = 0
pl.figure(1)
pl.clf()

for val in rmsValues_single:
    left = j
    base = -3
    pl.bar(left=left, height=log10(val)-base, width=0.05,
            bottom = base, color = pylabcolourstring[i])
    j += 0.05
    i += 1

pl.suptitle('RMS Levels of All Effects - Single Source', fontsize=24)
pl.savefig('/users/kkundert/winter_talk/all_effects_rms_single')
pl.show()

i = 0
j = 0
pl.figure(2)
pl.clf()

for val in rmsValues_multi:
    left = j
    base = -3
    pl.bar(left=left, height=log10(val)-base, width=0.05,
            bottom = base, color = pylabcolourstring[i])
    j += 0.05
    i += 1

pl.suptitle('RMS Levels of All Effects - Multi Source', fontsize=24)
pl.savefig('/users/kkundert/winter_talk/all_effects_rms_multi')
pl.show()

i = 0
j = 0
pl.figure(3)
pl.clf()

for corr in cor_vals_single:
    left = j
    base = 0.78
    pl.bar(left=left, height=corr[0]-base, width=0.05,
            bottom = base, color = pylabcolourstring[i])
    j += 0.05
    i += 1

#pl.legend()
pl.suptitle('Image Fidelity of All Effects - Single Source', fontsize=24)
pl.savefig('/users/kkundert/winter_talk/all_effects_fidelity_single')
pl.show()

i = 0
j = 0
k = 1.5
l = 3.0
pl.figure(4)
pl.clf()

for corr in cor_vals_multi:
    left0 = j
    left1 = k
    left2 = l
    base = 0.4
    pl.bar(left=left0, height=corr[0]-base, width=0.05,
            bottom = base, color = pylabcolourstring[i])
    pl.bar(left=left1, height=corr[1]-base, width=0.05,
            bottom = base, color = pylabcolourstring[i])
    pl.bar(left=left2, height=corr[2]-base, width=0.05,
            bottom = base, color = pylabcolourstring[i])
    j += 0.05
    k += 0.05
    l += 0.05
    i += 1

#pl.legend()
pl.suptitle('Image Fidelity of All Effects - Multi Source', fontsize=24)
pl.savefig('/users/kkundert/winter_talk/all_effects_fidelity_multi')
pl.show()
