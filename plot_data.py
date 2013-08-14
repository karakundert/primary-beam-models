from glob import glob
from matplotlib import axes
# choose the comparable data set

import re
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


for level in ['centered', 'half-power,' 'low-power']:
    print "\n###", level.upper()
    filenames = glob('eccentricity/Data*-%s' % level)
    #? xvals = xrange(len(filenames))
    rmsValues = []
    rmsValuesNear = []
    rmsValuesOff = []

    for filename in filenames:
        with open(filename+"/all_stats.txt") as f:
            contents = f.read()
            match = pattern.search(contents)
            assert match
            rmsValues.append(float(match.group('rms')))
        with open(filename+"/near_source_stats.txt") as f:
            contents = f.read()
            match = pattern.search(contents)
            assert match
            rmsValuesNear.append(float(match.group('rms')))
        with open(filename+"/off_source_stats.txt") as f:
            contents = f.read()
            match = pattern.search(contents)
            assert match
            rmsValuesOff.append(float(match.group('rms')))

    print rmsValues
    print rmsValuesNear
    print rmsValuesOff

    noise = 50.0
    rotate = 45.0
    eccentricity = 5.0
    phase = 2.0

    pb_diam = 6.8
    phs_off = phase / pb_diam

    x_axis = linspace(0,eccentricity,len(filenames))
    pl.clf()
    pl.title("RMS Levels in Eccentricity Simulation - %s" % level)
    p1 = pl.plot(x_axis,rmsValues,'b')
    p2 = pl.plot(x_axis,rmsValuesNear,'g')
    p3 = pl.plot(x_axis,rmsValuesOff,'r')
    pl.xlabel("Percent Increase in Semimajor Axis")
    pl.legend([p1,p2,p3],["Whole Sky", "Near Source", "Off Source"], loc = 2)
    pl.savefig("%s.png" % level)
