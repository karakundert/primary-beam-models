from glob import glob
from matplotlib import axes
# choose the comparable data set
filenames = glob('noise/Data*')

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

rmsValues = []
rmsValuesNear = []
rmsValuesOff = []
xvals = xrange(len(filenames))
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

x_axis = arange(0,noise,noise/len(filenames))
pl.clf()
pl.title("RMS Levels in Noise Simulation")
p1 = pl.plot(x_axis,rmsValues,'b')
p2 = pl.plot(x_axis,rmsValuesNear,'g')
p3 = pl.plot(x_axis,rmsValuesOff,'r')
#pl.xlim(0,noise,autoscale=True)
pl.xlabel("Percent of Amplitude")
pl.legend([p1,p2,p3],["Whole Sky", "Near Source", "Off Source"], loc = 2)
