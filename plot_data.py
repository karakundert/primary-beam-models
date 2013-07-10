from glob import glob
# choose the comparable data set
filenames = glob('rotate/Data*')

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

pl.clf()
pl.title("RMS Levels in Rotation Simulation")
p1 = pl.plot(xvals,rmsValues,'b')
p2 = pl.plot(xvals,rmsValuesNear,'g')
p3 = pl.plot(xvals,rmsValuesOff,'r')
pl.legend([p1,p2,p3],["Whole Sky", "Near Source", "Off Source"], loc = 2)
