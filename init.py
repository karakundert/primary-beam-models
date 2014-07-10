source_dir = "/lustre/kkundert/Code/"

execfile(source_dir+"run_test.py")

def rnp():
    run_test(tests=['no_perturbation'], run=True, plot=False)

def rnpnm():
    run_test(tests=['no_perturbation_no_make'], run=True, plot=False)

def rn():
    run_test(tests=['noise'], run=True, plot=False)

def rnwp():
    run_test(tests=['noise'], run=True, plot=True)

def pn():
    run_test(tests=['noise'], run=False, plot=True)

def rp():
    run_test(tests=['pointing'], run=True, plot=False)

def rpnm():
    run_test(tests=['pointing_no_make'], run=True, plot=False)

def rr():
    run_test(tests=['rot_aper'], run=True, plot=False)
