execfile('run_test.py')

def rnp():
    run_test(tests=['no_perturbation'], run=True, plot=False)

def rnpnm():
    run_test(tests=['no_perturbation_no_make'], run=True, plot=False)

def rn():
    run_test(tests=['noise'], run=True, plot=False)

def pn():
    run_test(tests=['noise'], run=False, plot=True)

def rp():
    run_test(tests=['pointing'], run=True, plot=False)

def rpnm():
    run_test(tests=['pointing_no_make'], run=True, plot=False)

