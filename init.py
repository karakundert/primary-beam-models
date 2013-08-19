execfile('makesources.py')
execfile('run_test.py')

def rn():
    run_test(tests=['noise'], run=True, plot=True)

def pn():
    run_test(tests=['noise'], run=False, plot=True)

def rp():
    run_test(tests=['pointing'], run=True, plot=True)

def pp():
    run_test(tests=['pointing'], run=False, plot=True)
