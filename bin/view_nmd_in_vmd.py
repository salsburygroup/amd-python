def viewNMDinVMD(filename):
    """Start VMD in the current Python session and load NMD data."""

    vmd = pathVMD()
    if vmd:
        os.system('{0} -e {1}'.format(vmd, abspath(filename)))