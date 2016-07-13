def writeNMD(filename, modes, atoms):
    """Return *filename* that contains *modes* and *atoms* data in NMD format
    described in :ref:`nmd-format`.  :file:`.nmd` extension is appended to
    filename, if it does not have an extension.

    .. note::
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`.Vector` instance is given, it will be normalized
          before it is written. It's length before normalization will be
          written as the scaling factor of the vector."""

    if not isinstance(modes, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, '
                        'not {0}'.format(type(modes)))
    if modes.numAtoms() != atoms.numAtoms():
        raise Exception('number of atoms do not match')
    out = openFile(addext(filename, '.nmd'), 'w')

    #out.write('#!{0} -e\n'.format(VMDPATH))
    out.write('nmwiz_load {0}\n'.format(abspath(filename)))
    name = modes.getTitle()
    name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = str(atoms)
        name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = splitext(split(filename)[1])[0]
    out.write('name {0}\n'.format(name))
    try:
        coords = atoms.getCoords()
    except:
        raise ValueError('coordinates could not be retrieved from atoms')
    if coords is None:
        raise ValueError('atom coordinates are not set')

    try:
        data = atoms.getNames()
        if data is not None:
            out.write('atomnames {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnames()
        if data is not None:
            out.write('resnames {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnums()
        if data is not None:
            out.write('resids ')
            data.tofile(out, ' ')
            out.write('\n')
    except:
        pass
    try:
        data = atoms.getChids()
        if data is not None:
            out.write('chainids {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getSegnames()
        if data is not None:
            out.write('segnames {0}\n'.format(' '.join(data)))
    except:
        pass

    try:
        data = atoms.getBetas()
        if data is not None:
            out.write('bfactors ')
            data.tofile(out, ' ', '%.2f')
            out.write('\n')
    except:
        pass

    format = '{0:.3f}'.format
    out.write('coordinates ')
    coords.tofile(out, ' ', '%.3f')
    out.write('\n')
    count = 0
    if isinstance(modes, Vector):
        out.write('mode 1 {0:.2f} '.format(abs(modes)))
        modes.getNormed()._getArray().tofile(out, ' ', '%.3f')
        out.write('\n')
        count += 1
    else:
        if isinstance(modes, Mode):
            modes = [modes]
        for mode in modes:
            if mode.getEigval() < ZERO:
                continue
            out.write('mode {0} {1:.2f} '.format(
                       mode.getIndex()+1, mode.getVariance()**0.5))
            arr = mode._getArray().tofile(out, ' ', '%.3f')
            out.write('\n')
            count += 1
    if count == 0:
        LOGGER.warning('No normal mode data was written. '
                       'Given modes might have 0 eigenvalues.')
    out.close()
    return filename