


def subsys_EDA(inputstyle,
               submf,
               dm,
               subsys,
               atomlist=[],
               spinlist=[],
               chglist=[],
               backlabel=[],
               backlist=[],
               method=['m062x', 'cc-pvtz', 'cart', 'charge'],
               coords=[],
               charges=[],
               edatype=[],
               verbose=4):
    if verbose > 4:
        logger.slog(f,"verbose = %d", verbose)
    return EDA(inputstyle, submf, dm, subsys, atomlist, spinlist, chglist,
               backlabel, backlist, method, coords, charges, edatype, verbose)


def EDA(inputstyle,
        submf,
        dm,
        subsys,
        atomlist=[],
        spinlist=[],
        chglist=[],
        backlabel=[],
        backlist=[],
        method=['m062x', 'cc-pvtz', 'cart', 'charge'],
        coords=[],
        charges=[],
        edatype='debug',
        #rctype='gebf',
        verbose=4):
    r"""
    Kwargs:
        edatype: 'full'  -> include all terms
                 'debug' -> include 'c'-related terms,
                            and 'e/e/e', 'e/e/e/e' if available
                 'smart' -> include 'debug' terms,
                            drop zero and less accurate terms
        verbose: 4 -> normal
                 5 -> print complex-frag energy
    """
    if inputstyle == 'frg':
        _cen = subsys[0]
        _env = subsys[1]
        _subsys = subsys[0] + subsys[1]
        _subsys.sort()

        s = ""
        for i in subsys[0]:
            s += (str(i) + ' ')
        cen_str = "/".join(s.split())
        #_subsys = count_from_zero(np.array(_subsys))

        ss = ""
        for i in _subsys:
            ss += (str(i) + ' ')
        sub_str = ",".join(ss.split())
    elif inputstyle == 'atomlist':
        _subsys = subsys
        cen_str = ''
        ss = ""
        for i in subsys:
            ss += (str(i) + ' ')
        sub_str = ",".join(ss.split())

    logger.slog(f,"## EDA in subsystem constructed by frag " + sub_str + " ##\n")

    EDA_data = []
    E = 0.0
    if verbose > 4:
        logger.mlog(f,"atomlist", atomlist)
    #logger.log(f,verbose)
    if verbose > 4:
        logger.slog(f,"## frag energy in real atom part ##\n")
    method_real = [item for item in method if item is not 'charge']
    # method for real atom part, removing 'charge'
    for i in _subsys:
        E = eda2_2.get_energy(inputstyle, submf, dm, [_subsys.index(i)],
                              atomlist, spinlist, chglist, method_real)
        if i in _cen:
            EDA_data.append([E[0], [[i], [], []], _cen])
        elif i in _env:
            EDA_data.append([E[0], [[], [i], []], _cen])
        else:
            logger.slog(f,"Error: label of frag not found in _cen or _env")
        if verbose > 4:
            logger.mlog(f, "  ", E[0], i, cen_stri, sep='\t\t')

    for i in _subsys:
        for j in _subsys:
            if j > i:
                E = eda2_2.get_energy(
                    inputstyle, submf, dm,
                    [_subsys.index(i), _subsys.index(j)], atomlist, spinlist,
                    chglist, method_real)
                if i in _cen:
                    if j in _cen:
                        EDA_data.append([E[0], [[i, j], [], []], _cen])
                    else:
                        EDA_data.append([E[0], [[i], [j], []], _cen])
                else:
                    if j in _cen:
                        EDA_data.append([E[0], [[j], [i], []], _cen])
                    else:
                        EDA_data.append([E[0], [[], [i, j], []], _cen])
                if verbose > 4:
                    logger.mlog(f,
                        "  ", E[0],
                        (i, j),
                        cen_str,
                        sep='\t\t')

    def ijkl_label(ijkl, _cen, _env):
        c = []  # center
        e = []  # environment
        b = []  # background
        for item in ijkl:
            if item in _cen:
                c.append(item)
            elif item in _env:
                e.append(item)
            else:
                b.append(item)
        return [c, e, b]

    for i in _subsys:
        for j in _subsys:
            for k in _subsys:
                if j > i and k > j:
                    E = eda2_2.get_energy(
                        inputstyle, submf, dm,
                        [_subsys.index(i),
                         _subsys.index(j),
                         _subsys.index(k)], atomlist, spinlist, chglist,
                        method_real)
                    EDA_data.append(
                        [E[0], ijkl_label([i, j, k], _cen, _env), _cen])
                    if verbose > 4:
                        logger.mlog(f,
                            "  ", E[0],
                            (i, j, k),
                            cen_str,
                            sep='\t\t')

    for i in _subsys:
        for j in _subsys:
            for k in _subsys:
                for l in _subsys:
                    if j > i and k > j and l > k:
                        E = eda2_2.get_energy(inputstyle, submf, dm, [
                            _subsys.index(i),
                            _subsys.index(j),
                            _subsys.index(k),
                            _subsys.index(l)
                        ], atomlist, spinlist, chglist, method_real)
                        EDA_data.append(
                            [E[0],
                             ijkl_label([i, j, k, l], _cen, _env), _cen])
                        if verbose > 4:
                            logger.log(f,
                                "  ", E[0],
                                (i, j, k, l),
                                cen_str,
                                sep='\t\t')

    logger.slog(f,"## inter energy in real atom part ##")
    EDA_RR_data = eda2_2.data2inter(EDA_data)  # real-real interaction
    logger.slog(f,"---------------------------------------------------------")
    logger.slog(f,'        Energy       interaction     center')
    for item in EDA_RR_data:
        logger.mlog(f,"  ", item[0], item[1], item[2], sep='\t\t')
    logger.slog(f,"---------------------------------------------------------\n\n")

    if 'charge' in method:
        logger.log(f,"coords", coords)
        logger.slog(f,"## iter energy between real atom frags and charges ##")
        logger.slog(f,"---------------------------------------------------------")
        logger.slog(f,'        Energy       interaction     center')
        EDA_RC_data = []  # real atom-charge interaction
        logger.mlog(f, "backlabel", backlabel)
        for i in _subsys:
            if edatype == 'subsys':
                if i not in _cen:
                    continue
            for j in backlabel:
                backrange = backlist[backlabel.index(j)]
                logger.mlog(f, "atomlist/backrange of subsys", atomlist[_subsys.index(i)], backrange)
                E = eda2_2.RC_inter(inputstyle, submf, dm, [_subsys.index(i)],
                                    atomlist, spinlist, chglist, 'qmmm',
                                    coords[np.ix_(backrange)],
                                    charges[np.ix_(backrange)])
                EDA_RC_data.append([E, ijkl_label([i, j], _cen, _env), _cen])
                logger.mlog(f,
                    "  ", E, (i, j), "*", cen_str, sep='\t\t')

        for i in _subsys:
            for j in _subsys:
                if j > i:
                    if edatype == 'subsys':
                        if i not in _cen and j not in _cen:
                            continue
                    for k in backlabel:
                        backrange = backlist[backlabel.index(k)]
                        E = eda2_2.RC_inter(
                            inputstyle, submf, dm,
                            [_subsys.index(i),
                             _subsys.index(j)], atomlist, spinlist, chglist,
                            'qmmm', coords[np.ix_(backrange)],
                            charges[np.ix_(backrange)])
                        EDA_RC_data.append(
                            [E, ijkl_label([i, j, k], _cen, _env), _cen])
                        logger.mlog(f,
                            "  ", E,
                            (i, j, k), "*",
                            cen_str,
                            sep='\t\t')
        for i in _subsys:
            if edatype == 'subsys':
                if i not in _cen:
                    continue
            for j in backlabel:
                for k in backlabel:
                    if k > j:
                        #backrange = backlist[backlabel.index(j)] \
                        #            + backlist[backlabel.index(k)]
                        #E = eda.RC_inter(inputstyle, submf, dm,
                        #                 [_subsys.index(i)], atomlist,
                        #                 spinlist, chglist, 'qmmm',
                        #                 coords[np.ix_(backrange)],
                        #                 charges[np.ix_(backrange)])
                        E = 0.0
                        EDA_RC_data.append(
                            [E, ijkl_label([i, j, k], _cen, _env), _cen])
                        logger.mlog(f,
                            "  ", E,
                            (i, j, k), "*",
                            cen_str,
                            sep='\t\t')

        for i in _subsys:
            for j in _subsys:
                for k in _subsys:
                    if edatype == 'subsys':
                        if i not in _cen and j not in _cen and k not in _cen:
                            continue
                    if j > i and k > j:
                        for l in backlabel:
                            #backrange = backlist[backlabel.index(l)]
                            #E = eda.RC_inter(inputstyle, submf, dm, [
                            #    _subsys.index(i),
                            #    _subsys.index(j),
                            #    _subsys.index(k)], atomlist, spinlist, chglist,
                            #                 'qmmm', coords[np.ix_(backrange)],
                            #                 charges[np.ix_(backrange)])
                            E = 0.0
                            EDA_RC_data.append([
                                E,
                                ijkl_label([i, j, k, l], _cen, _env), _cen
                            ])
                            logger.mlog(f,
                                "  ", E,
                                (i, j, k, l), "*",
                                cen_str,
                                sep='\t\t')
        for i in _subsys:
            for j in _subsys:
                if j > i:
                    if edatype == 'subsys':
                        if i not in _cen and j not in _cen:
                            continue
                    for k in backlabel:
                        for l in backlabel:
                            if l > k:
                                #backrange = backlist[backlabel.index(k)] \
                                #            + backlist[backlabel.index(l)]
                                #E = eda.RC_inter(
                                #    inputstyle, submf, dm,
                                #    [_subsys.index(i),
                                #     _subsys.index(j)], atomlist, spinlist,
                                #    chglist, 'qmmm', coords[np.ix_(backrange)],
                                #    charges[np.ix_(backrange)])
                                E = 0.0
                                EDA_RC_data.append([
                                    E,
                                    ijkl_label([i, j, k, l], _cen, _env), _cen
                                ])
                                logger.mlog(f,
                                    "  ", E,
                                    (i, j, k, l), "*",
                                    cen_str,
                                    sep='\t\t')
        for i in _subsys:
            if edatype == 'subsys':
                if i not in _cen:
                    continue
            for j in backlabel:
                for k in backlabel:
                    for l in backlabel:
                        if k > j and l > k:
                            #backrange = backlist[backlabel.index(j)] \
                            #            + backlist[backlabel.index(k)] \
                            #            + backlist[backlabel.index(l)]
                            #E = eda.RC_inter(inputstyle, submf, dm,
                            #                 [_subsys.index(i)], atomlist,
                            #                 spinlist, chglist, 'qmmm',
                            #                 coords[np.ix_(backrange)],
                            #                 charges[np.ix_(backrange)])
                            E = 0.0
                            EDA_RC_data.append([
                                E,
                                ijkl_label([i, j, k, l], _cen, _env), _cen
                            ])
                            logger.mlog(f,
                                "  ", E,
                                (i, j, k, l), "*",
                                cen_str,
                                sep='\t\t')
        logger.slog(f,"---------------------------------------------------------\n\n")
    if 'charge' in method:
        return EDA_RR_data + EDA_RC_data
    else:
        return EDA_RR_data


def NN_input(inputstyle,
             submf,
             dm,
             subsys,
             atomlist=[],
             spinlist=[],
             chglist=[],
             method=[
                 'm062x',
                 'cc-pvtz',
                 'cart',
             ]):
    EDA_RR_data = EDA(
        inputstyle,
        submf,
        dm,
        subsys,
        atomlist,
        spinlist,
        chglist,
        method=method)
    sumE = 0.0
    for item in EDA_RR_data:
        sumE += item[0]
    logger.mlog(f,"## tot_EDA_sum for subsys ", subsys, " = %18.10f ##" % sumE)
    EDA_NN = []
    for item in subsys:
        E = 0.0
        for jtem in EDA_RR_data:
                if len(jlabel) == 2:
                    E += jtem[0] * 0.5
                if len(jlabel) == 3:
                    E += jtem[0] * (1 / 3)
                if len(jlabel) == 4:
                    E += jtem[0] * 0.25
        EDA_NN.append([E, item])
    return EDA_NN