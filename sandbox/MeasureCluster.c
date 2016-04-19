01484 static int vmd_measure_cluster(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
01485   int i,j;
01486   // initialize optional arguments to default values
01487   int algorithm=0; // will be MEASURE_CLUSTER_QT when finished
01488   int likeness=MEASURE_DIST_FITRMSD;
01489   int numcluster=5;
01490   double cutoff=1.0;
01491   float *weights=NULL;
01492   int selupdate=0;
01493   int first=0, last=-1, step=1;
01494   int rc;
01495 
01496   // argument error message
01497   const char *argerrmsg = "<sel> [num <#clusters>] [distfunc <flag>] "
01498     "[cutoff <cutoff>] [first <first>] [last <last>] [step <step>] "
01499     "[selupdate <bool>] [weight <weights>]";
01500 
01501   // Two atom selections and optional keyword/value pairs.
01502   if ((argc < 2) || (argc > 19) || ((argc-1) % 2 == 0) )  {
01503     Tcl_WrongNumArgs(interp, 2, objv-1, (char *)argerrmsg);
01504     return TCL_ERROR;
01505   }
01506 
01507   // check atom selection
01508   AtomSel *sel = tcl_commands_get_sel(interp, Tcl_GetStringFromObj(objv[1], NULL));
01509   if (!sel) {
01510     Tcl_AppendResult(interp, "measure cluster: invalid atom selection", NULL);
01511     return TCL_ERROR;
01512   }
01513   if (!app->molecule_valid_id(sel->molid())) return MEASURE_ERR_NOMOLECULE;
01514 
01515   // parse optional arguments
01516   for (i=2; i<argc; i+=2) {
01517     const char *opt = Tcl_GetStringFromObj(objv[i], NULL);
01518     if (i==(argc-1)) {
01519       Tcl_WrongNumArgs(interp, 2, objv-1, (char *)argerrmsg);
01520       return TCL_ERROR;
01521     }
01522     if (!strcmp(opt, "num")) {
01523       if (Tcl_GetIntFromObj(interp, objv[i+1], &numcluster) != TCL_OK)
01524         return TCL_ERROR;
01525       if (numcluster < 1) {
01526         Tcl_AppendResult(interp, "measure cluster: invalid 'num' value (cannot be smaller than 1)", NULL);
01527         return TCL_ERROR;
01528       }
01529     } else if (!strcmp(opt, "cutoff")) {
01530       if (Tcl_GetDoubleFromObj(interp, objv[i+1], &cutoff) != TCL_OK)
01531         return TCL_ERROR;
01532       if (cutoff <= 0.0) {
01533         Tcl_AppendResult(interp, "measure cluster: invalid 'cutoff' value (should be larger than 0.0)", NULL);
01534         return TCL_ERROR;
01535       }
01536     } else if (!strcmp(opt, "distfunc")) {
01537       char *argstr = Tcl_GetStringFromObj(objv[i+1], NULL);
01538       if (!strcmp(argstr,"rmsd")) {
01539         likeness = MEASURE_DIST_RMSD;
01540       } else if (!strcmp(argstr,"fitrmsd")) {
01541         likeness = MEASURE_DIST_FITRMSD;
01542       } else if (!strcmp(argstr,"rgyrd")) {
01543         likeness = MEASURE_DIST_RGYRD;
01544       } else {
01545         Tcl_AppendResult(interp, "measure cluster: unknown distance function (supported are 'rmsd', 'rgyrd' and 'fitrmsd')", NULL);
01546         return TCL_ERROR;
01547       }
01548     } else if (!strcmp(opt, "selupdate")) {
01549       if (Tcl_GetBooleanFromObj(interp, objv[i+1], &selupdate) != TCL_OK)
01550         return TCL_ERROR;
01551     } else if (!strcmp(opt, "weight")) {
01552       // NOTE: we cannot use tcl_get_weights here, since we have to 
01553       // get a full (all atoms) list of weights as we may be updating
01554       // the selection in the process. Also we don't support explicit
01555       // lists of weights for now.
01556       const char *weight_string = Tcl_GetStringFromObj(objv[i+1], NULL);
01557       weights = new float[sel->num_atoms];
01558 
01559       if (!weight_string || !strcmp(weight_string, "none")) {
01560         for (j=0; j<sel->num_atoms; j++) weights[j]=1.0f;
01561       } else {
01562         // if a selection string was given, check the symbol table
01563         SymbolTable *atomSelParser = app->atomSelParser; 
01564         // weights must return floating point values, so the symbol must not 
01565         // be a singleword, so macro is NULL.
01566         atomsel_ctxt context(atomSelParser,
01567                              app->moleculeList->mol_from_id(sel->molid()), 
01568                              sel->which_frame, NULL);
01569 
01570         int fctn = atomSelParser->find_attribute(weight_string);
01571         if (fctn >= 0) {
01572           // the keyword exists, so get the data
01573           // first, check to see that the function returns floats.
01574           // if it doesn't, it makes no sense to use it as a weight
01575           if (atomSelParser->fctns.data(fctn)->returns_a != SymbolTableElement::IS_FLOAT) {
01576             Tcl_AppendResult(interp, "weight attribute must have floating point values", NULL);
01577             delete [] weights;
01578             return MEASURE_ERR_BADWEIGHTPARM;  // can't understand weight parameter 
01579           }
01580 
01581           double *tmp_data = new double[sel->num_atoms];
01582           int *all_on = new int[sel->num_atoms];
01583           for (j=0; j<sel->num_atoms; j++) all_on[j]=1;
01584 
01585           atomSelParser->fctns.data(fctn)->keyword_double(
01586             &context, sel->num_atoms, tmp_data, all_on);
01587           
01588           for (j=0; j<sel->num_atoms; j++) weights[j] = (float)tmp_data[j];
01589           
01590           // clean up.
01591           delete [] tmp_data;
01592           delete [] all_on;
01593         }
01594       }
01595     } else if (!strcmp(opt, "first")) {
01596       if (Tcl_GetIntFromObj(interp, objv[i+1], &first) != TCL_OK)
01597         return TCL_ERROR;
01598     } else if (!strcmp(opt, "last")) {
01599       if (Tcl_GetIntFromObj(interp, objv[i+1], &last) != TCL_OK)
01600         return TCL_ERROR;
01601     } else if (!strcmp(opt, "step")) {
01602       if (Tcl_GetIntFromObj(interp, objv[i+1], &step) != TCL_OK)
01603         return TCL_ERROR;
01604     } else { // unknown keyword.
01605       Tcl_AppendResult(interp, "unknown keyword '", opt, "'. usage: measure cluster ", argerrmsg,NULL);
01606       return TCL_ERROR;
01607     }
01608   }
01609 
01610   // set default for weights if not already defined
01611   if (!weights) {
01612     weights = new float[sel->num_atoms];
01613     for (j=0; j<sel->num_atoms; j++) weights[j]=1.0f;
01614   }
01615   
01616   // Allocate temporary result storage. we add one more cluster 
01617   // slot for collecting unclustered frames in an additional "cluster".
01618   // NOTE: the individual cluster lists are going to
01619   //       allocated in the ancilliary code.
01620   int  *clustersize = new int  [numcluster+1];
01621   int **clusterlist = new int *[numcluster+1];
01622   
01623   // do the cluster analysis
01624   rc = measure_cluster(sel, app->moleculeList, numcluster, algorithm, likeness, cutoff, 
01625                        clustersize, clusterlist, first, last, step, selupdate, weights);
01626 
01627   if (weights) delete [] weights;
01628 
01629   // XXX: this needs a 'case' structure to provide more meaninful error messages.
01630   if (rc != MEASURE_NOERR) { 
01631     Tcl_AppendResult(interp, "measure cluster: error during cluster analysis calculation.", NULL);
01632     return TCL_ERROR;
01633   }
01634 
01635   // convert the results of the lowlevel call to tcl lists
01636   // and build a list from them as return value.
01637   Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
01638 
01639   for (i=0; i <= numcluster; ++i) { 
01640     int j;
01641     Tcl_Obj *tcl_clist  = Tcl_NewListObj(0, NULL);
01642     for (j=0; j < clustersize[i]; ++j) {
01643       Tcl_ListObjAppendElement(interp, tcl_clist, Tcl_NewIntObj(clusterlist[i][j]));
01644     }
01645     Tcl_ListObjAppendElement(interp, tcl_result, tcl_clist);
01646   }
01647   Tcl_SetObjResult(interp, tcl_result);
01648 
01649   // free temporary result storage
01650   for (i=0; i <= numcluster; ++i)
01651     delete[] clusterlist[i];
01652 
01653   delete[] clusterlist;
01654   delete[] clustersize;
01655 
01656   return TCL_OK;
01657 }
01658 
01660 static int vmd_measure_clustsize(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
01661   int i;
01662   // initialize optional arguments to default values
01663   double cutoff=3.0; // arbitrary. took hbond cutoff.
01664   char *storesize=NULL;
01665   char *storenum=NULL;
01666   int usepbc=0;
01667   int minsize=2;
01668   int numshared=1;
01669   int rc;
01670 
01671   // argument error message
01672   const char *argerrmsg = "<sel> [cutoff <float>] [minsize <num>] [numshared <num>] "
01673     "[usepbc <bool>] [storesize <fieldname>] [storenum <fieldname>]";
01674 
01675   // Two atom selections and optional keyword/value pairs.
01676   if ((argc < 2) || (argc > 13) || ((argc-1) % 2 == 0) )  {
01677     Tcl_WrongNumArgs(interp, 2, objv-1, (char *)argerrmsg);
01678     return TCL_ERROR;
01679   }
01680 
01681   // check atom selection
01682   AtomSel *sel = tcl_commands_get_sel(interp, Tcl_GetStringFromObj(objv[1], NULL));
01683   if (!sel) {
01684     Tcl_AppendResult(interp, "measure clustsize: invalid atom selection", NULL);
01685     return TCL_ERROR;
01686   }
01687   if (!app->molecule_valid_id(sel->molid())) return MEASURE_ERR_NOMOLECULE;
01688 
01689   // parse optional arguments
01690   for (i=2; i<argc; i+=2) {
01691     const char *opt = Tcl_GetStringFromObj(objv[i], NULL);
01692     if (i==(argc-1)) {
01693       Tcl_WrongNumArgs(interp, 2, objv-1, (char *)argerrmsg);
01694       return TCL_ERROR;
01695     } else if (!strcmp(opt, "cutoff")) {
01696       if (Tcl_GetDoubleFromObj(interp, objv[i+1], &cutoff) != TCL_OK)
01697         return TCL_ERROR;
01698       if (cutoff <= 0.0) {
01699         Tcl_AppendResult(interp, "measure clustsize: invalid 'cutoff' value", NULL);
01700         return TCL_ERROR;
01701       }
01702     } else if (!strcmp(opt, "minsize")) {
01703       if (Tcl_GetIntFromObj(interp, objv[i+1], &minsize) != TCL_OK)
01704         return TCL_ERROR;
01705       if (minsize < 2) {
01706         Tcl_AppendResult(interp, "measure clustsize: invalid 'minsize' value", NULL);
01707         return TCL_ERROR;
01708       }
01709     } else if (!strcmp(opt, "numshared")) {
01710       if (Tcl_GetIntFromObj(interp, objv[i+1], &numshared) != TCL_OK)
01711         return TCL_ERROR;
01712       if (numshared < 0) {
01713         Tcl_AppendResult(interp, "measure clustsize: invalid 'numshared' value", NULL);
01714         return TCL_ERROR;
01715       }
01716     } else if (!strcmp(opt, "usepbc")) {
01717       if (Tcl_GetBooleanFromObj(interp, objv[i+1], &usepbc) != TCL_OK)
01718         return TCL_ERROR;
01719     } else if (!strcmp(opt, "storenum")) {
01720       storenum = Tcl_GetStringFromObj(objv[i+1], NULL);
01721     } else if (!strcmp(opt, "storesize")) {
01722       storesize = Tcl_GetStringFromObj(objv[i+1], NULL);
01723     } else { // unknown keyword.
01724       Tcl_AppendResult(interp, "unknown keyword '", opt, "'. usage: measure clustsize ", argerrmsg,NULL);
01725       return TCL_ERROR;
01726     }
01727   }
01728 
01729   if (usepbc) { 
01730     Tcl_AppendResult(interp, "measure clustsize: does not support periodic boundaries yet.", NULL);
01731     return TCL_ERROR;
01732   }
01733 
01734   // allocate temporary result storage
01735   // NOTE: the individual cluster lists are going to
01736   //       allocated in the ancilliary code.
01737   int num_selected=sel->selected;
01738   int *clustersize = new int[num_selected];
01739   int *clusternum= new int [num_selected];
01740   int *clusteridx= new int [num_selected];
01741   for (i=0; i < num_selected; i++) {
01742     clustersize[i] = 0;
01743     clusternum[i]  = -1;
01744     clusteridx[i]  = -1;
01745   }
01746   
01747   // do the cluster analysis
01748   rc = measure_clustsize(sel, app->moleculeList, cutoff, 
01749                          clustersize, clusternum, clusteridx,
01750                          minsize, numshared, usepbc);
01751 
01752   // XXX: this needs a 'case' structure to provide more meaninful error messages.
01753   if (rc != MEASURE_NOERR) { 
01754     Tcl_AppendResult(interp, "measure clustsize: error during cluster size analysis calculation.", NULL);
01755     return TCL_ERROR;
01756   }
01757 
01758 
01759   if (storenum || storesize) {
01760     // field names were given to store the results. check the keywords and so on.
01761     SymbolTable *atomSelParser = app->atomSelParser; 
01762     atomsel_ctxt context(atomSelParser, 
01763                          app->moleculeList->mol_from_id(sel->molid()), 
01764                          sel->which_frame, NULL);
01765 
01766     // the keyword exists, set the data
01767     if (storenum) {
01768       int fctn = atomSelParser->find_attribute(storenum);
01769       if (fctn >= 0) {
01770         if (atomSelParser->fctns.data(fctn)->returns_a == SymbolTableElement::IS_FLOAT) {
01771           double *tmp_data = new double[sel->num_atoms];
01772           int j=0;
01773           for (int i=0; i<sel->num_atoms; i++) {
01774             if (sel->on[i])
01775               tmp_data[i] = (double) clusternum[j++];
01776           }
01777           atomSelParser->fctns.data(fctn)->set_keyword_double(&context, 
01778                                                               sel->num_atoms,
01779                                                               tmp_data, sel->on);
01780           delete[] tmp_data;
01781           
01782         } else if (atomSelParser->fctns.data(fctn)->returns_a == SymbolTableElement::IS_INT) {
01783           int *tmp_data = new int[sel->num_atoms];
01784           int j=0;
01785           for (int i=0; i<sel->num_atoms; i++) {
01786             if (sel->on[i])
01787               tmp_data[i] = clusternum[j++];
01788           }
01789           atomSelParser->fctns.data(fctn)->set_keyword_int(&context, 
01790                                                            sel->num_atoms,
01791                                                            tmp_data, sel->on);
01792           delete[] tmp_data;
01793         } else {
01794           Tcl_AppendResult(interp, "measure clustsize: storenum field must accept numbers", NULL);
01795           return TCL_ERROR;
01796         }
01797       } else {
01798         Tcl_AppendResult(interp, "measure clustsize: invalid field name for storenum", NULL);
01799         return TCL_ERROR;
01800       }
01801     }
01802 
01803     // the keyword exists, set the data
01804     if (storesize) {
01805       int fctn = atomSelParser->find_attribute(storesize);
01806       if (fctn >= 0) {
01807         if (atomSelParser->fctns.data(fctn)->returns_a == SymbolTableElement::IS_FLOAT) {
01808           double *tmp_data = new double[sel->num_atoms];
01809           int j=0;
01810           for (int i=0; i<sel->num_atoms; i++) {
01811             if (sel->on[i])
01812               tmp_data[i] = (double) clustersize[j++];
01813           }
01814           atomSelParser->fctns.data(fctn)->set_keyword_double(&context, 
01815                                                               sel->num_atoms,
01816                                                               tmp_data, sel->on);
01817           delete[] tmp_data;
01818           
01819         } else if (atomSelParser->fctns.data(fctn)->returns_a == SymbolTableElement::IS_INT) {
01820           int *tmp_data = new int[sel->num_atoms];
01821           int j=0;
01822           for (int i=0; i<sel->num_atoms; i++) {
01823             if (sel->on[i])
01824               tmp_data[i] = clustersize[j++];
01825           }
01826           atomSelParser->fctns.data(fctn)->set_keyword_int(&context, 
01827                                                            sel->num_atoms,
01828                                                            tmp_data, sel->on);
01829           delete[] tmp_data;
01830         } else {
01831           Tcl_AppendResult(interp, "measure clustsize: storenum field must accept numbers", NULL);
01832           return TCL_ERROR;
01833         }
01834       } else {
01835         Tcl_AppendResult(interp, "measure clustsize: invalid field name for storesize", NULL);
01836         return TCL_ERROR;
01837       }
01838     }
01839   } else {
01840     // convert the results of the lowlevel call to tcl lists
01841     // and build a list from them as return value.
01842     Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
01843 
01844     Tcl_Obj *tcl_ilist  = Tcl_NewListObj(0, NULL);
01845     Tcl_Obj *tcl_clist  = Tcl_NewListObj(0, NULL);
01846     Tcl_Obj *tcl_nlist  = Tcl_NewListObj(0, NULL);
01847     for (i=0; i<num_selected; i++) { 
01848       Tcl_ListObjAppendElement(interp, tcl_ilist, Tcl_NewIntObj(clusteridx[i]));
01849       Tcl_ListObjAppendElement(interp, tcl_clist, Tcl_NewIntObj(clustersize[i]));
01850       Tcl_ListObjAppendElement(interp, tcl_nlist, Tcl_NewIntObj(clusternum[i]));
01851     }
01852     Tcl_ListObjAppendElement(interp, tcl_result, tcl_ilist);
01853     Tcl_ListObjAppendElement(interp, tcl_result, tcl_clist);
01854     Tcl_ListObjAppendElement(interp, tcl_result, tcl_nlist);
01855     Tcl_SetObjResult(interp, tcl_result);
01856   }
01857   
01858   delete[] clustersize;
01859   delete[] clusternum;
01860   delete[] clusteridx;
01861 
01862   return TCL_OK;
01863 }
01864 
01865 static int vmd_measure_hbonds(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
01866   
01867   // Cutoff, angle, and either one or two atom selections
01868   if (argc != 4 && argc != 5) {
01869     Tcl_WrongNumArgs(interp, 2, objv-1, (char *)"<cutoff> <angle> <selection1> [<selection2>]");
01870     return TCL_ERROR;
01871   }
01872   AtomSel *sel1 = tcl_commands_get_sel(interp, Tcl_GetStringFromObj(objv[3],NULL));
01873   if (!sel1) {
01874     Tcl_AppendResult(interp, "measure hbonds: invalid first atom selection", NULL);
01875     return TCL_ERROR;
01876   }
01877 
01878   AtomSel *sel2 = NULL;
01879   if (argc == 5) {
01880     sel2 = tcl_commands_get_sel(interp, Tcl_GetStringFromObj(objv[4],NULL));
01881     if (!sel2) {
01882       Tcl_AppendResult(interp, "measure hbonds: invalid second atom selection", NULL);
01883       return TCL_ERROR;
01884     }
01885   }
01886   if (sel2 && sel2->molid() != sel1->molid()) {
01887     Tcl_AppendResult(interp, "measure hbonds: error, atom selections must come from same molecule.", NULL);
01888     return TCL_ERROR;
01889   }
01890   double cutoff;
01891   if (Tcl_GetDoubleFromObj(interp, objv[1], &cutoff) != TCL_OK) 
01892     return TCL_ERROR;
01893 
01894   double maxangle;
01895   if (Tcl_GetDoubleFromObj(interp, objv[2], &maxangle) != TCL_OK) 
01896     return TCL_ERROR;
01897   
01898   const float *pos = sel1->coordinates(app->moleculeList);
01899   if (!pos) {
01900     Tcl_AppendResult(interp, "measure bondsearch: error, molecule contains no coordinates", NULL);
01901     return TCL_ERROR;
01902   }
01903 
01904   // XXX the actual code for measuring hbonds doesn't belong here, it should
01905   //     be moved into Measure.[Ch] where it really belongs.  This file
01906   //     only implements the Tcl interface, and should not be doing the
01907   //     hard core math, particularly if we want to expose the same
01908   //     feature via other scripting interfaces.  Also, having a single
01909   //     implementation avoids having different Tcl/Python bugs in the 
01910   //     long-term.  Too late to do anything about this now, but should be
01911   //     addressed for the next major version when time allows.
01912 
01913   // XXX This code is close, but not identical to the HBonds code in 
01914   //     DrawMolItem.  Is there any good reason they aren't identical?
01915   //     This version does a few extra tests that the other does not.
01916 
01917   Molecule *mol = app->moleculeList->mol_from_id(sel1->molid());
01918 
01919   const int *A = sel1->on;
01920   const int *B = sel2 ? sel2->on : sel1->on;
01921  
01922   GridSearchPair *pairlist = vmd_gridsearch2(pos, sel1->num_atoms, A, B, (float) cutoff, sel1->num_atoms * 27);
01923   GridSearchPair *p, *tmp;
01924   float donortoH[3], Htoacceptor[3];
01925   Tcl_Obj *donlist = Tcl_NewListObj(0, NULL);
01926   Tcl_Obj *hydlist = Tcl_NewListObj(0, NULL);
01927   Tcl_Obj *acclist = Tcl_NewListObj(0, NULL);
01928   for (p=pairlist; p != NULL; p=tmp) {
01929     MolAtom *a1 = mol->atom(p->ind1); 
01930     MolAtom *a2 = mol->atom(p->ind2); 
01931     
01932     // neither the donor nor acceptor may be hydrogens
01933     if (mol->atom(p->ind1)->atomType == ATOMHYDROGEN ||
01934         mol->atom(p->ind2)->atomType == ATOMHYDROGEN) {
01935       tmp = p->next;
01936       free(p);
01937       continue;
01938     } 
01939     if (!a1->bonded(p->ind2)) {
01940       int b1 = a1->bonds;
01941       int b2 = a2->bonds;
01942       const float *coor1 = pos + 3*p->ind1;
01943       const float *coor2 = pos + 3*p->ind2;
01944       int k;
01945       // first treat sel1 as donor
01946       for (k=0; k<b1; k++) {
01947         const int hindex = a1->bondTo[k];
01948         if (mol->atom(hindex)->atomType == ATOMHYDROGEN) {         
01949           const float *hydrogen = pos + 3*hindex;
01950           vec_sub(donortoH,hydrogen,coor1);
01951           vec_sub(Htoacceptor,coor2,hydrogen);
01952           if (angle(donortoH, Htoacceptor)  < maxangle ) {
01953             Tcl_ListObjAppendElement(interp, donlist, Tcl_NewIntObj(p->ind1));
01954             Tcl_ListObjAppendElement(interp, acclist, Tcl_NewIntObj(p->ind2));
01955             Tcl_ListObjAppendElement(interp, hydlist, Tcl_NewIntObj(hindex));
01956           }
01957         }
01958       }
01959       // if only one atom selection was given, treat sel2 as a donor as well
01960       if (!sel2) {
01961         for (k=0; k<b2; k++) {
01962           const int hindex = a2->bondTo[k];
01963           if (mol->atom(hindex)->atomType == ATOMHYDROGEN) {
01964             const float *hydrogen = pos + 3*hindex;
01965             vec_sub(donortoH,hydrogen,coor2);
01966             vec_sub(Htoacceptor,coor1,hydrogen);
01967             if (angle(donortoH, Htoacceptor)  < maxangle ) {
01968               Tcl_ListObjAppendElement(interp, donlist, Tcl_NewIntObj(p->ind2));
01969               Tcl_ListObjAppendElement(interp, acclist, Tcl_NewIntObj(p->ind1));
01970               Tcl_ListObjAppendElement(interp, hydlist, Tcl_NewIntObj(hindex));
01971             }
01972           }
01973         }
01974       } 
01975     }
01976     tmp = p->next;
01977     free(p);
01978   }
01979   Tcl_Obj *result = Tcl_NewListObj(0, NULL);
01980   Tcl_ListObjAppendElement(interp, result, donlist);
01981   Tcl_ListObjAppendElement(interp, result, acclist);
01982   Tcl_ListObjAppendElement(interp, result, hydlist);
01983   Tcl_SetObjResult(interp, result);
01984   return TCL_OK;
01985 }