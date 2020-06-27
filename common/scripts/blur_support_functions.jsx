function blur_image(blur_focal_distance, dm_invert, iris_radius)
{

    try
    {
        // =======================================================
        // This section hides the selected layer which should be the RGB image    
        var idHd = charIDToTypeID( "Hd  " );
            var desc4 = new ActionDescriptor();
            var idnull = charIDToTypeID( "null" );
                var list1 = new ActionList();
                    var ref2 = new ActionReference();
                    var idLyr = charIDToTypeID( "Lyr " );
                    var idOrdn = charIDToTypeID( "Ordn" );
                    var idTrgt = charIDToTypeID( "Trgt" );
                    ref2.putEnumerated( idLyr, idOrdn, idTrgt );
                list1.putReference( ref2 );
            desc4.putList( idnull, list1 );
        executeAction( idHd, desc4, DialogModes.NO );

        // =======================================================
        // this section selects the depthmap image and makes it not visible
        var idslct = charIDToTypeID( "slct" );
            var desc5 = new ActionDescriptor();
            var idnull = charIDToTypeID( "null" );
                var ref3 = new ActionReference();
                var idLyr = charIDToTypeID( "Lyr " );
                //ref3.putName( idLyr, "disp1.png" );
            ref3.putIndex( idLyr, 1 );
            desc5.putReference( idnull, ref3 );
            var idMkVs = charIDToTypeID( "MkVs" );
            desc5.putBoolean( idMkVs, false );
            var idLyrI = charIDToTypeID( "LyrI" );
                var list2 = new ActionList();
                list2.putInteger( 3 );
            desc5.putList( idLyrI, list2 );
        executeAction( idslct, desc5, DialogModes.NO );

        // =======================================================
        // This selects the depth map
        var idslct = charIDToTypeID( "slct" );
            var desc6 = new ActionDescriptor();
            var idnull = charIDToTypeID( "null" );
                var ref4 = new ActionReference();
                var idChnl = charIDToTypeID( "Chnl" );
                var idChnl = charIDToTypeID( "Chnl" );
                var idBl = charIDToTypeID( "Bl  " );
                ref4.putEnumerated( idChnl, idChnl, idBl );
            desc6.putReference( idnull, ref4 );
        executeAction( idslct, desc6, DialogModes.NO );
        
        // =======================================================
        // This section duplicates the blue channel selected and renames it alpha    
        var idDplc = charIDToTypeID( "Dplc" );
            var desc10 = new ActionDescriptor();
            var idnull = charIDToTypeID( "null" );
                var ref5 = new ActionReference();
                var idChnl = charIDToTypeID( "Chnl" );
                var idOrdn = charIDToTypeID( "Ordn" );
                var idTrgt = charIDToTypeID( "Trgt" );
                ref5.putEnumerated( idChnl, idOrdn, idTrgt );
            desc10.putReference( idnull, ref5 );
            var idNm = charIDToTypeID( "Nm  " );
            desc10.putString( idNm, """alpha""" );
        executeAction( idDplc, desc10, DialogModes.NO );

        // =======================================================
        // this code selects the RGB image
        var idslct = charIDToTypeID( "slct" );
            var desc13 = new ActionDescriptor();
            var idnull = charIDToTypeID( "null" );
                var ref6 = new ActionReference();
                var idLyr = charIDToTypeID( "Lyr " );
                ref6.putIndex( idLyr, 2 );
            desc13.putReference( idnull, ref6 );
            var idMkVs = charIDToTypeID( "MkVs" );
            desc13.putBoolean( idMkVs, false );
            var idLyrI = charIDToTypeID( "LyrI" );
                var list4 = new ActionList();
                list4.putInteger( 1 );
            desc13.putList( idLyrI, list4 );
        executeAction( idslct, desc13, DialogModes.NO );
    
        // =======================================================
        // this code turns the RGB image back on
        var idShw = charIDToTypeID( "Shw " );
            var desc12 = new ActionDescriptor();
            var idnull = charIDToTypeID( "null" );
                var list4 = new ActionList();
                    var ref7 = new ActionReference();
                    var idLyr = charIDToTypeID( "Lyr " );
                    var idOrdn = charIDToTypeID( "Ordn" );
                    var idTrgt = charIDToTypeID( "Trgt" );
                    ref7.putEnumerated( idLyr, idOrdn, idTrgt );
                list4.putReference( ref7 );
            desc12.putList( idnull, list4 );
        executeAction( idShw, desc12, DialogModes.NO );

        // =======================================================
        // this does the lens blur
        var idBokh = charIDToTypeID( "Bokh" );
            var desc14 = new ActionDescriptor();
            var idBkDi = charIDToTypeID( "BkDi" );
            var idBtDi = charIDToTypeID( "BtDi" );
            var idBeIa = charIDToTypeID( "BeIa" );
            desc14.putEnumerated( idBkDi, idBtDi, idBeIa );
            var idBkDc = charIDToTypeID( "BkDc" );
            desc14.putInteger( idBkDc, 0 );
            var idBkDp = charIDToTypeID( "BkDp" );
            desc14.putInteger( idBkDp, blur_focal_distance );
            var idBkDs = charIDToTypeID( "BkDs" );
            desc14.putBoolean( idBkDs, dm_invert );
            var idBkIs = charIDToTypeID( "BkIs" );
            var idBtIs = charIDToTypeID( "BtIs" );
            var idBeSeight = charIDToTypeID( "BeS8" );
            desc14.putEnumerated( idBkIs, idBtIs, idBeSeight );
            var idBkIb = charIDToTypeID( "BkIb" );
            desc14.putDouble( idBkIb, iris_radius );
            var idBkIc = charIDToTypeID( "BkIc" );
            desc14.putInteger( idBkIc, 0 );
            var idBkIr = charIDToTypeID( "BkIr" );
            desc14.putInteger( idBkIr, 0 );
            var idBkSb = charIDToTypeID( "BkSb" );
            desc14.putDouble( idBkSb, 0.000000 );
            var idBkSt = charIDToTypeID( "BkSt" );
            desc14.putInteger( idBkSt, 255 );
            var idBkNa = charIDToTypeID( "BkNa" );
            desc14.putInteger( idBkNa, 0 );
            var idBkNt = charIDToTypeID( "BkNt" );
            var idBtNt = charIDToTypeID( "BtNt" );
            var idBeNu = charIDToTypeID( "BeNu" );
            desc14.putEnumerated( idBkNt, idBtNt, idBeNu );
            var idBkNm = charIDToTypeID( "BkNm" );
            desc14.putBoolean( idBkNm, false );
        executeAction( idBokh, desc14, DialogModes.NO );
    }
    catch(e)
    {
        alert(e + ":" + e.line);
    }
}

function save_image(filename)
{
    // =======================================================
    var idsave = charIDToTypeID( "save" );
        var desc9 = new ActionDescriptor();
        var idAs = charIDToTypeID( "As  " );
            var desc10 = new ActionDescriptor();
            var idMthd = charIDToTypeID( "Mthd" );
            var idPNGMethod = stringIDToTypeID( "PNGMethod" );
            var idquick = stringIDToTypeID( "quick" );
            desc10.putEnumerated( idMthd, idPNGMethod, idquick );
            var idPGIT = charIDToTypeID( "PGIT" );
            var idPGIT = charIDToTypeID( "PGIT" );
            var idPGIN = charIDToTypeID( "PGIN" );
            desc10.putEnumerated( idPGIT, idPGIT, idPGIN );
            var idPNGf = charIDToTypeID( "PNGf" );
            var idPNGf = charIDToTypeID( "PNGf" );
            var idPGAd = charIDToTypeID( "PGAd" );
            desc10.putEnumerated( idPNGf, idPNGf, idPGAd );
            var idCmpr = charIDToTypeID( "Cmpr" );
            desc10.putInteger( idCmpr, 6 );
        var idPNGF = charIDToTypeID( "PNGF" );
        desc9.putObject( idAs, idPNGF, desc10 );
        var idIn = charIDToTypeID( "In  " );
        desc9.putPath( idIn, new File( filename ) );
        var idDocI = charIDToTypeID( "DocI" );
        desc9.putInteger( idDocI, 199 );
        var idCpy = charIDToTypeID( "Cpy " );
        desc9.putBoolean( idCpy, true );
        var idAlpC = charIDToTypeID( "AlpC" );
        desc9.putBoolean( idAlpC, false );
        var idsaveStage = stringIDToTypeID( "saveStage" );
        var idsaveStageType = stringIDToTypeID( "saveStageType" );
        var idsaveBegin = stringIDToTypeID( "saveBegin" );
        desc9.putEnumerated( idsaveStage, idsaveStageType, idsaveBegin );
    executeAction( idsave, desc9, DialogModes.NO );

}

function close_stack()
{
    // =======================================================
    var idCls = charIDToTypeID( "Cls " );
        var desc60 = new ActionDescriptor();
        var idSvng = charIDToTypeID( "Svng" );
        var idYsN = charIDToTypeID( "YsN " );
        var idN = charIDToTypeID( "N   " );
        desc60.putEnumerated( idSvng, idYsN, idN );
        var idDocI = charIDToTypeID( "DocI" );
        desc60.putInteger( idDocI, 366 );
        var idforceNotify = stringIDToTypeID( "forceNotify" );
        desc60.putBoolean( idforceNotify, true );
    executeAction( idCls, desc60, DialogModes.NO );

}
