/*
@@@BUILDINFO@@@ blur_script_main.jsx 1.0
*/

//
// This script loads a pair of images (RGB image and depth map) into photoshop
// and then blurs them base on the given parameters
//

/*

// BEGIN__HARVEST_EXCEPTION_ZSTRING

<javascriptresource>
<name>$$$/JavaScripts/Blur Script/Menu=Blur Images...</name>
<eventid>6F17BFA7-EFC8-40EA-B850-7B95ED8EA713</eventid>
</javascriptresource>

// END__HARVEST_EXCEPTION_ZSTRING

*/

// debug level: 0-2 (0:disable, 1:break on error, 2:break at beginning)
//$.level = (Window.version.search("d") != -1) ? 1 : 0;	// This chokes bridge
$.level = 2;

// debugger; // launch debugger on next line

// on localized builds we pull the $$$/Strings from a .dat file
$.localize = true;

// Put header files in a "Stack Scripts Only" folder.  The "...Only" tells
// PS not to place it in the menu.  For that reason, we do -not- localize that
// portion of the folder name.
var g_StackScriptFolderPath = app.path + "/"+ localize("$$$/ScriptingSupport/InstalledScripts=Presets/Scripts") + "/"
										+ localize("$$$/private/LoadStack/StackScriptOnly=Stack Scripts Only/");

// this is the root location where the support scripts are
var support_scripts_root = "D:/Projects/dfd/";

// load this file in to get the images
$.evalFile(support_scripts_root + "common/scripts/data_listing.jsx");

// this file contains all of the code to blur the image
$.evalFile(support_scripts_root + "common/scripts/blur_support_functions.jsx");

// load the support files 
$.evalFile(g_StackScriptFolderPath + "LatteUI.jsx");
$.evalFile(g_StackScriptFolderPath + "StackSupport.jsx");
$.evalFile(g_StackScriptFolderPath + "CreateImageStack.jsx");

/************************************************************/
// loadLayers routines

loadLayers = new ImageStackCreator( localize("$$$/AdobePlugin/Shared/LoadStack/Process/Name=Load Layers"),
										  localize('$$$/AdobePlugin/Shared/LoadStack/Auto/untitled=Untitled' ) );

// LoadLayers is less restrictive than MergeToHDR
loadLayers.mustBeSameSize			= false;	// Images' height & width don't need to match
loadLayers.mustBeUnmodifiedRaw		= false;	// Exposure adjustements in Camera raw are allowed
loadLayers.mustNotBe32Bit			= false;	// 32 bit images
loadLayers.createSmartObject		= false;	// If true, option to create smart object is checked.

// Override the default to use "Auto" alignment.
loadLayers.alignStack = function( stackDoc )
{
	selectAllLayers(stackDoc, 2);
	alignLayersByContent( "Auto" );
}

loadLayers.stackLayers = function()
{
	var result, i, stackDoc = null;
	
	stackDoc = this.loadStackLayers();
	if (! stackDoc)
		return;
	
	// Nuke the "destination" layer that got created (M2HDR holdover)
	stackDoc.layers[this.pluginName].remove();
	
	// Stack 'em up.
	if (this.createSmartObject)
	{
		selectAllLayers( stackDoc );
		executeAction( knewPlacedLayerStr, new ActionDescriptor(), DialogModes.NO );
	}
}

loadLayers.intoStack = function(filelist, alignFlag)
{
	if (typeof(alignFlag) == 'boolean')
		loadLayers.useAlignment = alignFlag;
		
	if (filelist.length < 2)
	{
		alert(localize("$$$/AdobeScripts/Shared/LoadLayers/AtLeast2=At least two files must be selected to create a stack."), this.pluginName, true );
		return;
	}
	var j;
	this.stackElements = new Array();
	for (j in filelist)
	{
		var f = filelist[j];
		this.stackElements.push( new StackElement( (typeof(f) == 'string') ? File(f) : f ) );
	}
		
	if (this.stackElements.length > 1)
		this.stackLayers();
}

loadLayers.getImages = function()
{
    // cycle through the image pairs listed in the image_file_pairs coming from the data_listing.jsx file
     for (idx in image_file_pairs)
	{
		var image_pairs = [File(data_directory + image_file_pairs[idx][0]), File(data_directory + image_file_pairs[idx][1])] ;
        
        // load the image pairs into the stack
        loadLayers.intoStack(image_pairs, false);
        
        var save_file_name = image_file_pairs[idx][0].substr(0,image_file_pairs[idx][0].lastIndexOf('.')) + "_ps_blur.png";
        
        // blur the image
        blur_image(255, false, 10.0);
        
        // save the image
        save_image(data_directory + save_file_name);
        
        // close the stack
        close_stack();
    }
}

if (typeof(loadLayersFromScript) == 'undefined')
{
    loadLayers.getImages();
}