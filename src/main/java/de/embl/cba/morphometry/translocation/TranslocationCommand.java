package de.embl.cba.morphometry.translocation;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Rois;
import de.embl.cba.morphometry.Utils;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.display.DisplayService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.io.File;
import java.util.ArrayList;

@Plugin(type = Command.class, menuPath = "Plugins>Measure Plasma Membrane Translocation" )
public class TranslocationCommand<T extends RealType<T> & NativeType< T > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter
	public DisplayService displayService;

	@Parameter
	public ImagePlus imagePlus;

	@Parameter( style = "directory" )
	public File outputDirectory;

	@Parameter
	public boolean showIntermediateResults = false;

	public void run()
	{

		final RoiManager rm = RoiManager.getInstance();

		if ( rm == null || rm.getRoisAsArray().length == 0 )
		{
			Logger.error( "This plugin requires rectangular ROIs added to the ROIManager." );
			return;
		}

		final Roi[] rois = rm.getRoisAsArray();
		final ArrayList< FinalInterval > intervals = Rois.asIntervals( rois );

		new TranslocationComputer(
				Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 ),
				intervals,
				opService );

	}

}
