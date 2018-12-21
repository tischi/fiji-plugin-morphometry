package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import ij.ImagePlus;
import ij.io.FileSaver;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;
import java.util.ArrayList;


@Plugin(type = Command.class, menuPath = "Plugins>Tracking>Microglia Tracking" )
public class MicrogliaTrackingCommand < T extends RealType<T> & NativeType< T > > implements Command
{
	@Parameter
	public UIService uiService;

	@Parameter
	public DatasetService datasetService;

	@Parameter
	public LogService logService;

	@Parameter
	public OpService opService;

	@Parameter
	public StatusService statusService;

	@Parameter( required = false )
	public ImagePlus imagePlus;

	MicrogliaTrackingSettings settings = new MicrogliaTrackingSettings();

	@Parameter
	public File inputFile;

	@Parameter( style = "directory" )
	public File outputDirectory;

	@Parameter ( label = "Microglia channel index", min = "1")
	public long microgliaChannelIndexOneBased = settings.microgliaChannelIndexOneBased;

	@Parameter ( label = "Minimal time frame to be processed", min = "1" )
	public long tMin = 1;

	@Parameter ( label = "Maximal time frame to be processed", min = "1" )
	public long tMax = 1000000;

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;


	public void run()
	{
		processFile( inputFile );
	}

	private void processFile( File file )
	{
		ImagePlus imagePlus = ImageIO.openWithBioFormats( file.getAbsolutePath() );

		if ( imagePlus == null )
		{
			Utils.error( "Could not open image: " + file );
		}

		final MicrogliaTracking microgliaTracking = new MicrogliaTracking( imagePlus, showIntermediateResults, opService, microgliaChannelIndexOneBased, tMin, tMax );

		microgliaTracking.run();

		final ArrayList< RandomAccessibleInterval< T > > labelings = microgliaTracking.getLabelings();

		final ImagePlus labelsImp = Utils.labelingsAsImagePlus( labelings );

		final String outputPath = outputDirectory + File.separator + inputFile.getName().split( "\\." )[ 0 ] + "-labelMasks.tif";

		new FileSaver( labelsImp ).saveAsTiff( outputPath );

		logService.info( "Results saved: " + outputPath );

	}


}
