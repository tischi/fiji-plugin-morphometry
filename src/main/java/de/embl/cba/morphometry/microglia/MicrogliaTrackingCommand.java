package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import ij.ImagePlus;
import ij.io.FileSaver;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.io.File;
import java.util.ArrayList;


@Plugin(type = Command.class, menuPath = "Plugins>Tracking>Microglia Tracking" )
public class MicrogliaTrackingCommand < T extends RealType<T> & NativeType< T > > implements Command
{
	@Parameter
	public LogService logService;

	@Parameter
	public OpService opService;

	MicrogliaTrackingSettings settings = new MicrogliaTrackingSettings();

	@Parameter( label = "Input time series (must be 2D and single channel)")
	public File inputIntensitiesFile;

//	@Parameter( label = "Channel to be segmented and tracked", min = "1")
//	public long microgliaChannelIndexOneBased = settings.microgliaChannelIndexOneBased;

	@Parameter( label = "Output directory", style = "directory" )
	public File outputDirectory;

	@Parameter( label = "Minimal time frame to be processed", min = "1" )
	public long tMinOneBased = 1;

	@Parameter( label = "Maximal time frame to be processed", min = "1" )
	public long tMaxOneBased = 1000000;

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;


	public void run()
	{
		processFile( inputIntensitiesFile );
	}

	private void processFile( File file )
	{
		ImagePlus imagePlus = ImageIO.openWithBioFormats( file.getAbsolutePath() );

		if ( imagePlus == null )
		{
			Utils.error( "Could not open image: " + file );
			return;
		}

		saveLabels( computeLabels( imagePlus ) );

	}

	private ArrayList< RandomAccessibleInterval< T > > computeLabels( ImagePlus imagePlus )
	{
		final MicrogliaTracking microgliaTracking = new MicrogliaTracking( imagePlus, showIntermediateResults, opService, 1, tMinOneBased, tMaxOneBased );

		microgliaTracking.run();

		return (ArrayList< RandomAccessibleInterval< T > > ) microgliaTracking.getLabelings();
	}

	private void saveLabels( ArrayList< RandomAccessibleInterval< T > > labelings )
	{
		final ImagePlus labelsImp = Utils.labelingsAsImagePlus( labelings );

		final String outputPath = outputDirectory + File.separator + inputIntensitiesFile.getName().split( "\\." )[ 0 ] + "-labelMasks.tif";

		new FileSaver( labelsImp ).saveAsTiff( outputPath );

		Logger.log( "Segmented and tracked label images saved: " + outputPath );
	}

}
