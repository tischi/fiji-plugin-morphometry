package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaSegmentationAndTracking;
import de.embl.cba.morphometry.microglia.MicrogliaSegmentationAndTrackingSettings;
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


@Plugin(type = Command.class, menuPath = "Plugins>Tracking>Microglia Segmentation And Tracking" )
public class MicrogliaSegmentationAndTrackingCommand< T extends RealType<T> & NativeType< T > > implements Command
{
	@Parameter
	public LogService logService;

	@Parameter
	public OpService opService;

	MicrogliaSegmentationAndTrackingSettings settings = new MicrogliaSegmentationAndTrackingSettings();

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
	private ImagePlus imagePlus;
	private ArrayList< RandomAccessibleInterval< T > > intensities;

	public void run()
	{
		processFile( inputIntensitiesFile );
	}

	private void processFile( File file )
	{
		openIntensitiesAsFrameList( file );

		saveLabels( computeLabels( intensities, Utils.get2dCalibration( imagePlus ) ) );
	}

	private void openIntensitiesAsFrameList( File file )
	{
		imagePlus = ImageIO.openWithBioFormats( file.getAbsolutePath() );

		if ( imagePlus == null )
		{
			Logger.error( "Could not open image: " + file );
			return;
		}

		if ( imagePlus.getNChannels() > 1 )
		{
			Logger.error( "Only single channel files are supported. " +
					"Please use [ Image > Color > Split Channels ] and [ File > Save as..] to " +
					"save the channel that you want to segment and track as a single file.");
			return;
		}

		intensities = Utils.get2DImagePlusMovieAsFrameList(
				imagePlus,
				1,
				tMinOneBased,
				Math.min( tMaxOneBased, imagePlus.getNFrames() ) );
	}

	private ArrayList< RandomAccessibleInterval< T > > computeLabels(
			ArrayList< RandomAccessibleInterval< T > > intensities,
			double[] calibration )
	{
		final MicrogliaSegmentationAndTracking microgliaSegmentationAndTracking =
				new MicrogliaSegmentationAndTracking(
						intensities,
						calibration,
						showIntermediateResults,
						opService );

		microgliaSegmentationAndTracking.run();

		return (ArrayList< RandomAccessibleInterval< T > > ) microgliaSegmentationAndTracking.getLabelings();
	}

	private void saveLabels( ArrayList< RandomAccessibleInterval< T > > labelings )
	{
		final ImagePlus labelsImp = Utils.labelingsAsImagePlus( labelings );

		final String outputPath = outputDirectory + File.separator + inputIntensitiesFile.getName().split( "\\." )[ 0 ] + "-labelMasks.tif";

		new FileSaver( labelsImp ).saveAsTiff( outputPath );

		Logger.log( "Segmented and tracked label images saved: " + outputPath );
	}

}
