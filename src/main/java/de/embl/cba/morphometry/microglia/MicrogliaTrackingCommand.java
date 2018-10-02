package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.segmentation.SimpleSegmenter;
import de.embl.cba.morphometry.splitting.TrackingSplitter;
import de.embl.cba.morphometry.tracking.MaximalOverlapTracker;
import ij.IJ;
import ij.ImagePlus;
import javafx.scene.control.RadioMenuItem;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;
import java.util.ArrayList;

import static de.embl.cba.morphometry.microglia.Constants.INTENSITIES;
import static de.embl.cba.morphometry.microglia.Constants.SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING;


@Plugin(type = Command.class, menuPath = "Plugins>Tracking>Microglia Tracking" )
public class MicrogliaTrackingCommand<T extends RealType<T> & NativeType< T > > implements Command
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


	MicrogliaSettings settings = new MicrogliaSettings();


	public static final String PROCESS_DIRECTORY = "Process all files in directory";
	public static final String PROCESS_CURRENT_IMAGE = "Current image";

	@Parameter( choices = { PROCESS_DIRECTORY }) // , PROCESS_CURRENT_IMAGE
	public String inputModality = PROCESS_DIRECTORY;

	@Parameter( style = "directory" )
	public File inputDirectory;

	@Parameter
	public String fileNameEndsWith = ".lif";

	@Parameter ( label = "Microglia channel index", min = "1")
	public long microgliaChannelIndexOneBased = settings.microgliaChannelIndexOneBased;

	@Parameter ( label = "Minimal time frame to be processed", min = "1" )
	public long tMin = 0;

	@Parameter ( label = "Maximal time frame to be processed", min = "1" )
	public long tMax = 100;

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;


	public void run()
	{
		if ( inputModality.equals( PROCESS_CURRENT_IMAGE ) && imagePlus != null )
		{
			//
		}
		else if ( inputModality.equals( PROCESS_DIRECTORY ) )
		{
			processDirectory();
		}

		Utils.log( "Done!" );

	}

	public void processDirectory()
	{
		final File directory = inputDirectory;

		String[] files = directory.list();

		for( String file : files )
		{
			if ( Utils.acceptFile( fileNameEndsWith, file ) )
			{
				processFile( new File ( directory + File.separator + file ) );
			}
		}
	}

	private void processFile( File file )
	{
		final ImagePlus imagePlus = ImageIO.openWithBioFormats( file.getAbsolutePath() );

		if ( imagePlus == null )
		{
			Utils.error( "Could not open image: " + file );
		}

		configureSettings( imagePlus );

		ArrayList< RandomAccessibleInterval< T > > intensities =
				Algorithms.createMaximumProjectedIntensitiesAssumingImagePlusDimensionOrder(
						( RandomAccessibleInterval ) ImageJFunctions.wrapReal( imagePlus ),
						microgliaChannelIndexOneBased - 1,
						tMin, tMax );

		ArrayList<  RandomAccessibleInterval< T > > masks = new ArrayList<>();
		for ( long t = 0; t <= ( tMax - tMin ) ; ++t )
		{
			final SimpleSegmenter simpleSegmenter = new SimpleSegmenter( intensities.get( ( int ) t ), settings );
			simpleSegmenter.run();
			masks.add( simpleSegmenter.getMask() );
		}

		// ImageJFunctions.show( Views.stack( masks ), "Simple segmentation" );

		final TrackingSplitter splitter = new TrackingSplitter( masks, intensities, settings );
		splitter.run();
		masks = splitter.getSplitMasks();

		final MaximalOverlapTracker maximalOverlapTracker = new MaximalOverlapTracker( masks );
		maximalOverlapTracker.run();
		final ArrayList< RandomAccessibleInterval< T > > labelings = maximalOverlapTracker.getLabelings();

		createOutput( intensities, labelings );

	}

	private void createOutput( ArrayList< RandomAccessibleInterval< T > > intensities, ArrayList< RandomAccessibleInterval< T > > labelings )
	{
		Utils.showAsIJ1Movie( labelings, SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING );
		IJ.run("3-3-2 RGB", "");

		Utils.showAsIJ1Movie( intensities, INTENSITIES );
		IJ.wait( 500 );
		IJ.run("16-bit", "");

		IJ.wait( 500 );
		IJ.run("Merge Channels...", "c1=intensities c2=[Simple segmentation - Tracking splitting - Simple tracking] create");
	}


	public void configureSettings( ImagePlus imagePlus )
	{
		settings.inputCalibration = Utils.get2dCalibration( imagePlus ) ;
		settings.workingVoxelSize = settings.inputCalibration[ 0 ];
		settings.maxPossibleValueInDataSet = Math.pow( 2, imagePlus.getBitDepth() ) - 1.0;
		settings.maxShortAxisDist = 6;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = Double.MAX_VALUE;
		settings.watershedSeedsGlobalDistanceThreshold = 2.5;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = null; //new File( path ).getParentFile();
		settings.inputDataSetName = "test";
		settings.returnEarly = true;
		settings.skeletonMaxLength = 600 * settings.workingVoxelSize;
		settings.minimalObjectSize = 200;  // um2
		settings.minimalObjectCenterDistance = 6;
		settings.maximalWatershedLength = 10;
		settings.closingRadius = 3;
		settings.showIntermediateResults = false;
		settings.opService = opService;
		settings.microgliaChannelIndexOneBased = microgliaChannelIndexOneBased;
		settings.tMin = tMin;
		tMax = Math.min( tMax, imagePlus.getNFrames() );
		settings.tMax = tMax;
	}


}
