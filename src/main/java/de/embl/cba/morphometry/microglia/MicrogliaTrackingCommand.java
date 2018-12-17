package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.segmentation.SimpleSegmenter;
import de.embl.cba.morphometry.splitting.TrackingSplitter;
import de.embl.cba.morphometry.tracking.MaximalOverlapTracker;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
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

	MicrogliaTrackingSettings settings = new MicrogliaTrackingSettings();

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
	public long tMin = 1;

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
		ImagePlus imagePlus = ImageIO.openWithBioFormats( file.getAbsolutePath() );

		if ( imagePlus == null )
		{
			Utils.error( "Could not open image: " + file );
		}

		configureSettings( imagePlus );

		ArrayList< RandomAccessibleInterval< T > > intensities = createMaximumProjection( imagePlus );

		imagePlus = null; System.gc();

		ArrayList< RandomAccessibleInterval< T > > masks = createBinaryMasks( intensities );

		masks = splitTouchingObjects( intensities, masks );

		final ArrayList< RandomAccessibleInterval< T > > labelings = createTrackingBasedLabels( masks );

		createOutput( intensities, labelings );

	}

	private ArrayList< RandomAccessibleInterval< T > > createTrackingBasedLabels( ArrayList< RandomAccessibleInterval< T > > masks )
	{
		final MaximalOverlapTracker maximalOverlapTracker = new MaximalOverlapTracker( masks );
		maximalOverlapTracker.run();
		return (ArrayList< RandomAccessibleInterval< T > > ) maximalOverlapTracker.getLabelings();
	}

	private ArrayList< RandomAccessibleInterval< T > > splitTouchingObjects( ArrayList< RandomAccessibleInterval< T > > intensities, ArrayList< RandomAccessibleInterval< T > > masks )
	{
		final TrackingSplitter splitter = new TrackingSplitter( masks, intensities, settings );
		splitter.run();
		masks = splitter.getSplitMasks();
		return masks;
	}

	private ArrayList< RandomAccessibleInterval< T > > createBinaryMasks( ArrayList< RandomAccessibleInterval< T > > intensities )
	{
		ArrayList<  RandomAccessibleInterval< T > > masks = new ArrayList<>();
		for ( long t = 0; t < intensities.size() ; ++t )
		{
			Utils.log("Creating mask for frame " + ( t + 1 ) );
			final SimpleSegmenter simpleSegmenter = new SimpleSegmenter( intensities.get( ( int ) t ), settings );
			simpleSegmenter.run();
			masks.add( simpleSegmenter.getMask() );
		}
		return masks;
	}


	private ArrayList< RandomAccessibleInterval< T > > createMaximumProjection( ImagePlus imagePlus )
	{
		ArrayList< RandomAccessibleInterval< T > > intensities =
				Algorithms.createMaximumIntensityProjectionsAssumingImagePlusDimensionOrder(
						( RandomAccessibleInterval ) ImageJFunctions.wrapReal( imagePlus ),
						microgliaChannelIndexOneBased - 1,
						settings.tMin, settings.tMax );


		return intensities;
	}

	private void createOutput( ArrayList< RandomAccessibleInterval< T > > intensities, ArrayList< RandomAccessibleInterval< T > > labelings )
	{
		ImagePlus labelImagePlus = Utils.createIJ1Movie( labelings, SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING );
		labelImagePlus.setLut( Utils.getGoldenAngleLUT() );
		labelImagePlus.setTitle( SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING );
		labelImagePlus.show();
		IJ.run( labelImagePlus, "Enhance Contrast", "saturated=0.35");
		IJ.wait( 1000 );
		labelings = null; System.gc();

		Utils.createIJ1Movie( intensities, INTENSITIES ).show();
		IJ.wait( 1000 );
		IJ.run("16-bit", "");
		IJ.run("Enhance Contrast", "saturated=0.35");
		intensities = null; System.gc();

		IJ.wait( 1000 );
		IJ.run("Merge Channels...", "c1=intensities c2=[" + SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING + "] create keep");
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
		settings.minimalTrackingSplittingObjectArea = 20; // this can be very small, um2
		settings.minimalObjectCenterDistance = 6;
		settings.maximalWatershedLength = 10;
		settings.closingRadius = 3;
		settings.showIntermediateResults = showIntermediateResults;
		settings.opService = opService;
		settings.microgliaChannelIndexOneBased = microgliaChannelIndexOneBased;
		settings.tMin = tMin - 1;
		tMax = Math.min( tMax, imagePlus.getNFrames() );
		settings.tMax = tMax - 1;
	}


}
