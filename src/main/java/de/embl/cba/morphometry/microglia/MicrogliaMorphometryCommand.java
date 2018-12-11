package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.skeleton.SkeletonCreator;
import de.embl.cba.tables.InteractiveTablePanel;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imagej.table.GenericTable;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static de.embl.cba.morphometry.microglia.Constants.INTENSITIES;
import static de.embl.cba.morphometry.microglia.Constants.SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING;


@Plugin(type = Command.class, menuPath = "Plugins>Tracking>Microglia Morphometry" )
public class MicrogliaMorphometryCommand<T extends RealType<T> & NativeType< T > > implements Command
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

	@Parameter
	public File inputImageFile;

	@Parameter ( style = "save" )
	public File outputTableFile;

	MicrogliaMorphometrySettings settings = new MicrogliaMorphometrySettings<>();

	@Parameter ( label = "Label map channel index", min = "1")
	public long labelMapChannelIndexOneBased = settings.labelMapChannelIndex;

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;

	private ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimepointList;
	private ArrayList< RandomAccessibleInterval< BitType > > skeletons;
	private RandomAccessibleInterval< IntType > labelMaps;

	public void run()
	{
		processFile( inputImageFile );
	}

	private void processFile( File file )
	{
		// TODO: refactor out of command into separate class
		final RandomAccessibleInterval inputImage = openInputImage( file );

		labelMaps = Views.dropSingletonDimensions( Views.hyperSlice( inputImage, Constants.CHANNEL, settings.labelMapChannelIndex ) );

		createSkeletons( );

		initObjectMeasurements( );

		performMeasurements( );

		final ArrayList< String > measurements = Measurements.asTableRows( measurementsTimepointList );

		Measurements.saveMeasurements( outputTableFile, measurements );

		showResults( file, measurements );
	}

	private void showResults( File file, ArrayList< String > measurements )
	{
		final ImagePlus imagePlus = IJ.openImage( file.getAbsolutePath() );
		imagePlus.show();

		final GenericTable table = Measurements.createGenericTableFromTableRows( measurements );

		int[] xyzt = new int[ 4 ];
		xyzt[ 0 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + "X" + Measurements.SEP + Measurements.PIXEL_UNITS );
		xyzt[ 1 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + "Y" + Measurements.SEP + Measurements.PIXEL_UNITS );
		xyzt[ 2 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + "Z" + Measurements.SEP + Measurements.PIXEL_UNITS );
		xyzt[ 3 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + Measurements.TIME + Measurements.SEP +  Measurements.FRAME_UNITS );

		final InteractiveTablePanel interactiveTablePanel = new InteractiveTablePanel( table );
		interactiveTablePanel.setCoordinateColumns( xyzt );
		interactiveTablePanel.setImagePlus( imagePlus );
	}

	private RandomAccessibleInterval openInputImage( File file )
	{
		final ImagePlus imagePlus = ImageIO.openWithBioFormats( file.getAbsolutePath() );

		if ( imagePlus == null )
		{
			Utils.error( "Could not open image: " + file );
		}

		configureSettings( imagePlus );

		return ImageJFunctions.wrapReal( imagePlus );
	}

	private void performMeasurements( )
	{
		for ( int t = 0; t < labelMaps.dimension( 2 ); ++t )
		{
			final HashMap< Integer, Map< String, Object > > measurements = measurementsTimepointList.get( t );

			final ImgLabeling< Integer, IntType > imgLabeling = Utils.labelMapAsImgLabelingRobert( Views.hyperSlice( labelMaps, 2, t ) );

			Measurements.measurePositions(
					measurements,
					imgLabeling,
					null);

			// Volumes ( = areas )
			Measurements.measureVolumes(
					measurements,
					imgLabeling);

			// Surfaces ( = perimeters )
			Measurements.measureSurface(
					measurements,
					imgLabeling,
					opService );

			// TODO: move to skeletonAnalyzer?
			Measurements.measureSkeletons(
					measurements,
					imgLabeling,
					skeletons.get( t ),
					opService );

			// Form factor could be calculated later, e.g. in R

			// Analyze Skeletons: length, branch-points, branches
			// avg branch-length = length / branches

			// Measure: distance travelled

			// Also,we are presently using MtrackJ to calculate velocity, distance travelled and displacement.
			// => I would recommend you do this in Excel as this is downstream analysis.
			// => What I can work on is a tool to upload your extended table again and view it on top of the objects

			// With the segmented microglia movie generated with your plugin, can we do automatic tracking?
			// The cells are already tracked.

			//	2-The next  challenge would be to measure phagocytosis of green particles and quantify "black holes" as we discussed last summer.


		}
	}

	private void createSkeletons( )
	{
		final SkeletonCreator skeletonCreator =
				new SkeletonCreator( Utils.labelMapsAsMasks( labelMaps ), settings );
		skeletonCreator.run();
		skeletons = skeletonCreator.getSkeletons();
	}

	private void initObjectMeasurements( )
	{
		measurementsTimepointList = new ArrayList<>();
		for ( int t = 0; t < labelMaps.dimension( 2 ); ++t )
		{
			final HashMap< Integer, Map< String, Object > > objectMeasurements = new HashMap<>();
			measurementsTimepointList.add( objectMeasurements );
		}
	}

	private void createOutput( ArrayList< RandomAccessibleInterval< T > > intensities, ArrayList< RandomAccessibleInterval< T > > labelings )
	{
		ImagePlus labelImagePlus = Utils.createIJ1Movie( labelings, SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING );
		labelImagePlus.setLut( Utils.getGoldenAngleLUT() );
		labelImagePlus.setTitle( SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING );
		labelImagePlus.show();
		IJ.run( labelImagePlus, "Enhance Contrast", "saturated=0.35");
		IJ.wait( 1000 );

		Utils.createIJ1Movie( intensities, INTENSITIES ).show();
		IJ.wait( 1000 );
		IJ.run("16-bit", "");
		IJ.run("Enhance Contrast", "saturated=0.35");

		IJ.wait( 1000 );
		IJ.run("Merge Channels...", "c1=intensities c2=[" + SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING + "] create keep");
	}


	public void configureSettings( ImagePlus imagePlus )
	{
		settings.inputCalibration = Utils.get2dCalibration( imagePlus ) ;
		settings.opService = opService;
		settings.labelMapChannelIndex = labelMapChannelIndexOneBased - 1;
	}


}
