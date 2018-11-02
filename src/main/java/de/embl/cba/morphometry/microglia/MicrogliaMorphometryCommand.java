package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.ObjectMeasurements;
import de.embl.cba.morphometry.skeleton.Skeleton;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
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
import java.io.FileNotFoundException;
import java.io.PrintWriter;
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
		final RandomAccessibleInterval inputImage = openInputImage( file );

		labelMaps = Views.dropSingletonDimensions( Views.hyperSlice( inputImage, Constants.CHANNEL, settings.labelMapChannelIndex ) );

		createSkeletons( );

		initObjectMeasurements( );

		performMeasurements( );

		final ArrayList< String > lines = ObjectMeasurements.printMeasurements( measurementsTimepointList );

		logMeasurements( lines );

		saveMeasurements( outputTableFile, lines );
	}

	public void saveMeasurements( File file, ArrayList<String> lines )
	{
		try (PrintWriter out = new PrintWriter( file ) )
		{
			for ( String line : lines )
			{
				out.println( line );
			}

			Utils.log( "\nSaved table to: " + file );
		}
		catch ( FileNotFoundException e )
		{
			e.printStackTrace();
		}
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

	public static void logMeasurements( ArrayList< String > lines )
	{
		Utils.log( " ");
		Utils.log( "----------- RESULTS (Tab delimited => Copy to Excel) -------------");
		Utils.log( " ");

		for ( String line : lines )
		{
			Utils.log( line );
		}
	}


	private void performMeasurements( )
	{
		for ( int t = 0; t < labelMaps.dimension( 2 ); ++t )
		{
			final HashMap< Integer, Map< String, Object > > measurements = measurementsTimepointList.get( t );

			final ImgLabeling< Integer, IntType > imgLabeling = Utils.labelMapAsImgLabeling( Views.hyperSlice( labelMaps, 2, t ) );

			ObjectMeasurements.measurePositions(
					measurements,
					imgLabeling,
					null);

			ObjectMeasurements.measureSizes(
					measurements,
					imgLabeling );

			ObjectMeasurements.measureSumIntensities(
					measurements,
					imgLabeling,
					skeletons.get( t ),
					"Skeleton" );


		}
	}

	private void createSkeletons( )
	{
		final Skeleton skeleton = new Skeleton( labelMaps, settings );
		skeleton.run();
		skeletons = skeleton.getSkeletons();

		//ImagePlus skeletonImagePlus = Utils.createIJ1Movie( skeletons, "Skeletons" );
		//skeletonImagePlus.show();

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
