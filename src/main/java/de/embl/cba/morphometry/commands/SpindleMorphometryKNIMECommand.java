package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.spindle.SpindleMorphometry;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import net.imagej.Dataset;
import net.imagej.ops.OpService;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.util.HashMap;
import java.util.Map;


@Plugin( type = Command.class, menuPath = "DataSet Spindle Morphometry", headless = true )
public class SpindleMorphometryKNIMECommand< R extends RealType< R > & NativeType< R > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter
	public LogService logService;

	@Parameter( label = "Input image" )
	public Dataset inputImage;

	@Parameter ( label = "DNA channel [zero-based index]" )
	public long dnaChannelIndexOneBased = 1;

	@Parameter ( label = "Spindle channel [zero-based index]" )
	public long spindleChannelIndexOneBased = 0;

	@Parameter ( label = "Spindle length", type = ItemIO.OUTPUT )
	public double output01;

	@Parameter ( label = "Debug image", type = ItemIO.OUTPUT )
	public Img< R > outputImage;

	public boolean showIntermediateResults = false;
	private SpindleMorphometrySettings settings = new SpindleMorphometrySettings();
	private String imageName;

	public void run()
	{
		setSettings();
		imageName = "myImage";
		process( inputImage );
	}

	private void setSettings()
	{
		settings.showIntermediateResults = false;
		settings.outputDirectory = new File( "" );
		settings.workingVoxelSize = 0.25;
		settings.maxShortAxisDist = 6;
		settings.derivativeDelta = 3.0;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
		settings.watershedSeedsGlobalDistanceThreshold = 2.0;
		settings.interestPointsRadius = 0.5;
	}

	private void process( Dataset image )
	{
		logStart();

		setImageSpecificSettings( image );

		settings.dnaImage = Views.hyperSlice( image, 2, dnaChannelIndexOneBased - 1 );
		settings.tubulinImage = Views.hyperSlice( image, 2, spindleChannelIndexOneBased - 1 );

		SpindleMorphometry morphometry = new SpindleMorphometry( settings, opService );
		morphometry.run();

		outputImage = (Img) morphometry.getRandomAccessibleIntervalOutput();

		getMeasurements( morphometry );

		logEnd();
	}


	private void logEnd()
	{
		Logger.log( "Done!" );
	}

	private void logStart()
	{
		logService.info( "# Spindle Morphometry" );

		logService.info( "Processing file " + imageName );
		Logger.log( " " );
		Logger.log( "# Spindle Morphometry" );
		Logger.log( "Processing file " + imageName );
	}

	private void getMeasurements( SpindleMorphometry morphometry )
	{
		final HashMap<Integer, Map< String, Object > > objectMeasurements = morphometry.getObjectMeasurements();
		final JTable jTable = Measurements.asTable( objectMeasurements );

		output01 = (Double) ( jTable.getModel().getValueAt( 0, 2 ) );
	}

	private void setImageSpecificSettings( Dataset image )
	{
		final double[] calibration = new double[ image.numDimensions() ];
		for ( int d = 0; d < 3; d++ )
		{
			calibration[ d ] = 1.0;
		}
		settings.inputCalibration = calibration;
		settings.inputDataSetName = imageName;
	}

}
