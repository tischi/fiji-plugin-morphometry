package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.spindle.SpindleMorphometry;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import de.embl.cba.tables.TableUtils;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.util.HashMap;
import java.util.Map;


@Plugin(type = Command.class, menuPath = "Plugins>Morphometry>Spindle Morphometry" )
public class SpindleMorphometryCommand< R extends RealType< R > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter
	public File inputImageFile;

	@Parameter ( style = "directory" )
	public File outputDirectory;

	@Parameter ( label = "DNA channel [one-based index]" )
	public long dnaChannelIndexOneBased = 2;

	@Parameter ( label = "Spindle channel [one-based index]" )
	public long spindleChannelIndexOneBased = 1;

	@Parameter
	public boolean showIntermediateResults = false;

	private String imageName;
	private SpindleMorphometrySettings settings = new SpindleMorphometrySettings();

	public void run()
	{
		setSettingsFromUI();
		processFile( inputImageFile );
	}

	private void setSettingsFromUI()
	{
		settings.showIntermediateResults = showIntermediateResults;
		settings.workingVoxelSize = 0.25;
		settings.maxShortAxisDist = 6;
		settings.derivativeDelta = 3.0;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
		settings.watershedSeedsGlobalDistanceThreshold = 2.0;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = outputDirectory;
	}

	private void processFile( File file )
	{
		logStart();

		imageName = inputImageFile.getName().replace( ".tif", "" );

		final ImagePlus imagePlus = IJ.openImage( file.toString() );
		setSettingsFromImagePlus( imagePlus );

		final RandomAccessibleInterval< R > rai = ImageJFunctions.wrapReal( imagePlus );
		final RandomAccessibleInterval< R > dapi = Views.hyperSlice( rai, 2, dnaChannelIndexOneBased - 1 );
		final RandomAccessibleInterval< R > tubulin = Views.hyperSlice( rai, 2, spindleChannelIndexOneBased - 1 );

		settings.dnaImage = dapi;
		settings.tubulinImage = tubulin;

		SpindleMorphometry morphometry = new SpindleMorphometry( settings, opService );
		morphometry.run();

		getAndSaveOutputImage( morphometry );

		computeAndSaveMeasurements( morphometry );

		logEnd();

	}

	private void getAndSaveOutputImage( SpindleMorphometry morphometry )
	{
		save( morphometry.getOutputImage() );
	}

	private void logEnd()
	{
		Logger.log( "Done!" );
	}

	private void logStart()
	{
		Logger.log( " " );
		Logger.log( "# Spindle Morphometry" );
		Logger.log( "Processing file " + imageName );
	}

	private void computeAndSaveMeasurements( SpindleMorphometry morphometry )
	{
		final HashMap<Integer, Map< String, Object > > objectMeasurements = morphometry.getObjectMeasurements();
		final JTable jTable = Measurements.asTable( objectMeasurements );


		final File tableOutputFile = new File( getOutputDirectory() + "measurements.txt" );

		Logger.log( "Saving " + tableOutputFile );
		TableUtils.saveTable( jTable, tableOutputFile );
	}

	private String getOutputDirectory()
	{
		return outputDirectory
				+ File.separator
				+ imageName
				+ File.separator;
	}

	private void setSettingsFromImagePlus( ImagePlus imagePlus )
	{
		settings.inputCalibration = Utils.getCalibration( imagePlus );
		settings.imagePlusCalibration = imagePlus.getCalibration();
		settings.inputDataSetName = imagePlus.getTitle();
	}

	private void save( ImagePlus imagePlus )
	{
		final String path = getOutputDirectory() + imageName + "-out.tif";
		new File( path ).getParentFile().mkdirs();
		Logger.log( "Saving: " + path );
		IJ.saveAsTiff( imagePlus, path );
	}

}
