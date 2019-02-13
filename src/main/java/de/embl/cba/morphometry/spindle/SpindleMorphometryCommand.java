package de.embl.cba.morphometry.spindle;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.tables.TableUtils;
import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.process.LUT;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.awt.*;
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

	SpindleMorphometrySettings settings = new SpindleMorphometrySettings();

	@Parameter
	double dapiMaskErosion = settings.erosionOfDapiMaskInCalibratedUnits;

	@Parameter
	public long dapiChannelIndexOneBased = 2;

	@Parameter
	public long tubulinChannelIndexOneBased = 1;

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;
	private String imageTitle;


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
		settings.erosionOfDapiMaskInCalibratedUnits = dapiMaskErosion;
	}

	private void processFile( File file )
	{
		imageTitle = inputImageFile.getName().replace( ".tif", "" );

		logStart();

		final ImagePlus imagePlus = IJ.openImage( file.toString() );

		setSettingsFromImagePlus( imagePlus );

		final RandomAccessibleInterval< R > rai = ImageJFunctions.wrapReal( imagePlus );

		final RandomAccessibleInterval< R > dapi = Views.hyperSlice( rai, 2, dapiChannelIndexOneBased - 1 );
		final RandomAccessibleInterval< R > tubulin = Views.hyperSlice( rai, 2, tubulinChannelIndexOneBased - 1 );

		settings.dapiImage = dapi;
		settings.tubulinImage = tubulin;

		SpindleMorphometry morphometry = new SpindleMorphometry( settings, opService );
		morphometry.run();

		getAndSaveOutputImage( morphometry );

		computeAndSaveMeasurements( morphometry );

		logEnd();

	}

	private void getAndSaveOutputImage( SpindleMorphometry morphometry )
	{
		final ImagePlus imp = morphometry.getOutputImage();

		final CompositeImage compositeImage = new CompositeImage( imp );

		compositeImage.setC(1);
		IJ.run(compositeImage, "Blue", "");

		compositeImage.setC(2);
		IJ.run(compositeImage, "Green", "");

		compositeImage.setC(3);
		IJ.run(compositeImage, "Yellow", "");

		compositeImage.setDisplayMode( CompositeImage.COMPOSITE );
		save( compositeImage );
	}

	private void logEnd()
	{
		Logger.log( "Done!" );
	}

	private void logStart()
	{
		Logger.log( " " );
		Logger.log( "# Spindle Morphometry" );
		Logger.log( "Processing file " + imageTitle );
	}

	private void computeAndSaveMeasurements( SpindleMorphometry morphometry )
	{
		final HashMap<Integer, Map< String, Object > > objectMeasurements = morphometry.getObjectMeasurements();
		final JTable jTable = Measurements.asTable( objectMeasurements );


		final File tableOutputFile = new File( getOutputDirectory() + "measurements.txt" );

		TableUtils.saveTable( jTable, tableOutputFile );
	}

	private String getOutputDirectory()
	{
		return outputDirectory
				+ File.separator
				+ imageTitle
				+ File.separator;
	}

	private void setSettingsFromImagePlus( ImagePlus imagePlus )
	{
		settings.inputCalibration = Utils.getCalibration( imagePlus );
		settings.imagePlusCalibration = imagePlus.getCalibration();
		settings.maxPossibleValueInDataSet = Math.pow( 2, imagePlus.getBitDepth() ) - 1.0;
		settings.inputDataSetName = imagePlus.getTitle();
	}

	private void save( ImagePlus imagePlus )
	{
		final String path = getOutputDirectory() + imageTitle + "-out.tif";
		new File( path ).getParentFile().mkdirs();
		Logger.log( "Saving: " + path );
		IJ.saveAsTiff( imagePlus, path );
	}

}
