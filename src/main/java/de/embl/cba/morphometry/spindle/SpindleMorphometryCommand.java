package de.embl.cba.morphometry.spindle;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.tables.InteractiveTablePanel;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imagej.table.GenericTable;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;
import java.util.HashMap;
import java.util.Map;


@Plugin(type = Command.class, menuPath = "Plugins>Morphometry>Spindle Morphometry" )
public class SpindleMorphometryCommand< R extends RealType< R > > implements Command
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
		settings.derivativeDelta = 1.0;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
		settings.watershedSeedsGlobalDistanceThreshold = 2.0;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = outputDirectory;
		settings.erosionOfDapiMaskInCalibratedUnits = dapiMaskErosion;
	}

	private void processFile( File file )
	{
		final ImagePlus imagePlus = IJ.openImage( file.toString() );

		setSettingsFromImagePlus( imagePlus );

		final RandomAccessibleInterval< R > rai = ImageJFunctions.wrapReal( imagePlus );

		final RandomAccessibleInterval< R > dapi = Views.hyperSlice( rai, 2, dapiChannelIndexOneBased - 1 );
		final RandomAccessibleInterval< R > tubulin = Views.hyperSlice( rai, 2, tubulinChannelIndexOneBased - 1 );

		settings.dapiImage = dapi;
		settings.tubulinImage = tubulin;

		SpindleMorphometry morphometry = new SpindleMorphometry( settings, opService );
		morphometry.run();

		final HashMap<Integer, Map< String, Object > > objectMeasurements = morphometry.getObjectMeasurements();

		// TODO: get rid of genericTable
		final GenericTable genericTable = Measurements.createGenericTable( objectMeasurements );
		final InteractiveTablePanel interactiveTablePanel = new InteractiveTablePanel( genericTable );

		Utils.log( "Done!" );

	}

	private void setSettingsFromImagePlus( ImagePlus imagePlus )
	{
		settings.inputCalibration = Utils.getCalibration( imagePlus );
		settings.maxPossibleValueInDataSet = Math.pow( 2, imagePlus.getBitDepth() ) - 1.0;
		settings.inputDataSetName = imagePlus.getTitle();
	}

}
