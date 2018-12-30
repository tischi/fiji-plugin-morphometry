package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import ij.ImagePlus;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
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
import java.util.HashMap;
import java.util.Map;


@Plugin(type = Command.class, menuPath = "Plugins>Tracking>Microglia Morphometry" )
public class MicrogliaMorphometryCommand<T extends RealType<T> & NativeType< T > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter
	public File inputLabelMaskFile;

	@Parameter ( style = "directory" )
	public File outputDirectory;

	@Parameter
	public boolean showIntermediateResults = false;


	public void run()
	{

		final ImagePlus imagePlus = ImageIO.openWithBioFormats( inputLabelMaskFile.toString() );

		final ArrayList< RandomAccessibleInterval< T > > labelMasks = Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final MicrogliaMorphometry microgliaMorphometry = new MicrogliaMorphometry( labelMasks, opService );

		microgliaMorphometry.run();

		final ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimepointList = microgliaMorphometry.getMeasurementsTimepointList();

		final ArrayList< String > measurementsAsRows = Measurements.asTableRows( measurementsTimepointList, "\t" );

		Measurements.saveRowsToFile( outputDirectory, measurementsAsRows );

	}
}
