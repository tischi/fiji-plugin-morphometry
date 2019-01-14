package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.display.DisplayService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

@Plugin(type = Command.class, menuPath = "Plugins>Measurement>Measure Microglia Morphometry" )
public class MicrogliaMorphometryCommand<T extends RealType<T> & NativeType< T > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter
	public DisplayService displayService;

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

		// TODO: save skeletons!

		final ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimepointList = microgliaMorphometry.getMeasurementsTimepointList();

		final ArrayList< String > measurementsAsRows = Measurements.measurementsAsTableRowsStringList( measurementsTimepointList, "\t" );

		Measurements.saveRowsToFile( new File( outputDirectory.toString() + File.separator + imagePlus.getTitle() + ".csv"), measurementsAsRows );

	}
}
