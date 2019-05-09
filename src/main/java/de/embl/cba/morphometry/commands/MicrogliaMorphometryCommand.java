package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.util.*;

import static de.embl.cba.tables.Tables.addRelativeImagePathColumn;
import static de.embl.cba.tables.Tables.saveTable;

@Plugin(type = Command.class, menuPath = "Plugins>Morphometry>Microglia Morphometry" )
public class MicrogliaMorphometryCommand < T extends RealType< T > & NativeType< T > >
		implements Command
{
	public static final String DELIM = "\t";

	@Parameter
	public OpService opService;

//	@Parameter ( style = "Input directory" )
	public File inputDirectory;

	@Parameter ( label = "File containing object label masks" )
	public File labelMaskFile;

//	@Parameter ( label = "File containing intensities (optional)" )
	public File intensityFile;

	@Parameter ( style = "directory" )
	public File outputDirectory;

	@Parameter
	public boolean showIntermediateResults = false;
	private ImagePlus imagePlus;


	public void run()
	{
		Logger.log( "MicrogliaMorphometryCommand: Running..." );

//		fetchFilesFromFolder();

		intensityFile = new File( labelMaskFile.toString().replace( "-labelMasks",  "" ) );

		imagePlus = ImageIO.openWithBioFormats( labelMaskFile.toString() );

		final String dataSetID = imagePlus.getTitle();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks =
				Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final MicrogliaMorphometry microgliaMorphometry = new MicrogliaMorphometry(
				labelMasks, opService );

		microgliaMorphometry.run();

		saveMeasurements( dataSetID, microgliaMorphometry );

		Logger.log( "MicrogliaMorphometryCommand: Done!" );

		// TODO: save skeletons

	}

	public void saveMeasurements( String dataSetID,
								  MicrogliaMorphometry< T > microgliaMorphometry )
	{
		final ArrayList< HashMap< Integer, Map< String, Object > > >
				measurementsTimepointList = microgliaMorphometry.getMeasurementsTimepointList();

		Measurements.addCalibration( measurementsTimepointList, imagePlus );

		final JTable table = Measurements.asTable( measurementsTimepointList );

		final File tableOutputFile = new File(
				outputDirectory.toString() + File.separator + dataSetID + ".csv" );

		addRelativeImagePathColumn( table,
				outputDirectory, labelMaskFile, "LabelMasks" );

		addRelativeImagePathColumn( table,
				outputDirectory, intensityFile, "Intensities" );

		Logger.log( "Saving results table: " + tableOutputFile );
		saveTable( table, tableOutputFile );
	}

}
