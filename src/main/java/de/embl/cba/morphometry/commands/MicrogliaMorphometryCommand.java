package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import de.embl.cba.tables.Tables;
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

	@Parameter ( label = "Label masks" )
	public File labelMaskFile;

	@Parameter ( label = "Raw image data" )
	private File intensityFile;

	@Parameter ( style = "directory" )
	public File outputDirectory;

	@Parameter
	public boolean showIntermediateResults = false;


	private ImagePlus labelMaskImagePlus;
	private File tableOutputFile;
	private String dataSetID;
	private ImagePlus intensityImagePlus;


	public void run()
	{
		Logger.log( "" );
		Logger.log( "" );
		Logger.log( "Microglia Morphometry Command" );
		Logger.log( "Analyzing: " + labelMaskFile );

//		fetchFilesFromFolder();

		intensityFile = new File( labelMaskFile.toString().replace( "-labelMasks",  "" ) );

		final MicrogliaMorphometry< T > microgliaMorphometry =
				new MicrogliaMorphometry(
						openLabelMasks(),
						openIntensities(),
						opService );

		microgliaMorphometry.run();

		saveResults( dataSetID, microgliaMorphometry );

		Logger.log( "Done!" );

	}

	public ArrayList< RandomAccessibleInterval< T > > openLabelMasks()
	{
		labelMaskImagePlus = Utils.openWithBioFormats( labelMaskFile.toString() );
		dataSetID = labelMaskImagePlus.getTitle();
		return Utils.get2DImagePlusMovieAsFrameList( labelMaskImagePlus, 1 );
	}

	public ArrayList< RandomAccessibleInterval< T > > openIntensities()
	{
		intensityImagePlus = Utils.openWithBioFormats( intensityFile.toString() );
		return Utils.get2DImagePlusMovieAsFrameList( intensityImagePlus, 1 );
	}

	private void saveResults( String dataSetID,
							  MicrogliaMorphometry< T > microgliaMorphometry )
	{
		final ArrayList< HashMap< Integer, Map< String, Object > > >
				measurementsTimepointList = microgliaMorphometry.getMeasurementsTimepointList();

		Measurements.addCalibration( measurementsTimepointList, labelMaskImagePlus );

		final JTable table = Measurements.asTable( measurementsTimepointList );

		tableOutputFile = new File(
				outputDirectory.toString() + File.separator + dataSetID + ".csv" );

		addRelativeImagePathColumn( table,
				outputDirectory, labelMaskFile, "LabelMasks" );

		addRelativeImagePathColumn( table,
				outputDirectory, intensityFile, "Intensities" );

		saveSkeletons( dataSetID, microgliaMorphometry, table );

		Logger.log( "Saving results table: " + tableOutputFile );

		saveTable( table, tableOutputFile );
	}

	public File getTableOutputFile()
	{
		return tableOutputFile;
	}

	private void saveSkeletons(
			String dataSetID,
			MicrogliaMorphometry< T > microgliaMorphometry,
			JTable table )
	{
		final File file = new File( dataSetID + "-skeletons.tif" );

		final String imageName = dataSetID + "-skeletons";

		Utils.saveRAIListAsMovie(
				microgliaMorphometry.getSkeletons(),
				labelMaskImagePlus.getCalibration(),
				outputDirectory.toString() + File.separator + file,
				imageName );

		Tables.addColumn( table, "Path_Skeletons", file );
	}

	private void saveAnnotations(
			String dataSetID,
			MicrogliaMorphometry< T > microgliaMorphometry,
			JTable table )
	{
		final File file = new File( dataSetID + "-annotations.tif" );

		final String imageName = dataSetID + "-annotations";

		Utils.saveRAIListAsMovie(
				microgliaMorphometry.getAnnotations(),
				labelMaskImagePlus.getCalibration(),
				outputDirectory.toString() + File.separator + file,
				imageName );

		Tables.addColumn( table, "Path_Annotations", file );
	}

}
