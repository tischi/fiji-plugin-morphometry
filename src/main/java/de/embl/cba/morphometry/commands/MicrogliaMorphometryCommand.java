package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import de.embl.cba.tables.FileUtils;
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
import java.util.ArrayList;
import java.util.List;

import static de.embl.cba.tables.TableUtils.addRelativeImagePathColumn;
import static de.embl.cba.tables.TableUtils.saveTable;

@Plugin(type = Command.class, menuPath = "Plugins>Morphometry>Measure Microglia Morphometry" )
public class MicrogliaMorphometryCommand<T extends RealType<T> & NativeType< T > > implements Command
{
	public static final String DELIM = "\t";

	@Parameter
	public OpService opService;

//	@Parameter ( style = "Input directory" )
	public File inputDirectory;
//
	@Parameter ( label = "File containing object label masks" )
	public File labelMaskFile;

//	@Parameter ( label = "File containing intensities (optional)" )
	public File intensityFile;

	@Parameter ( style = "directory" )
	public File outputDirectory;

	@Parameter
	public boolean showIntermediateResults = false;


	public void run()
	{
		Logger.log( "MicrogliaMorphometryCommand: Running..." );

//		fetchFilesFromFolder();

		intensityFile = new File( labelMaskFile.toString().replace( "-labelMasks",  "" ) );

		final ImagePlus imagePlus = ImageIO.openWithBioFormats( labelMaskFile.toString() );

		final String dataSetID = imagePlus.getTitle();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks = Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final MicrogliaMorphometry microgliaMorphometry = new MicrogliaMorphometry( labelMasks, opService );

		microgliaMorphometry.run();

		saveMeasurements( dataSetID, microgliaMorphometry );

		Logger.log( "MicrogliaMorphometryCommand: Done!" );

		// TODO: save skeletons

	}

	public void fetchFilesFromFolder()
	{
		final List< File > fileList = FileUtils.getFileList( inputDirectory, ".*\\.tif", false );

		for ( File file : fileList )
		{
			if ( file.toString().contains( "-labelMask" ) )
			{
				labelMaskFile = file;
			}
			else
			{
				intensityFile = file;
			}
		}
	}

	public void saveMeasurements( String dataSetID, MicrogliaMorphometry microgliaMorphometry )
	{
		final JTable table = Measurements.asTable( microgliaMorphometry.getMeasurementsTimepointList() );

		final File tableOutputFile = new File( outputDirectory.toString() + File.separator + dataSetID + ".csv" );

		addRelativeImagePathColumn( table, outputDirectory, labelMaskFile, "LabelMasks" );

		addRelativeImagePathColumn( table, outputDirectory, intensityFile, "Intensities" );

		Logger.log( "Saving results table: " + tableOutputFile );
		saveTable( table, tableOutputFile );
	}

}
