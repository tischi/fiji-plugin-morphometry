package de.embl.cba.morphometry.commands;

import com.jgoodies.forms.layout.FormLayout;
import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import de.embl.cba.tables.TableUtils;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.display.DisplayService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static de.embl.cba.tables.TableUtils.addRelativeImagePathColumn;
import static de.embl.cba.tables.TableUtils.saveTable;

@Plugin(type = Command.class, menuPath = "Plugins>Morphometry>Measure Microglia Morphometry" )
public class MicrogliaMorphometryCommand<T extends RealType<T> & NativeType< T > > implements Command
{
	public static final String DELIM = "\t";
	@Parameter
	public OpService opService;

	@Parameter
	public DisplayService displayService;

	@Parameter ( label = "File containing object label masks" )
	public File inputLabelMaskFile;

	@Parameter ( label = "File containing intensities (optional)" )
	public File inputIntensityFile;

	@Parameter ( style = "directory" )
	public File outputDirectory;

	@Parameter
	public boolean showIntermediateResults = false;


	public void run()
	{

		final ImagePlus imagePlus = ImageIO.openWithBioFormats( inputLabelMaskFile.toString() );

		final String dataSetID = imagePlus.getTitle();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks = Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final MicrogliaMorphometry microgliaMorphometry = new MicrogliaMorphometry( labelMasks, opService );

		microgliaMorphometry.run();

		saveMeasurements( dataSetID, microgliaMorphometry );

		// TODO: save skeletons

	}

	public void saveMeasurements( String dataSetID, MicrogliaMorphometry microgliaMorphometry )
	{
		final JTable table = Measurements.asTable( microgliaMorphometry.getMeasurementsTimepointList() );

		final File tableOutputFile = new File( outputDirectory.toString() + File.separator + dataSetID + ".csv" );

		addRelativeImagePathColumn( table, tableOutputFile, inputLabelMaskFile, "LabelMasks" );
		addRelativeImagePathColumn( table, tableOutputFile, inputIntensityFile, "Intensities" );

		saveTable( table, tableOutputFile );
	}

}
