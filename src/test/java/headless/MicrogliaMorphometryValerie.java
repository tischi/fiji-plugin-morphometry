package headless;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import de.embl.cba.tables.Tables;
import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import javax.swing.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class MicrogliaMorphometryValerie
{

	public static < T extends RealType< T > & NativeType< T > > void main( String[] args )
	{
		final net.imagej.ImageJ ij = new net.imagej.ImageJ();

		final File tableOutputFile = new File( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/microglia/MAX_5C-crop-t1-3-measurements.csv" );

		final String labelMaskInputFile = "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/microglia/MAX_5C-crop-t1-3-labelMasks.tif";

		final ImagePlus imagePlus = ImageIO.openWithBioFormats( labelMaskInputFile );

		final ArrayList< RandomAccessibleInterval< T > > labelMasks =
				Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final MicrogliaMorphometry microgliaMorphometry
				= new MicrogliaMorphometry( labelMasks, ij.op() );

		microgliaMorphometry.run();

		final ArrayList< HashMap< Integer, Map< String, Object > > >
				measurements = microgliaMorphometry.getMeasurementsTimepointList();

		final JTable jTable = Measurements.asTable( measurements );

		Tables.saveTable( jTable,
				tableOutputFile );

	}
}
