import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaTracking;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import net.imagej.ops.DefaultOpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class TestMicrogliaTracking
{
	public static < T extends RealType< T > & NativeType< T > > void main ( String[] args )
	{
		final net.imagej.ImageJ ij = new net.imagej.ImageJ();

		final String inputFile = TestMicrogliaTracking.class.getResource( "microglia/MAX_5C-crop-t1-3.tif" ).getFile().toString();

		ImagePlus imagePlus = ImageIO.openWithBioFormats( inputFile );

		final MicrogliaTracking microgliaTracking = new MicrogliaTracking(
				imagePlus,
				false,
				ij.op(),
				2, 1, 100000 );

		microgliaTracking.run();

		final ArrayList< RandomAccessibleInterval< T > > labelings = microgliaTracking.getLabelings();

		final ImagePlus labelMask = Utils.labelingsAsImagePlus( labelings );

		labelMask.show();

		new FileSaver( labelMask ).saveAsTiff( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/microglia/MAX_5C-crop-t1-3-labelMasks.tif" );

	}
}
