import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaSegmentationAndTracking;
import ij.ImagePlus;
import ij.io.FileSaver;
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

		final MicrogliaSegmentationAndTracking microgliaSegmentationAndTracking =
				new MicrogliaSegmentationAndTracking(
						Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 ),
						Utils.getCalibration( imagePlus ),
				false, ij.op() );

		microgliaSegmentationAndTracking.run();

		final ArrayList< RandomAccessibleInterval< T > > labelings = microgliaSegmentationAndTracking.getLabelings();

		final ImagePlus labelMask = Utils.labelingsAsImagePlus( labelings );

		labelMask.show();

		new FileSaver( labelMask ).saveAsTiff( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/microglia/MAX_5C-crop-t1-3-labelMasks.tif" );

	}
}
