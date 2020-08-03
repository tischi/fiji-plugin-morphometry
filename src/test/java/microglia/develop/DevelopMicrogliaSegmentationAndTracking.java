package microglia.develop;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaSegmentationAndTracking;
import de.embl.cba.morphometry.microglia.MicrogliaSettings;
import ij.ImagePlus;
import ij.io.FileSaver;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;
import java.util.ArrayList;

public class DevelopMicrogliaSegmentationAndTracking
{
	public static < T extends RealType< T > & NativeType< T > > void main ( String[] args )
	{
		final net.imagej.ImageJ ij = new net.imagej.ImageJ();
		ij.ui().showUI();

		final String inputFile = "/Users/tischer/Downloads/Christian_threshold3/pg20-3C12h/MAX_pg20_3C12h.tif";

		ImagePlus imagePlus = Utils.openWithBioFormats( inputFile );

		final MicrogliaSettings settings = new MicrogliaSettings();
		settings.opService = ij.op();
		settings.calibration = imagePlus.getCalibration();
		settings.outputDirectory = new File( "" );
		settings.showIntermediateResults = true;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 0.5;
		Logger.showDebugInformation = true;

		final MicrogliaSegmentationAndTracking microgliaSegmentationAndTracking =
				new MicrogliaSegmentationAndTracking(
						Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 ),
						settings );

		microgliaSegmentationAndTracking.run();

		final ArrayList< RandomAccessibleInterval< T > > labelings = microgliaSegmentationAndTracking.getLabelings();

		final ImagePlus labelMask = Utils.labelingsAsImagePlus( labelings );

		labelMask.show();
//
//		new FileSaver( labelMask ).saveAsTiff( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/microglia/MAX_5C-crop-t1-3-labelMasks.tif" );

	}
}
