package playground;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaSettings;
import de.embl.cba.morphometry.tracking.SemiAutomatedTrackingSplitter;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class TestTrackingSplitter
{
	static class MasksAndIntensities
	{
		ImagePlus masks;
		ImagePlus intensities;
	}


	public static <T extends RealType< T > & NativeType< T > > void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final MasksAndIntensities masksAndIntensities = openLargeData();

		final ArrayList< RandomAccessibleInterval< T > > masks = Utils.get2DImagePlusMovieAsFrameList( masksAndIntensities.masks, 1 );
		final ArrayList< RandomAccessibleInterval< T > > intensities = Utils.get2DImagePlusMovieAsFrameList( masksAndIntensities.intensities, 1);

		MicrogliaSettings settings = new MicrogliaSettings();
		settings = MicrogliaSettings.configureSettings( settings );
		settings.calibration2D = Utils.getCalibration( masksAndIntensities.intensities );
		settings.showIntermediateResults = false;
		settings.opService = ij.op();

		final SemiAutomatedTrackingSplitter splitter = new SemiAutomatedTrackingSplitter( Utils.asMasks( masks ), intensities, settings );
		splitter.run();
		final ArrayList labelings = splitter.getLabelings();
		Utils.asImagePlusMovie( labelings, "Labels" ).show();
	}


	public static MasksAndIntensities openSmallData()
	{
		final MasksAndIntensities masksAndIntensities = new MasksAndIntensities();
		masksAndIntensities.masks = IJ.openImage( TestTrackingSplitter.class.getResource( "microglia/MAX_5C-crop-t1-3-masks.tif" ).getFile().toString() );
		masksAndIntensities.intensities = IJ.openImage( TestTrackingSplitter.class.getResource( "microglia/MAX_5C-crop-t1-3-intensities.tif" ).getFile().toString() );
		return masksAndIntensities;
	}

	public static MasksAndIntensities openLargeData()
	{
		final MasksAndIntensities masksAndIntensities = new MasksAndIntensities();
		masksAndIntensities.masks = IJ.openImage( "/Users/tischer/Documents/valerie-blanche-petegnief-CORBEL-microglia-quantification--data/5C/MAX_5C-cellMasks.tif" );
		masksAndIntensities.intensities = IJ.openImage( "/Users/tischer/Documents/valerie-blanche-petegnief-CORBEL-microglia-quantification--data/5C/MAX_5C-intensities.tif" );
		return masksAndIntensities;
	}
}
