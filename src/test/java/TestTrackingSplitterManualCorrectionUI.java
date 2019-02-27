import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.tracking.TrackingSplitterManualCorrectionUI;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class TestTrackingSplitterManualCorrectionUI
{
	public static < T extends RealType< T > & NativeType< T > > void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		ij.ui().showUI();

		final ImagePlus imagePlus = IJ.openImage(
				TestTrackingSplitterManualCorrectionUI.class.getResource( "microglia/MAX_5C-crop-t1-3-labelMasks.tif" ).getFile().toString() );

		final ArrayList< RandomAccessibleInterval< T > > labels =
				Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final TrackingSplitterManualCorrectionUI trackingSplitterManualCorrectionUI =
				new TrackingSplitterManualCorrectionUI(
						labels,
						100000L,
						"",
						null );

		while ( ! trackingSplitterManualCorrectionUI.isThisFrameFinished() ){
			Utils.wait( 100 );
		}

		final ArrayList< RandomAccessibleInterval< T > > correctedLabels
				= trackingSplitterManualCorrectionUI.getLabelings();

	}

}
