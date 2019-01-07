import de.embl.cba.morphometry.Rois;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.translocation.TranslocationCommand;
import de.embl.cba.morphometry.translocation.TranslocationComputer;
import de.embl.cba.morphometry.translocation.TranslocationResult;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class TranslocationExample01
{
	public static  < T extends RealType< T > & NativeType< T > >
	void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		final ImagePlus imagePlus = IJ.openImage(
				TranslocationExample01.class.getResource(
						"translocation/test01.zip" ).getFile() );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", TranslocationExample01.class.getResource( "translocation/test01-rois.zip" ).getFile() );
		final ArrayList< FinalInterval > intervals = Rois.asIntervals( rm.getRoisAsArray() );

		ArrayList< RandomAccessibleInterval< T > > intensities = Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );


		final TranslocationComputer computer = new TranslocationComputer(
				intensities,
				intervals,
				ij.op() );

		final ArrayList< TranslocationResult > results = computer.getResults();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks =
				TranslocationCommand.createLabelMasks( intensities, results );

		TranslocationCommand.showLabelMasksAndIntensities( intensities, labelMasks );

		TranslocationCommand.plotTranslocations( results );

	}
}
