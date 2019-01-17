import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public class TranslocationCellEdgeDetection
{
	public static < R extends RealType< R > & NativeType< R > >
	void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		ij.ui().showUI();

		final ImagePlus imagePlus = IJ.openImage(
				TranslocationTestCommand.class.getResource(
						"translocation/test01-singleSlice.zip" ).getFile() );

		final RandomAccessibleInterval< R > intensity = ImageJFunctions.wrapReal( imagePlus );

		final OpService opService = ij.op();
		final RandomAccessibleInterval< R > gauss =
				opService.filter().gauss(
						intensity,
						1 );

		final RandomAccessibleInterval< R > erode = Algorithms.erode( gauss, 3 );

		RandomAccessibleInterval< R > gradient = Utils.createEmptyCopy( erode );

		LoopBuilder.setImages( gradient, gauss, erode ).forEachPixel( ( g, i, e ) ->
				{
					g.setReal( i.getRealDouble() - e.getRealDouble() );
				}
		);

//		final RandomAccessibleInterval< R > gradient = Algorithms.computeGradient( gauss,
//				new HyperSphereShape( 1 ) );


//		final RandomAccessibleInterval< R > gradient = Algorithms.computeGradient( gauss,
//				new HyperSphereShape( 1 ) );

		RandomAccessibleInterval< BitType > mask = ArrayImgs.bits(
				Intervals.dimensionsAsLongArray( intensity ) );

		opService.threshold().isoData(
				Views.iterable( mask ),
				Views.iterable( gradient ) );

		final RandomAccessibleInterval< BitType > thin = ArrayImgs.bits(
				Intervals.dimensionsAsLongArray( mask ) );

		opService.morphology().thinGuoHall( thin, mask );



		ImageJFunctions.show( intensity, "intensity" );
		ImageJFunctions.show( gauss, "gauss" );
		ImageJFunctions.show( gradient, "gradient" );
		ImageJFunctions.show( mask, "mask" );
		ImageJFunctions.show( thin, "thin" );



	}

}
