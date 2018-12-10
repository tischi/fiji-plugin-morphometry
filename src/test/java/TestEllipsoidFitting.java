import bdv.util.BdvFunctions;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.geometry.ellipsoids.Ellipsoids3DImageSuite;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidsMLJ;
import de.embl.cba.transforms.utils.Transforms;
import ij.IJ;
import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.converter.Converters;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.PHI;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.PSI;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.THETA;

public class TestEllipsoidFitting
{
	public static < T extends RealType< T > & NativeType< T > > void main( String[] args )
	{

		ArrayList< String > paths = new ArrayList<>(  );
		//paths.add( TestEllipsoidFitting.class.getResource( "3d-ellipsoid-2.zip" ).getPath()  );
		//paths.add( TestEllipsoidFitting.class.getResource( "3d-ellipsoid-3.zip" ).getPath()  );
		//paths.add( TestEllipsoidFitting.class.getResource( "3d-ellipsoid-4.zip" ).getPath()  );
		paths.add( TestEllipsoidFitting.class.getResource( "dapi_mask_2.zip" ).getPath()  );

		for ( String path : paths )
		{

			final ImagePlus imagePlus = IJ.openImage( path );
			final RandomAccessibleInterval< T > wrap = ImageJFunctions.wrapReal( imagePlus );

			//BdvFunctions.show( wrap, "input" ).getBdvHandle().getViewerPanel();

			// MorpholibJ
			final RandomAccessibleInterval< BitType > mask = Converters.convert( wrap, ( i, o ) -> o.set( i.getRealDouble() > 1 ? true : false ), new BitType() );

			System.out.println( "\nDataset: " + path.toString() );

			final RandomAccessibleInterval aligned = createMLJAligned( mask );

			Ellipsoids3DImageSuite.fitEllipsoid( imagePlus );

			System.out.println( "\nDataset (round 2): " + path.toString() );

			final RandomAccessibleInterval aligned2 = createMLJAligned( aligned );

			Ellipsoids3DImageSuite.fitEllipsoid( Utils.asImagePlus( aligned, "" ) );

			System.out.println( "\nDataset (round 3): " + path.toString() );

			final RandomAccessibleInterval aligned3 = createMLJAligned( aligned2 );

			Ellipsoids3DImageSuite.fitEllipsoid( Utils.asImagePlus( aligned2, "" ) );

		}

	}

	public static RandomAccessibleInterval createMLJAligned( RandomAccessibleInterval< BitType > mask )
	{
		final EllipsoidMLJ ellipsoidParameters = EllipsoidsMLJ.computeParametersFromBinaryImage( mask );
		System.out.println( ellipsoidParameters.toString() );
		final AffineTransform3D alignmentTransform = EllipsoidsMLJ.createAlignmentTransform( ellipsoidParameters );
		final RandomAccessibleInterval aligned = Transforms.createTransformedView( mask, alignmentTransform );
		//BdvFunctions.show( aligned, s ).getBdvHandle().getViewerPanel();
		return aligned;
	}
}
