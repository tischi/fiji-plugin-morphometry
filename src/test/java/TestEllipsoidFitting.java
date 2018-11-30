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

public class TestEllipsoidFitting
{
	public static < T extends RealType< T > & NativeType< T > > void main( String[] args )
	{

		ArrayList< String > paths = new ArrayList<>(  );
		paths.add( TestEllipsoidFitting.class.getResource( "3d-ellipsoid-2.zip" ).getPath()  );
		paths.add( TestEllipsoidFitting.class.getResource( "3d-ellipsoid-3.zip" ).getPath()  );

		for ( String path : paths )
		{

			final ImagePlus imagePlus = IJ.openImage( path );
			final RandomAccessibleInterval< T > wrap = ImageJFunctions.wrapReal( imagePlus );
			//BdvFunctions.show( wrap, "input" ).getBdvHandle().getViewerPanel();

			// MorpholibJ
			final RandomAccessibleInterval< BitType > mask = Converters.convert( wrap, ( i, o ) -> o.set( i.getRealDouble() > 1 ? true : false ), new BitType() );
			final EllipsoidMLJ ellipsoidParameters = EllipsoidsMLJ.computeParametersFromBinaryImage( mask );
			printAngles( ellipsoidParameters );

			final AffineTransform3D alignmentTransform = EllipsoidsMLJ.createAlignmentTransform( ellipsoidParameters );
			final RandomAccessibleInterval aligned = Transforms.createTransformedView( mask, alignmentTransform );
			//BdvFunctions.show( aligned, "MLJ aligned" ).getBdvHandle().getViewerPanel();

			// 3D ImageSuite
			Ellipsoids3DImageSuite.fitEllipsoid( imagePlus );
		}

	}

	public static void printAngles( EllipsoidMLJ ellipsoidParameters )
	{
		System.out.println( "\nMLJ angles:");
		System.out.println( ellipsoidParameters.eulerAnglesInDegrees[ 0 ] );
		System.out.println( ellipsoidParameters.eulerAnglesInDegrees[ 1 ] );
		System.out.println( ellipsoidParameters.eulerAnglesInDegrees[ 2 ] );
	}
}
