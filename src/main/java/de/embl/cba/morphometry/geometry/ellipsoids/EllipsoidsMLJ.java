package de.embl.cba.morphometry.geometry.ellipsoids;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.type.logic.BitType;
import net.imglib2.view.Views;

import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.PHI;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.PSI;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.THETA;
import static java.lang.Math.*;

public abstract class EllipsoidsMLJ
{

	// Adapted from MorpholibJ


	public static EllipsoidMLJ computeParametersFromBinaryImage( RandomAccessibleInterval< BitType > binaryImg )
	{
		double[] sums = new double[ 3 ];
		double[] sumSquares = new double[ 6 ];

		final long numPixels = computeSumsAndSumSquares( binaryImg, sums, sumSquares );

		final double[] center = computeCenter( sums, numPixels );

		final double[] moments = computeMoments( sumSquares, center, numPixels );

		final Matrix momentsMatrix = getMomentsMatrix( moments );

		EllipsoidMLJ ellipsoidParameters = new EllipsoidMLJ();

		ellipsoidParameters.center = center;

		SingularValueDecomposition svd = new SingularValueDecomposition( momentsMatrix );

		ellipsoidParameters.radii = computeRadii( svd.getS() );

		ellipsoidParameters.eulerAnglesInDegrees = computeEllipsoidAlignmentAngles( svd.getU() );

		return ellipsoidParameters;
	}

	public static
	AffineTransform3D createAlignmentTransform( EllipsoidMLJ ellipsoidParameters )
	{
		AffineTransform3D translation = new AffineTransform3D();
		translation.translate( ellipsoidParameters.center  );
		translation = translation.inverse();

		AffineTransform3D rotation = new AffineTransform3D();
		rotation.rotate( Z, - toRadians( ellipsoidParameters.eulerAnglesInDegrees[ PHI ] ) );
		rotation.rotate( Y, - toRadians( ellipsoidParameters.eulerAnglesInDegrees[ THETA ] ) );
		rotation.rotate( X, - toRadians( ellipsoidParameters.eulerAnglesInDegrees[ PSI ] ) );

		AffineTransform3D combinedTransform = translation.preConcatenate( rotation );

		return combinedTransform;
	}

	private static double[] computeMoments( double[] sumSquares, double[] center, long numPixels )
	{
		double[] moments = new double[ 6 ];

		for ( int d : XYZ )
		{
			moments[ d ] = sumSquares[ d ] / numPixels - center[ d ] * center[ d ];
		}

		moments[ XY ] = sumSquares[ XY ] / numPixels - center[ X ] * center[ Y ];
		moments[ XZ ] = sumSquares[ XZ ] / numPixels - center[ X ] * center[ Z ];
		moments[ YZ ] = sumSquares[ YZ ] / numPixels - center[ Y ] * center[ Z ];

		return moments;
	}

	private static long computeSumsAndSumSquares( RandomAccessibleInterval< BitType > binaryImg,
												  double[] sums, double[] sumSquares )
	{

		Cursor< BitType > cursor = Views.iterable( binaryImg ).localizingCursor();

		long[] xyzPosition = new long[ 3 ];

		BitType one = new BitType( true );

		long numPixels = 0;

		while ( cursor.hasNext() )
		{
			if ( cursor.next().valueEquals( one ) )
			{
				numPixels++;

				cursor.localize( xyzPosition );

				for ( int d : XYZ )
				{
					sums[ d ] += xyzPosition[ d ];
				}

				for ( int d : XYZ )
				{
					sumSquares[ d ] += ( xyzPosition[ d ] * xyzPosition[ d ] );
				}

				sumSquares[ XY ] += xyzPosition[ X ] * xyzPosition[ Y ];
				sumSquares[ YZ ] += xyzPosition[ Y ] * xyzPosition[ Z ];
				sumSquares[ XZ ] += xyzPosition[ X ] * xyzPosition[ Z ];

			}
		}

		return numPixels;
	}


	private static double[] computeCenter( double[] sum, long n )
	{
		double[] center = new double[ 3 ];

		for ( int i = 0; i < 3; ++i )
		{
			center[ i ] = sum[ i ] / n;
		}

		return center;
	}


	private static Matrix getMomentsMatrix( double[] moments )
	{
		// run run parameters for each region
		Matrix matrix = new Matrix( 3, 3 );

		for ( int d : XYZ )
		{
			matrix.set( d, d, moments[ d ] );
		}

		matrix.set( 0, 1, moments[ XY ] );
		matrix.set( 0, 2, moments[ XZ ] );

		matrix.set( 1, 0, moments[ XY ] );
		matrix.set( 1, 2, moments[ YZ ] );

		matrix.set( 2, 0, moments[ XZ ] );
		matrix.set( 2, 1, moments[ YZ ] );

		return matrix;
	}


	private static double[] computeRadii( Matrix sdvS )
	{

		double[] radii = new double[ 3 ];

		for ( int d : XYZ )
		{
			radii[ d ] = sqrt( 5 ) * sqrt( sdvS.get( d, d ) );
		}

		return radii;
	}

	private static double[] computeEllipsoidAlignmentAngles( Matrix mat )
	{

		double[] angles = new double[ 3 ];

		// extract |cos(theta)|
		double tmp = hypot( mat.get(0, 0), mat.get(1, 0) );
		double phi, theta, psi;

		// avoid dividing by 0
		if (tmp > 16 * Double.MIN_VALUE)
		{
			// normal case: theta <> 0
			psi     = atan2( mat.get(2, 1), mat.get(2, 2));
			theta   = atan2(-mat.get(2, 0), tmp);
			phi     = atan2( mat.get(1, 0), mat.get(0, 0));
		}
		else
		{
			// theta is around 0
			psi     = atan2(-mat.get(1, 2), mat.get(1,1));
			theta   = atan2(-mat.get(2, 0), tmp);
			phi     = 0;
		}

		angles[ EllipsoidMLJ.PHI ] = toDegrees( phi );
		angles[ EllipsoidMLJ.THETA ] = toDegrees( theta );
		angles[ EllipsoidMLJ.PSI ] = toDegrees( psi );

		return angles;
	}



}
