package de.embl.cba.morphometry.geometry.ellipsoids;

import de.embl.cba.morphometry.Angles;
import de.embl.cba.morphometry.Vectors;
import ij.ImagePlus;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Vector3D;
import mcib3d.image3d.ImageByte;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.util.LinAlgHelpers;

import static de.embl.cba.morphometry.Constants.X;
import static de.embl.cba.morphometry.Constants.Y;
import static de.embl.cba.morphometry.Constants.Z;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.PHI;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.PSI;
import static de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ.THETA;
import static java.lang.Math.toRadians;

public abstract class Ellipsoids3DImageSuite
{
	public static EllipsoidVectors fitEllipsoid( ImagePlus mask )
	{
		final ImageByte imageByte = new ImageByte( mask );
		Objects3DPopulation objects3DPopulation = new Objects3DPopulation( imageByte, 0 );
		final Object3D object = objects3DPopulation.getObject( 0 );

//		System.out.println( "\n3D suite vectors");
//		System.out.println( object.getVectorAxis( 0 ) );
//		System.out.println( object.getVectorAxis( 1 ) );
//		System.out.println( object.getVectorAxis( 2 ) );

		final EllipsoidVectors ellipsoidVectors = new EllipsoidVectors();
		ellipsoidVectors.shortestAxis = object.getVectorAxis( 0 );
		ellipsoidVectors.middleAxis = object.getVectorAxis( 1 );
		ellipsoidVectors.longestAxis = object.getVectorAxis( 2 );
		ellipsoidVectors.center = object.getCenterAsArray();

		return ellipsoidVectors;
	}


	public static AffineTransform3D createAlignmentTransform( EllipsoidVectors ellipsoidVectors )
	{
		AffineTransform3D translation = new AffineTransform3D();
		translation.translate( ellipsoidVectors.center  );
		translation = translation.inverse();

		final Vector3D zAxis = new Vector3D( new Point3D( 0, 0, 1 ) );
		final double angleToZAxis = Math.PI / 180.0 * ellipsoidVectors.shortestAxis.angleDegrees( zAxis );
		final Vector3D axis = ellipsoidVectors.shortestAxis.crossProduct( zAxis );

		final double[] q = new double[ 4 ];
		LinAlgHelpers.quaternionFromAngleAxis( axis.getArray(), angleToZAxis, q );

		final double[][] m = new double[ 3 ][ 4 ];
		LinAlgHelpers.quaternionToR( q, m );

		AffineTransform3D rotation = new AffineTransform3D();
		rotation.set( m );

		AffineTransform3D combinedTransform = translation.preConcatenate( rotation );

		return combinedTransform;
	}


}
