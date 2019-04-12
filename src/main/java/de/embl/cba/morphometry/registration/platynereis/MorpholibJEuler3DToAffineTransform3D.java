package de.embl.cba.morphometry.registration.platynereis;

import net.imglib2.realtransform.AffineTransform3D;

public class MorpholibJEuler3DToAffineTransform3D
{

	/**
	 *
	 * @return an affine transform performing the Euler transform, in millimeter units
	 */
	public static AffineTransform3D convert( double[] rotationCenterInMicrometer, double[] eulerAnglesInDegrees )
	{

		// TODO: move to image-transform-converters repo


		// convert
		//
		final double[] rotationCentrePositive = new double[ 3 ];
		final double[] rotationCentreNegative = new double[ 3 ];

		for ( int d = 0; d < 3; ++d )
		{
			rotationCentrePositive[ d ] = rotationCenterInMicrometer[ d ];
			rotationCentreNegative[ d ] = - rotationCenterInMicrometer[ d ];
		}


		final AffineTransform3D transform3D = new AffineTransform3D();

		// rotate around rotation centre
		//

		// make rotation centre the image centre
		transform3D.translate( rotationCentreNegative );

		// rotate
		for ( int d = 0; d < 3; ++d )
		{
			transform3D.rotate( d, eulerAnglesInDegrees[ d ] );
		}


		// move image centre back
		final AffineTransform3D translateBackFromRotationCentre = new AffineTransform3D();
		translateBackFromRotationCentre.translate( rotationCentrePositive );
		transform3D.preConcatenate( translateBackFromRotationCentre );

		// translate
		//
		//transform3D.translate( translationInMillimeters );

		return transform3D;
	}
}
