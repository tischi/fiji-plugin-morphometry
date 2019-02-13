package de.embl.cba.morphometry;

import net.imglib2.Point;
import net.imglib2.util.LinAlgHelpers;

import static de.embl.cba.morphometry.Constants.X;
import static de.embl.cba.morphometry.Constants.Y;
import static java.lang.Math.acos;
import static java.lang.Math.atan;
import static java.lang.Math.toDegrees;

public abstract class Angles
{
	public static double angle2DToCoordinateSystemsAxisInDegrees( Point point )
	{
		final double[] vector = Vectors.asDoubles( point );

		return angle2DToCoordinateSystemsAxisInDegrees( vector );
	}

	public static double angle2DToCoordinateSystemsAxisInDegrees( double[] vector )
	{

		double angleToZAxisInDegrees;

		if ( vector[ Y ] == 0 )
		{
			angleToZAxisInDegrees = Math.signum( vector[ X ] ) * 90;
		}
		else
		{
			angleToZAxisInDegrees = toDegrees( atan( vector[ X ] / vector[ Y ] ) );

			if ( vector[ Y ] < 0 )
			{
				angleToZAxisInDegrees += 180;
			}
		}

		return angleToZAxisInDegrees;
	}


	public static double angleOfSpindleAxisToXaxisInRadians( final double[] vector )
	{
		double[] xAxis = new double[]{ 1, 0, 0};

		double angleInRadians = angleInRadians( vector, xAxis );

		return angleInRadians;
	}

	public static double angleInRadians( double[] v1, double[] v2 )
	{
		final double dot = LinAlgHelpers.dot( v2, v1 );

		double signum = Math.signum( v1[ 1 ] );

		signum = signum == 0 ? 1 : signum; // TODO: ????

		final double angle = -1.0 * signum * acos( dot / ( LinAlgHelpers.length( v2 ) * LinAlgHelpers.length( v1 ) ) );

		return angle;
	}

}
