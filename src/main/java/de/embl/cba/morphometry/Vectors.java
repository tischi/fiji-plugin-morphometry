package de.embl.cba.morphometry;

import net.imglib2.Point;

public abstract class Vectors
{
	public static double[] asDoubles( Point point )
	{
		final double[] vector = new double[ 3 ];

		for ( int d = 0; d < 3; d++ )
		{
			vector[ d ] = point.getDoublePosition( d );
		}
		return vector;
	}
}
