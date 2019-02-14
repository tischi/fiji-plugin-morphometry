package de.embl.cba.morphometry.geometry;

import java.util.ArrayList;

import static java.lang.Math.abs;

public abstract class CurveAnalysis
{
	public static CoordinatesAndValues derivative( CoordinatesAndValues coordinatesAndValues, int di )
	{
		final CoordinatesAndValues derivative = new CoordinatesAndValues();

		for ( int i = di / 2 + 1; i < coordinatesAndValues.values.size() - di / 2 - 1; ++i )
		{
			derivative.values.add( coordinatesAndValues.values.get( i + di / 2 ) - coordinatesAndValues.values.get( i - di / 2 ) );
			derivative.coordinates.add( 0.5 * ( coordinatesAndValues.coordinates.get( i + di / 2 ) + coordinatesAndValues.coordinates.get( i - di / 2 ) ));
		}

		return derivative;
	}

	// TODO:
//	public static CoordinatesAndValues derivative( CoordinateToValue cv, int di )
//	{
//		final CoordinatesAndValues derivative = new CoordinatesAndValues();
//
//		for ( int i = di / 2 + 1; i < cv.size() - di / 2 - 1; ++i )
//		{
//			derivative.values.add( cv.values.get( i + di / 2 ) - cv.values.get( i - di / 2 ) );
//			derivative.coordinates.add( 0.5 * ( cv.coordinates.get( i + di / 2 ) + cv.coordinates.get( i - di / 2 ) ));
//		}
//
//		return derivative;
//	}

	public static ArrayList< Double > computeAbsoluteDerivatives( ArrayList< Double > values, int di )
	{
		final ArrayList< Double > derivatives = new ArrayList<>();

		for ( int i = di / 2 + 1; i < values.size() - di / 2 - 1; ++i )
		{
			derivatives.add( abs( values.get( i + di / 2 ) - values.get( i - di / 2 ) ) );
		}

		return derivatives;
	}

	public static double computeFWHM( CoordinatesAndValues coordinatesAndValues )
	{

		final IndexAndValue indexAndValue = maximumIndexAndValue( coordinatesAndValues );

		final int n = coordinatesAndValues.values.size();

		double halfMaxLoc1 = 0.0, halfMaxLoc2 = 0.0;

		for ( int i = indexAndValue.index; i < n; i++ )
		{
			if ( coordinatesAndValues.values.get( i ) < indexAndValue.value / 2.0 )
			{
				halfMaxLoc2 = 0.5 * ( coordinatesAndValues.values.get( i - 1 ) + coordinatesAndValues.values.get( i ) );
				break;
			}
		}

		for ( int i = indexAndValue.index; i >= 0; i-- )
		{
			if ( coordinatesAndValues.values.get( i ) < indexAndValue.value / 2.0 )
			{
				halfMaxLoc1 = 0.5 * ( coordinatesAndValues.values.get( i + 1 ) + coordinatesAndValues.values.get( i ) );
				break;
			}
		}

		return Math.abs( halfMaxLoc2 - halfMaxLoc1 );
	}


	public static IndexAndValue maximumIndexAndValue( CoordinatesAndValues coordinatesAndValues )
	{
		final int n = coordinatesAndValues.coordinates.size();

		final IndexAndValue indexAndValue = new IndexAndValue();
		indexAndValue.value = - Double.MAX_VALUE;

		for ( int i = 0; i < n; i++ )
		{
			if ( coordinatesAndValues.values.get( i ) > indexAndValue.value )
			{
				indexAndValue.value = coordinatesAndValues.values.get( i );
				indexAndValue.index = i;
			}
		}

		return indexAndValue;
	}

	public static double maxLocCoordinate( CoordinatesAndValues coordinatesAndValues )
	{
		final int n = coordinatesAndValues.coordinates.size();

		final IndexAndValue maxLocIndexAndValue = new IndexAndValue();
		maxLocIndexAndValue.value = - Double.MAX_VALUE;

		for ( int i = 0; i < n; i++ )
		{
			if ( coordinatesAndValues.values.get( i ) > maxLocIndexAndValue.value )
			{
				maxLocIndexAndValue.value = coordinatesAndValues.values.get( i );
				maxLocIndexAndValue.index = i;
			}
		}

		final Double maxLocCoordinate = coordinatesAndValues.values.get( maxLocIndexAndValue.index );

		return maxLocCoordinate;
	}


	public static double minLoc( CoordinatesAndValues coordinatesAndValues )
	{
		return minLoc( coordinatesAndValues, null );
	}

	public static double minLoc( CoordinatesAndValues coordinatesAndValues, double[] coordinateRangeMinMax )
	{
		final int n = coordinatesAndValues.coordinates.size();
		final ArrayList< Double > coordinates = coordinatesAndValues.coordinates;
		final ArrayList< Double > values = coordinatesAndValues.values;

		final IndexAndValue minLocIndexAndValue = new IndexAndValue();
		minLocIndexAndValue.value = Double.MAX_VALUE;

		for ( int i = 0; i < n; i++ )
		{

			if ( coordinateRangeMinMax != null )
			{
				if ( coordinates.get( i ) < coordinateRangeMinMax[ 0 ] ) continue;
				if ( coordinates.get( i ) > coordinateRangeMinMax[ 1 ] ) continue;
			}

			if ( values.get( i ) < minLocIndexAndValue.value )
			{
				minLocIndexAndValue.value = coordinatesAndValues.values.get( i );
				minLocIndexAndValue.index = i;
			}
		}

		final Double maxLocCoordinate = coordinatesAndValues.coordinates.get( minLocIndexAndValue.index );

		return maxLocCoordinate;
	}

	public static double maxLoc( CoordinatesAndValues coordinatesAndValues )
	{
		return maxLoc( coordinatesAndValues, null );
	}

	public static double maxLoc( CoordinatesAndValues coordinatesAndValues, double[] coordinateRangeMinMax )
	{
		final ArrayList< Double > coordinates = coordinatesAndValues.coordinates;
		final ArrayList< Double > values = coordinatesAndValues.values;
		final int n = values.size();

		double max = - Double.MAX_VALUE;
		double maxLoc = coordinates.get( 0 );

		for ( int i = 0; i < n; ++i )
		{
			if ( coordinateRangeMinMax != null )
			{
				if ( coordinates.get( i ) < coordinateRangeMinMax[ 0 ] ) continue;
				if ( coordinates.get( i ) > coordinateRangeMinMax[ 1 ] ) continue;
			}

			if ( values.get( i ) > max )
			{
				max = values.get( i );
				maxLoc = coordinates.get( i );
			}
		}

		return maxLoc;
	}

	public static double[] leftMaxAndRightMinLoc( CoordinatesAndValues coordinatesAndValues )
	{
		double[] rangeMinMax = new double[ 2 ];
		double[] locations = new double[ 2 ];

		// left
		rangeMinMax[ 0 ] = - Double.MAX_VALUE;
		rangeMinMax[ 1 ] = 0;
		locations[ 0 ] = maxLoc( coordinatesAndValues, rangeMinMax );

		// right
		rangeMinMax[ 0 ] = 0;
		rangeMinMax[ 1 ] = Double.MAX_VALUE;
		locations[ 1 ] = minLoc( coordinatesAndValues, rangeMinMax );

		return locations;
	}
}
