package de.embl.cba.morphometry.geometry;

import de.embl.cba.morphometry.geometry.CoordinatesAndValues;

import java.util.ArrayList;

import static java.lang.Math.abs;

public abstract class CurveAnalysis
{
	public static CoordinatesAndValues computeDerivatives( CoordinatesAndValues coordinatesAndValues, int di )
	{
		final CoordinatesAndValues derivative = new CoordinatesAndValues();

		for ( int i = di / 2 + 1; i < coordinatesAndValues.values.size() - di / 2 - 1; ++i )
		{
			derivative.values.add( coordinatesAndValues.values.get( i + di / 2 ) - coordinatesAndValues.values.get( i - di / 2 ) );
			derivative.coordinates.add( 0.5 * ( coordinatesAndValues.coordinates.get( i + di / 2 ) + coordinatesAndValues.coordinates.get( i - di / 2 ) ));
		}

		return derivative;
	}

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

		final IndexAndValue indexAndValue = computeMaximumIndexAndValue( coordinatesAndValues );

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


	public static IndexAndValue computeMaximumIndexAndValue( CoordinatesAndValues coordinatesAndValues )
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

}
