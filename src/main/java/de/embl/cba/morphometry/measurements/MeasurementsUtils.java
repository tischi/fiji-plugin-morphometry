package de.embl.cba.morphometry.measurements;

import de.embl.cba.morphometry.Utils;

import net.imagej.table.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class MeasurementsUtils
{

	public static ArrayList< String > printMeasurements( ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimePointList )
	{

		final Set< Integer > objectLabelsFirstTimePoint = measurementsTimePointList.get( 0 ).keySet();
		final Set< String > measurementNames = measurementsTimePointList.get( 0 ).get( objectLabelsFirstTimePoint.iterator().next() ).keySet();

		final ArrayList< String > lines = new ArrayList<>();

		String header = "Object_Label";

		header += "\t" + ObjectMeasurements.COORDINATE + ObjectMeasurements.SEP + "Time"+ ObjectMeasurements.SEP + ObjectMeasurements.FRAME_UNITS;

		for ( String measurementName : measurementNames )
		{
			header += "\t" + measurementName ;
		}

		lines.add( header );

		for ( int t = 0; t < measurementsTimePointList.size(); ++t )
		{
			final HashMap< Integer, Map< String, Object > > measurements = measurementsTimePointList.get( t );

			final Set< Integer > objectLabels = measurements.keySet();

			for ( int label : objectLabels )
			{
				final Map< String, Object > measurementsMap = measurements.get( label );

				String values = String.format( "%05d", label );

				values += "\t" + String.format( "%05d", t );

				for ( String measurementName : measurementNames )
				{
					values += "\t" + measurementsMap.get( measurementName );
				}

				lines.add( values );
			}
		}

		return lines;
	}

	public static void saveMeasurements( File file, ArrayList<String> lines )
	{
		try (PrintWriter out = new PrintWriter( file ) )
		{
			for ( String line : lines )
			{
				out.println( line );
			}

			Utils.log( "\nSaved table to: " + file );
		}
		catch ( FileNotFoundException e )
		{
			e.printStackTrace();
		}
	}

	public static GenericTable createTable( ArrayList< String > lines  )
	{

		final DefaultGenericTable table = new DefaultGenericTable();

		// we create columns
		final String[] headers = lines.get( 0 ).split( "\t" );
		final int numColumns = headers.length;

		for ( int columnIndex = 0; columnIndex < numColumns; ++columnIndex )
		{
			GenericColumn column = new GenericColumn(headers[ columnIndex ] );

			for ( int i = 1; i < lines.size(); ++i )
			{
				String measurement = lines.get( i ).split( "\t" )[ columnIndex ];
				column.add( measurement );
			}

			table.add( column );
		}

		return table;

	}
}
