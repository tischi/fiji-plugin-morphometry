package de.embl.cba.morphometry.measurements;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class ObjectMeasurements
{

	public static final String COORDINATE = "Coordinate";
	public static final String SIZE_PIXEL_UNITS = "Size_Pixels";
	public static final String PIXEL_UNITS = "Pixels";
	public static final String SUM_INTENSITY = "SumIntensity";
	public static final String GOBAL_BACKGROUND_INTENSITY = "GobalBackgroundIntensity";
	public static final String SEP = "_";
	private static final String FRAME_UNITS = "Frames";

	public static void measurePositions( HashMap<Integer, Map<String, Object>> objectMeasurements, ImgLabeling<Integer, IntType> imgLabeling, double[] calibration )
	{
		String[] XYZ = new String[]{"X","Y","Z"};

		String unit = "";
		if ( calibration == null )
		{
			unit = PIXEL_UNITS;
		}

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final int label = ( int ) ( labelRegion.getLabel() );

			final double[] position = new double[ labelRegion.numDimensions() ];

			labelRegion.getCenterOfMass().localize( position );

			for ( int d = 0; d < position.length; ++d )
			{
				if ( calibration != null ) position[ d ] *= calibration[ d ];
				addMeasurement( objectMeasurements, label, COORDINATE + SEP + XYZ[ d ] + SEP + unit, position[ d ] );
			}

		}
	}

	public static void measureSizes( HashMap<Integer, Map<String, Object>> objectMeasurements,
									 ImgLabeling<Integer, IntType> imgLabeling )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final int label = ( int ) ( labelRegion.getLabel() );
			addMeasurement( objectMeasurements, label, PIXEL_UNITS, labelRegion.size() );
		}
	}

	public static void addMeasurement( HashMap< Integer, Map< String, Object > > objectMeasurements, int objectLabel, String name, Object value )
	{
		if ( ! objectMeasurements.keySet().contains( objectLabel ) )
		{
			objectMeasurements.put( objectLabel, new HashMap<>(  ) );
		}

		objectMeasurements.get( objectLabel ).put( name, value );
	}

	public static < T extends RealType< T > & NativeType< T > >
	void measureSumIntensities( HashMap< Integer, Map< String, Object > > objectMeasurements,
								ImgLabeling< Integer, IntType > imgLabeling,
								RandomAccessibleInterval< T > image,
								String channel )
	{
		final RandomAccess< T > imageRandomAccess = image.randomAccess();

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );

		for ( LabelRegion labelRegion : labelRegions )
		{
			long sum = measureSumIntensity( imageRandomAccess, labelRegion );
			addMeasurement( objectMeasurements, (int) labelRegion.getLabel(), SUM_INTENSITY + SEP + channel, sum );
		}
	}

	private static < T extends RealType< T > & NativeType< T > >
	long measureSumIntensity( RandomAccess< T > imageRandomAccess, LabelRegion labelRegion )
	{
		final LabelRegionCursor cursor = labelRegion.cursor();

		long sum = 0;

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			imageRandomAccess.setPosition( cursor );
			sum += imageRandomAccess.get().getRealDouble();
		}
		return sum;
	}

	public static < T extends RealType< T > & NativeType< T > >
	long measureSize( RandomAccessibleInterval< IntType > labeling,
					  int label )
	{

		final Cursor< IntType > labelCursor = Views.iterable( labeling ).localizingCursor();
		long size = 0;

		while ( labelCursor.hasNext() )
		{
			long value = labelCursor.next().getInteger();

			if( value == label )
			{
				size++;
			}
		}

		return size;

	}

	public static < T extends RealType< T > & NativeType< T > >
	double measureBgCorrectedSumIntensity( RandomAccessibleInterval< IntType > labeling,
										   int label,
										   RandomAccessibleInterval< T > image )
	{

		final Cursor< IntType > labelCursor = Views.iterable( labeling ).localizingCursor();
		final RandomAccess< T > intensityAccess = image.randomAccess();

		double sum = 0;
		double sumBg = 0;
		long nObject = 0;
		long nBackground = 0;
		int value;

		while ( labelCursor.hasNext() )
		{

			value = labelCursor.next().getInteger();

			if( value == label )
			{
				intensityAccess.setPosition( labelCursor );
				sum += intensityAccess.get().getRealDouble();
				nObject++;
			}
			else if ( value == 0 )
			{
				intensityAccess.setPosition( labelCursor );
				sumBg += intensityAccess.get().getRealDouble();
				nBackground++;
			}

		}

		final double meanBg = sumBg / nBackground;
		return ( sum - nObject * meanBg );

	}

	public static void addGlobalBackgroundMeasurement( HashMap<Integer, Map<String, Object>> objectMeasurements, ImgLabeling<Integer, IntType> imgLabeling, double offset )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final int label = ( int ) ( labelRegion.getLabel() );
			addMeasurement( objectMeasurements, label, GOBAL_BACKGROUND_INTENSITY, offset );
		}
	}

	public static ArrayList< String > printMeasurements( ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimePointList )
	{

		final Set< Integer > objectLabelsFirstTimePoint = measurementsTimePointList.get( 0 ).keySet();
		final Set< String > measurementNames = measurementsTimePointList.get( 0 ).get( objectLabelsFirstTimePoint.iterator().next() ).keySet();

		final ArrayList< String > lines = new ArrayList<>();

		String header = "Object_Label";

		header += "\t" + COORDINATE + SEP + "Time"+ SEP +FRAME_UNITS;

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
}
