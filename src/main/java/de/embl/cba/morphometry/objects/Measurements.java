package de.embl.cba.morphometry.objects;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.HashMap;
import java.util.Map;

public class Measurements
{

	public static final String CALIBRATED_POSITION = "CalibratedPosition";
	public static final String SIZE_PIXEL_UNITS = "PixelSize";
	public static final String SUM_INTENSITY = "SumIntensity";
	public static final String GOBAL_BACKGROUND_INTENSITY = "GobalBackgroundIntensity";

	public static void measureObjectPosition( HashMap<Integer, Map<String, Object>> objectMeasurements, ImgLabeling<Integer, IntType> imgLabeling, double[] calibration )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final int label = ( int ) ( labelRegion.getLabel() );
			final double[] position = new double[ calibration.length ];
			labelRegion.getCenterOfMass().localize( position );
			for ( int d = 0; d < position.length; ++d )
			{
				position[ d ] *= calibration[ d ];
			}

			addMeasurement( objectMeasurements, label, CALIBRATED_POSITION, position );
		}
	}

	public static void measureObjectPixelSizes( HashMap<Integer, Map<String, Object>> objectMeasurements, ImgLabeling<Integer, IntType> imgLabeling )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final int label = ( int ) ( labelRegion.getLabel() );
			addMeasurement( objectMeasurements, label, SIZE_PIXEL_UNITS, labelRegion.size() );
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
	void measureObjectSumIntensities( HashMap< Integer, Map< String, Object > > objectMeasurements, ImgLabeling< Integer, IntType > imgLabeling, RandomAccessibleInterval< T > image, String channel )
	{

		final RandomAccess< T > imageRandomAccess = image.randomAccess();

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final LabelRegionCursor cursor = labelRegion.cursor();

			final int label = ( int ) ( labelRegion.getLabel() );

			long sum = 0;

			while ( cursor.hasNext() )
			{
				cursor.fwd();
				imageRandomAccess.setPosition( cursor );
				sum += imageRandomAccess.get().getRealDouble();
			}

			addMeasurement( objectMeasurements, label, SUM_INTENSITY + "_" + channel, sum );

		}

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
}
