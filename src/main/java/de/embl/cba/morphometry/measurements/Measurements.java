package de.embl.cba.morphometry.measurements;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.morphometry.skeleton.SkeletonAnalyzer;
import de.embl.cba.tables.Tables;
import ij.ImagePlus;
import ij.measure.Calibration;
import net.imagej.ops.OpService;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.geom.real.Polygon2D;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;

import javax.swing.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Measurements
{

	public static final String COORDINATE = "Coordinate";

	public static final String VOLUME = "Volume";
	public static final String AREA = "Area";
	public static final String LENGTH = "Length";

	public static final String PERIMETER = "Perimeter";
	public static final String SURFACE = "Surface";

	public static final String PIXEL_UNITS = "Pixels";
	public static final String SUM_INTENSITY = "SumIntensity";
	public static final String NUM_BOUNDARY_PIXELS = "NumBoundaryPixels";

	public static final String GLOBAL_BACKGROUND_INTENSITY = "GlobalBackgroundIntensity";
	public static final String SKELETON_LENGTH = "SkeletonLength";
	public static final String SKELETON_NUMBER_OF_BRANCHPOINTS = "SkeletonNumBranchPoints";

	public static final String SEP = "_";
	public static final String FRAME_UNITS = "Frames";
	public static final String TIME = "Time";
	public static final String VOXEL_SPACING = "VoxelSpacing";
	public static final String FRAME_INTERVAL = "FrameInterval";

	public static String getVolumeName( int numDimensions )
	{
		if ( numDimensions == 1 ) return LENGTH;
		if ( numDimensions == 2 ) return AREA;
		if ( numDimensions == 3 ) return VOLUME;

		return null;
	}

	public static String getSurfaceName( int numDimensions )
	{
		if ( numDimensions == 1 ) return LENGTH;
		if ( numDimensions == 2 ) return PERIMETER;
		if ( numDimensions == 3 ) return SURFACE;

		return null;
	}


	public static void measurePositions(
			HashMap< Integer, Map< String, Object > > objectMeasurements,
			ImgLabeling<Integer, IntType> imgLabeling,
			double[] calibration )
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

			final double[] position = new double[ 3 ];

			labelRegion.getCenterOfMass().localize( position );

			for ( int d = 0; d < position.length; ++d )
			{
				if ( calibration != null ) position[ d ] *= calibration[ d ];

				addMeasurement(
						objectMeasurements,
						label,
						COORDINATE + SEP + XYZ[ d ] + SEP + unit,
						position[ d ] );
			}
		}
	}

	public static void measureVolumes( HashMap<Integer, Map<String, Object>> objectMeasurements,
									   ImgLabeling<Integer, IntType> imgLabeling )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final int label = ( int ) ( labelRegion.getLabel() );
			addMeasurement( objectMeasurements, label, getVolumeName( labelRegion.numDimensions() ) + SEP + PIXEL_UNITS, labelRegion.size() );
		}
	}

	public static void measureSurface( HashMap<Integer, Map<String, Object>> objectMeasurements,
									   ImgLabeling<Integer, IntType> imgLabeling,
									   OpService opService )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			final int label = ( int ) ( labelRegion.getLabel() );

			final RandomAccessibleInterval< BitType > mask = Regions.labelRegionAsMask( labelRegion );

			// See: https://forum.image.sc/t/measure-surface-perimeter-in-imglib2/21213

			final Polygon2D contour = opService.geom().contour( mask, true );
			final double boundarySize = opService.geom().boundarySize( contour ).getRealDouble();

			addMeasurement( objectMeasurements, label, getSurfaceName( labelRegion.numDimensions() ) + SEP + PIXEL_UNITS, boundarySize );
		}
	}


	public static void measureSkeletons( HashMap<Integer, Map<String, Object>> objectMeasurements,
										 ImgLabeling<Integer, IntType> imgLabeling,
										 RandomAccessibleInterval< BitType > skeleton,
										 OpService opService )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );

		for ( LabelRegion labelRegion : labelRegions )
		{
			final RandomAccessibleInterval< BitType > regionSkeleton =
					Regions.getMaskedAndCropped( skeleton, labelRegion );

			final SkeletonAnalyzer skeletonAnalyzer =
					new SkeletonAnalyzer( regionSkeleton, opService );

			final int label = ( int ) ( labelRegion.getLabel() );

			addMeasurement( objectMeasurements,
					label,
					 SKELETON_LENGTH + SEP + PIXEL_UNITS,
					skeletonAnalyzer.getSkeletonLength() );

			addMeasurement( objectMeasurements,
					label,
					SKELETON_NUMBER_OF_BRANCHPOINTS + SEP + PIXEL_UNITS,
					skeletonAnalyzer.getNumBranchPoints() );

		}
	}


	public static void addMeasurement(
			HashMap< Integer, Map< String, Object > > objectMeasurements,
			int objectLabel, String name, Object value )
	{
		if ( ! objectMeasurements.keySet().contains( objectLabel ) )
		{
			objectMeasurements.put( objectLabel, new HashMap<>(  ) );
		}

		objectMeasurements.get( objectLabel ).put( name, value );
	}

	public static < T extends RealType< T > & NativeType< T > >
	double sumIntensity( RandomAccessibleInterval< T > image,
						 RandomAccessibleInterval< BitType > mask

	)
	{
		final Cursor< BitType > maskCursor = Views.iterable( mask ).cursor();
		final RandomAccess< T > intensityAccess = image.randomAccess();

		double sum = 0;

		while ( maskCursor.hasNext() )
		{
			if ( maskCursor.next().get() )
			{
				intensityAccess.setPosition( maskCursor );
				sum += intensityAccess.get().getRealDouble();
			}


		}

		return sum;
	}

	public static < T extends RealType< T > & NativeType< T > >
	double meanIntensity( RandomAccessibleInterval< T > image,
						 RandomAccessibleInterval< BitType > mask

	)
	{
		final Cursor< BitType > maskCursor = Views.iterable( mask ).cursor();
		final RandomAccess< T > intensityAccess = image.randomAccess();

		double sum = 0;
		long n = 0;

		while ( maskCursor.hasNext() )
		{
			if ( maskCursor.next().get() )
			{
				intensityAccess.setPosition( maskCursor );
				sum += intensityAccess.get().getRealDouble();
				n++;
			}


		}

		return sum / n;
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

	private static < T extends RealType< T > & NativeType< T > >
	long measureNumBoundaryPixels( LabelRegion labelRegion, long[] min, long[] max )
	{
		final LabelRegionCursor cursor = labelRegion.cursor();

		long numBoundaryPixels = 0;

		final int numDimensions = min.length;
		final long[] position = new long[ numDimensions ];
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( position );
			for ( int d = 0; d < numDimensions ; d++ )
			{
				if ( position[ d ] == min[ d ] || position[ d ] == max[ d ])
				{
					numBoundaryPixels++;
					break;
				}
			}
		}
		return numBoundaryPixels;
	}

	public static < I extends IntegerType< I >  >
	long measureSizeInPixels( RandomAccessibleInterval< I > labeling,
							  int label )
	{
		final Cursor< I > labelCursor = Views.iterable( labeling ).localizingCursor();
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

	public static
	long measureSizeInPixels( RandomAccessibleInterval< BitType > mask )
	{

		final Cursor< BitType > cursor = Views.iterable( mask ).cursor();
		long size = 0;

		while ( cursor.hasNext() )
		{
			if( cursor.next().get() )
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
			addMeasurement( objectMeasurements, label, GLOBAL_BACKGROUND_INTENSITY, offset );
		}
	}

	public static void saveRowsToFile( File file, ArrayList<String> lines )
	{
		try ( PrintWriter out = new PrintWriter( file ) )
		{
			for ( String line : lines )
			{
				out.println( line );
			}

			Logger.log( "\nSaved table to: " + file );
		}
		catch ( FileNotFoundException e )
		{
			e.printStackTrace();
		}
	}

	public static JTable asTable( HashMap< Integer, Map< String, Object > > objectMeasurements )
	{
		final ArrayList< HashMap< Integer, Map< String, Object > > > timepoints = new ArrayList<>();
		timepoints.add( objectMeasurements );
		return Tables.createJTableFromStringList( measurementsAsTableRowsStringList( timepoints, "\t" ), "\t" );
	}

	public static JTable asTable( ArrayList< HashMap< Integer, Map< String, Object > > > timepoints )
	{
		return Tables.createJTableFromStringList( measurementsAsTableRowsStringList( timepoints, "\t" ), "\t" );
	}

	public static ArrayList< String > measurementsAsTableRowsStringList(
			ArrayList< HashMap< Integer,
			Map< String, Object > > > measurementsTimePointList,
			String delim )
	{

		final Set< Integer > objectLabelsFirstTimePoint = measurementsTimePointList.get( 0 ).keySet();
		final Set< String > measurementNames = measurementsTimePointList.get( 0 ).get( objectLabelsFirstTimePoint.iterator().next() ).keySet();

		final ArrayList< String > lines = new ArrayList<>();

		String header = "Object_Label";

		header += delim + COORDINATE + SEP + TIME + SEP + FRAME_UNITS;

		for ( String measurementName : measurementNames )
		{
			header += delim + measurementName ;
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

				values += delim + String.format( "%05d", t + 1 ); // convert to one-based time points

				for ( String measurementName : measurementNames )
				{
					values += delim + measurementsMap.get( measurementName );
				}

				lines.add( values );
			}
		}

		return lines;
	}

	public static ArrayList< HashMap< Integer, Map< String, Object > > >  initMeasurements( int numTimepoints )
	{
		ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimepointList = new ArrayList<>();

		for ( int t = 0; t < numTimepoints; ++t )
		{
			final HashMap< Integer, Map< String, Object > > objectMeasurements = new HashMap<>();
			measurementsTimepointList.add( objectMeasurements );
		}

		return measurementsTimepointList;
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


	public static < T extends RealType< T > & NativeType< T > >
	void measureNumBoundaryPixels( HashMap< Integer, Map< String, Object > > objectMeasurements,
								   ImgLabeling< Integer, IntType > imgLabeling )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );

		final int numDimensions = imgLabeling.numDimensions();
		final long[] imgBoundaryMin = new long[ numDimensions ];
		final long[] imgBoundaryMax = new long[ numDimensions ];
		imgLabeling.min( imgBoundaryMin );
		imgLabeling.max( imgBoundaryMax );

		for ( LabelRegion labelRegion : labelRegions )
		{
			long numBoundaryPixels = measureNumBoundaryPixels( labelRegion, imgBoundaryMin, imgBoundaryMax );
			addMeasurement( objectMeasurements, (int) labelRegion.getLabel(), NUM_BOUNDARY_PIXELS, numBoundaryPixels );
		}
	}

	public static void addCalibration(
			ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimepointList,
			ImagePlus imagePlus )
	{
		final Calibration calibration = imagePlus.getCalibration();
		final double pixelWidth = calibration.pixelWidth;
		final double pixelHeight = calibration.pixelHeight;
		final double pixelDepth = calibration.pixelDepth;
		final String unit = calibration.getUnit();
		final double frameInterval = calibration.frameInterval;
		final String timeUnit = calibration.getTimeUnit();

		for ( int t = 0; t < measurementsTimepointList.size(); ++t )
		{
			final HashMap< Integer, Map< String, Object > > measurements =
					measurementsTimepointList.get( t );

			final Set< Integer > labels = measurements.keySet();

			for( int label : labels )
			{
				measurements.get( label ).put(
						VOXEL_SPACING + SEP + "X",
						pixelWidth );

				measurements.get( label ).put(
						VOXEL_SPACING + SEP + "Y",
						pixelHeight );

				measurements.get( label ).put(
						VOXEL_SPACING + SEP + "Z",
						pixelDepth );

				measurements.get( label ).put(
						VOXEL_SPACING + SEP + "Unit",
						unit );

				measurements.get( label ).put(
						FRAME_INTERVAL,
						frameInterval );

				measurements.get( label ).put(
						FRAME_INTERVAL + SEP + "Unit",
						timeUnit );
			}

		}
	}
}
