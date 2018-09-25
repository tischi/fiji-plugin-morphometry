package de.embl.cba.morphometry;

import de.embl.cba.morphometry.geometry.CentroidsParameters;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.Dataset;
import net.imagej.axis.LinearAxis;
import net.imglib2.*;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.LongArray;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.outofbounds.OutOfBounds;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.Intervals;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import org.scijava.log.LogService;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.viewing.BdvImageViewer.show;
import static java.lang.Math.*;

public class Utils
{
	public static int imagePlusChannelDimension = 2;

	public static void log( String message )
	{
		IJ.log( message );
	}


	public static void log( String message, LogService logService )
	{
		logService.info( message );
	}

	public static < T extends RealType< T > & NativeType< T > >
	CoordinatesAndValues computeAverageIntensitiesAlongAxis(
			RandomAccessibleInterval< T > rai, double maxAxisDist, int axis, double calibration )
	{
		final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

		for ( long coordinate = rai.min( axis ); coordinate <= rai.max( axis ); ++coordinate )
		{
			final IntervalView< T > intensitySlice = Views.hyperSlice( rai, axis, coordinate );
			coordinatesAndValues.coordinates.add( (double) coordinate * calibration );
			coordinatesAndValues.values.add( computeAverage( intensitySlice, maxAxisDist ) );
		}

		return coordinatesAndValues;
	}

	public static < T extends RealType< T > & NativeType< T > >
	CoordinatesAndValues computeMaximumIntensitiesAlongAxis(
			RandomAccessibleInterval< T > rai, double maxAxisDist, int axis, double calibration )
	{
		final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

		for ( long coordinate = rai.min( axis ); coordinate <= rai.max( axis ); ++coordinate )
		{
			final IntervalView< T > intensitySlice = Views.hyperSlice( rai, axis, coordinate );
			coordinatesAndValues.coordinates.add( (double) coordinate * calibration );
			coordinatesAndValues.values.add( computeMaximum( intensitySlice, maxAxisDist ) );
		}

		return coordinatesAndValues;
	}


	public static double sum( List<Double> a ){
		if (a.size() > 0) {
			double sum = 0;
			for (Double d : a) {
				sum += d;
			}
			return sum;
		}
		return 0;
	}
	public static double mean( List<Double> a ){
		double sum = sum( a );
		double mean = 0;
		mean = sum / ( a.size() * 1.0 );
		return mean;
	}

	public static double median( List<Double> a ){

		if ( a.size() == 1 ) return a.get( 0 );
		if ( a.size() == 0 ) return 0;

		int middle = a.size()/2;

		if (a.size() % 2 == 1) {
			return a.get(middle);
		} else {
			return (a.get(middle-1) + a.get(middle)) / 2.0;
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	CoordinatesAndValues computeAverageIntensitiesAlongAxisWithinMask( RandomAccessibleInterval< T > rai, RandomAccessibleInterval< BitType > mask, int axis, double calibration )
	{
		final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

		for ( long coordinate = rai.min( axis ); coordinate <= rai.max( axis ); ++coordinate )
		{
			final IntervalView< T > intensitySlice = Views.hyperSlice( rai, axis, coordinate );
			final IntervalView< BitType > maskSlice = Views.hyperSlice( mask, axis, coordinate );

			coordinatesAndValues.coordinates.add( (double) coordinate * calibration );
			coordinatesAndValues.values.add( computeAverage( intensitySlice, maskSlice ) );
		}

		return coordinatesAndValues;
	}

	public static < T extends RealType< T > & NativeType< T > >
	CoordinatesAndValues computeAverageIntensitiesAlongAxis(
			RandomAccessibleInterval< T > rai, int axis, double calibration )
	{
		final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

		for ( long coordinate = rai.min( axis ); coordinate <= rai.max( axis ); ++coordinate )
		{
			final IntervalView< T > intensitySlice = Views.hyperSlice( rai, axis, coordinate );
			coordinatesAndValues.coordinates.add( (double) coordinate * calibration );
			coordinatesAndValues.values.add( computeAverage( intensitySlice ) );
		}

		return coordinatesAndValues;
	}


	public static <T extends RealType<T> & NativeType< T > >
	ArrayList< RandomAccessibleInterval< T > > maskAllChannels(
			ArrayList< RandomAccessibleInterval< T > > channels,
			RandomAccessibleInterval< BitType > mask )
	{
		ArrayList< RandomAccessibleInterval< T > > maskedChannels = new ArrayList<>(  );

		long numChannels = channels.size();
		for ( int c = 0; c < numChannels; ++c )
		{
			final RandomAccessibleInterval< T > channel = Utils.copyAsArrayImg( channels.get( c ) );
			Utils.applyMask( channel, mask );
			maskedChannels.add( channel );
		}

		return maskedChannels;
	}

	public static CentroidsParameters computeCentroidsParametersAlongXAxis(
			RandomAccessibleInterval< BitType > rai,
			double calibration,
			double maxDistanceToCenter )
	{

		CentroidsParameters centroidsParameters = new CentroidsParameters();

		final double[] unitVectorInNegativeZDirection = new double[]{ 0, -1 };

		final double[] centralCentroid = computeCentroidPerpendicularToAxis( rai, X, 0 ); // this is not very robust, because the central one could be off

		for ( long coordinate = rai.min( X ); coordinate <= rai.max( X ); ++coordinate )
		{

			if ( Math.abs( coordinate * calibration ) < maxDistanceToCenter )
			{

				final double[] centroid = computeCentroidPerpendicularToAxis( rai, X, coordinate );
				final long numVoxels = computeNumberOfVoxelsPerpendicularToAxis( rai, X, coordinate );

				if ( centroid != null )
				{
					//				double[] centerDisplacementVector = subtract( centroid, centralCentroid ); // this is not very robust, because the central one could be off

					double[] centerDisplacementVector = centroid;
					double centerDisplacementLength = LinAlgHelpers.length( centerDisplacementVector );

					/**
					 *  centroid[ 0 ] is the y-axis coordinate
					 *  the sign of the y-axis coordinate determines the sign of the angle,
					 *  i.e. the direction of rotation
					 */
					final double angle = Math.signum( centerDisplacementVector[ 0 ] ) * 180 / Math.PI * acos( dotProduct( centerDisplacementVector, unitVectorInNegativeZDirection ) / centerDisplacementLength );

					centroidsParameters.distances.add( centerDisplacementLength * calibration );
					centroidsParameters.angles.add( angle );
					centroidsParameters.axisCoordinates.add( ( double ) coordinate * calibration );
					centroidsParameters.centroids.add( new RealPoint( coordinate * calibration, centroid[ 0 ] * calibration, centroid[ 1 ] * calibration ) );
					centroidsParameters.numVoxels.add( ( double ) numVoxels );
				}
			}
		}

		return centroidsParameters;

	}

	public static double vectorLength( double[] vector )
	{
		double norm = 0;

		for ( int d = 0; d < vector.length; ++d )
		{
			norm += vector[ d ] * vector[ d ];
		}

		norm = Math.sqrt( norm );

		return norm;
	}

	public static double dotProduct( double[] vector01, double[] vector02  )
	{
		double dotProduct = 0;

		for ( int d = 0; d < vector01.length; ++d )
		{
			dotProduct += vector01[ d ] * vector02[ d ];
		}

		return dotProduct;
	}

	public static double[] subtract( double[] vector01, double[] vector02  )
	{
		double[] subtraction = new double[ vector01.length ];

		for ( int d = 0; d < vector01.length; ++d )
		{
			subtraction[d] = vector01[ d ] - vector02[ d ];
		}

		return subtraction;
	}


	private static double[] computeCentroidPerpendicularToAxis( RandomAccessibleInterval< BitType > rai, int axis, long coordinate )
	{
		final IntervalView< BitType > slice = Views.hyperSlice( rai, axis, coordinate );

		int numHyperSliceDimensions = rai.numDimensions() - 1;

		final double[] centroid = new double[ numHyperSliceDimensions ];

		int numPoints = 0;

		final Cursor< BitType > cursor = slice.cursor();

		while( cursor.hasNext() )
		{
			if( cursor.next().get() )
			{
				for ( int d = 0; d < numHyperSliceDimensions; ++d )
				{
					centroid[ d ] += cursor.getLongPosition( d );
					numPoints++;
				}
			}
		}

		if ( numPoints > 0 )
		{
			for ( int d = 0; d < numHyperSliceDimensions; ++d )
			{
				centroid[ d ] /= 1.0D * numPoints;
			}
			return centroid;
		}
		else
		{
			return null;
		}
	}

	private static long computeNumberOfVoxelsPerpendicularToAxis( RandomAccessibleInterval< BitType > rai, int axis, long coordinate )
	{
		final IntervalView< BitType > slice = Views.hyperSlice( rai, axis, coordinate );

		int numHyperSliceDimensions = rai.numDimensions() - 1;

		final double[] centroid = new double[ numHyperSliceDimensions ];

		int numPoints = 0;

		final Cursor< BitType > cursor = slice.cursor();

		while ( cursor.hasNext() )
		{
			if ( cursor.next().get() )
			{
				numPoints++;
			}
		}

		return numPoints;
	}


	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createBlurredRai( RandomAccessibleInterval< T > rai, double sigma, double scaling )
	{
		ImgFactory< T > imgFactory = new ArrayImgFactory( rai.randomAccess().get()  );

		RandomAccessibleInterval< T > blurred = imgFactory.create( Intervals.dimensionsAsLongArray( rai ) );

		blurred = Views.translate( blurred, Intervals.minAsLongArray( rai ) );

		Gauss3.gauss( sigma / scaling, Views.extendBorder( rai ), blurred ) ;

		return blurred;
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createGaussFilteredArrayImg( RandomAccessibleInterval< T > rai, double[] sigmas )
	{
		ImgFactory< T > imgFactory = new ArrayImgFactory( rai.randomAccess().get()  );

		RandomAccessibleInterval< T > blurred = imgFactory.create( Intervals.dimensionsAsLongArray( rai ) );

		blurred = Views.translate( blurred, Intervals.minAsLongArray( rai ) );

		Gauss3.gauss( sigmas, Views.extendBorder( rai ), blurred ) ;

		return blurred;
	}


	public static  < T extends RealType< T > & NativeType< T > >
	void applyMask( RandomAccessibleInterval< T > rai, RandomAccessibleInterval< BitType > mask )
	{
		final Cursor< T > cursor = Views.iterable( rai ).cursor();
		final OutOfBounds< BitType > maskAccess = Views.extendZero( mask ).randomAccess();

  		while ( cursor.hasNext() )
		{
			cursor.fwd();
			maskAccess.setPosition( cursor );
			if ( ! maskAccess.get().get() )
			{
				cursor.get().setZero();
			}
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createAverageProjectionAlongAxis(
			RandomAccessibleInterval< T > rai, int d, double min, double max, double scaling )
	{
		Projection< T > projection = new Projection< T >(  rai, d,  new FinalInterval( new long[]{ (long) ( min / scaling) },  new long[]{ (long) ( max / scaling ) } ) );
		return projection.average();
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createSumProjectionAlongAxis(
			RandomAccessibleInterval< T > rai, int d, double min, double max, double scaling )
	{
		Projection< T > projection = new Projection< T >(  rai, d,  new FinalInterval( new long[]{ (long) ( min / scaling) },  new long[]{ (long) ( max / scaling ) } ) );
		return projection.sum();
	}

	public static < T extends RealType< T > & NativeType< T > > RandomAccessibleInterval< BitType > createBinaryImage( RandomAccessibleInterval< T > input, double doubleThreshold )
	{
		final ArrayImg< BitType, LongArray > binaryImage = ArrayImgs.bits( Intervals.dimensionsAsLongArray( input ) );

		T threshold = input.randomAccess().get().copy();
		threshold.setReal( doubleThreshold );

		final BitType one = new BitType( true );
		final BitType zero = new BitType( false );

		LoopBuilder.setImages( input, binaryImage ).forEachPixel( ( i, b ) ->
				{
					b.set( i.compareTo( threshold ) > 0 ?  one : zero );
				}
		);

		return binaryImage;

	}

	public static < T extends RealType< T > & NativeType< T > >
	AffineTransform3D createOrientationTransformation(
			RandomAccessibleInterval< T > rai, int longAxisDimension, double derivativeDelta, double calibration,
			boolean showPlots )
	{

		final CoordinatesAndValues coordinatesAndValues = computeAverageIntensitiesAlongAxis( rai, longAxisDimension, calibration );

		ArrayList< Double > absoluteDerivatives = Algorithms.computeAbsoluteDerivatives( coordinatesAndValues.values, (int) (derivativeDelta / calibration ));

		double maxLoc = computeMaxLoc( coordinatesAndValues.coordinates, absoluteDerivatives, null );

		System.out.println( "maxLoc = " + maxLoc );

		if ( showPlots )
		{
			Plots.plot( coordinatesAndValues.coordinates, coordinatesAndValues.values, "x", "intensity" );
			Plots.plot( coordinatesAndValues.coordinates, absoluteDerivatives, "x", "abs( derivative )" );
		}

		if ( maxLoc > 0 )
		{
			AffineTransform3D affineTransform3D = new AffineTransform3D();
			affineTransform3D.rotate( Z, toRadians( 180.0D ) );

			return affineTransform3D;
		}
		else
		{
			return new AffineTransform3D();
		}

	}

	public static double computeMaxLoc( CoordinatesAndValues coordinatesAndValues )
	{
		return computeMaxLoc( coordinatesAndValues.coordinates, coordinatesAndValues.values, null );
	}

	public static double computeMaxLoc( ArrayList< Double > coordinates, ArrayList< Double > values, double[] coordinateRangeMinMax )
	{
		double max = Double.MIN_VALUE;
		double maxLoc = coordinates.get( 0 );

		for ( int i = 0; i < values.size(); ++i )
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

	public static double[] getCalibration( Dataset dataset )
	{
		double[] calibration = new double[ 3 ];

		for ( int d : XYZ )
		{
			calibration[ d ] = ( ( LinearAxis ) dataset.getImgPlus().axis( d ) ).scale();
		}

		return calibration;
	}

	public static double[] getCalibration( ImagePlus imp )
	{
		double[] calibration = new double[ 3 ];

		calibration[ X ] = imp.getCalibration().pixelWidth;
		calibration[ Y ] = imp.getCalibration().pixelHeight;
		calibration[ Z ] = imp.getCalibration().pixelDepth;

		return calibration;
	}

	public static double[] get2dCalibration( ImagePlus imp )
	{
		double[] calibration = new double[ 2 ];

		calibration[ X ] = imp.getCalibration().pixelWidth;
		calibration[ Y ] = imp.getCalibration().pixelHeight;

		return calibration;
	}

	public static void correctCalibrationForSubSampling( double[] calibration, int subSampling )
	{
		for ( int d : XYZ )
		{
			calibration[ d ] *= subSampling;
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	double computeAverage( final RandomAccessibleInterval< T > rai )
	{
		final Cursor< T > cursor = Views.iterable( rai ).cursor();

		double average = 0;

		while ( cursor.hasNext() )
		{
			average += cursor.next().getRealDouble();
		}

		average /= Views.iterable( rai ).size();

		return average;
	}

	public static < T extends RealType< T > & NativeType< T > >
	double computeAverage( final RandomAccessibleInterval< T > rai, double maxAxisDist )
	{
		final Cursor< T > cursor = Views.iterable( rai ).cursor();

		double average = 0;
		long n = 0;
		double[] position = new double[ rai.numDimensions() ];

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( position );
			if ( Utils.vectorLength( position ) <= maxAxisDist)
			{
				average += cursor.get().getRealDouble();
				++n;
			}
		}

		average /= n;

		return average;
	}

	public static < T extends RealType< T > & NativeType< T > >
	double computeMaximum( final RandomAccessibleInterval< T > rai, double maxAxisDist )
	{
		final Cursor< T > cursor = Views.iterable( rai ).cursor();

		double max = - Double.MAX_VALUE;
		double[] position = new double[ rai.numDimensions() ];

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( position );
			if ( Utils.vectorLength( position ) <= maxAxisDist)
			{
				if( cursor.get().getRealDouble() > max )
				{
					max = cursor.get().getRealDouble();
				}
			}
		}

		return max;
	}


	public static < T extends RealType< T > & NativeType< T > >
	double computeAverage( final RandomAccessibleInterval< T > rai, final RandomAccessibleInterval< BitType > mask )
	{
		final Cursor< BitType > cursor = Views.iterable( mask ).cursor();
		final RandomAccess< T > randomAccess = rai.randomAccess();

		randomAccess.setPosition( cursor );

		double average = 0;
		long n = 0;

		while ( cursor.hasNext() )
		{
			if ( cursor.next().get() )
			{
				randomAccess.setPosition( cursor );
				average += randomAccess.get().getRealDouble();
				++n;
			}
		}

		average /= n;

		return average;
	}


	public static < T extends RealType< T > & NativeType< T > >
	List< RealPoint > computeMaximumLocation( RandomAccessibleInterval< T > blurred, int sigmaForBlurringAverageProjection )
	{
		Shape shape = new HyperSphereShape( sigmaForBlurringAverageProjection );

		List< RealPoint > points = Algorithms.findLocalMaximumValues( blurred, shape );

		return points;
	}

	public static List< RealPoint > asRealPointList( Point maximum )
	{
		List< RealPoint > realPoints = new ArrayList<>();
		final double[] doubles = new double[ maximum.numDimensions() ];
		maximum.localize( doubles );
		realPoints.add( new RealPoint( doubles) );

		return realPoints;
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > copyAsArrayImg( RandomAccessibleInterval< T > orig )
	{
		RandomAccessibleInterval< T > copy = new ArrayImgFactory( orig.randomAccess().get() ).create( orig );
		copy = Transforms.getWithAdjustedOrigin( orig, copy );
		LoopBuilder.setImages( copy, orig ).forEachPixel( ( c, o ) -> c.set( o ) );

		return copy;
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > copyAsEmptyArrayImg( RandomAccessibleInterval< T > rai )
	{
		RandomAccessibleInterval< T > copy = new ArrayImgFactory( rai.randomAccess().get() ).create( rai );
		copy = Transforms.getWithAdjustedOrigin( rai, copy );
		return copy;
	}

	public static < T extends RealType< T > & NativeType< T > >
	long[] getCenterLocation( RandomAccessibleInterval< T > rai )
	{
		int numDimensions = rai.numDimensions();

		long[] center = new long[ numDimensions ];

		for ( int d = 0; d < numDimensions; ++d )
		{
			center[ d ] = ( rai.max( d ) - rai.min( d ) ) / 2 + rai.min( d );
		}

		return center;

	}


	public static < T extends RealType< T > & NativeType< T > >
	boolean isBoundaryPixel( Cursor< T > cursor, RandomAccessibleInterval< T > rai )
	{
		int numDimensions = rai.numDimensions();
		final long[] position = new long[ numDimensions ];
		cursor.localize( position );


		for ( int d = 0; d < numDimensions; ++d )
		{
			if ( position[ d ] == rai.min( d ) ) return true;
			if ( position[ d ] == rai.max( d ) ) return true;
		}

		return false;
	}

	public static < T extends RealType< T > & NativeType< T > >
	boolean isLateralBoundaryPixel( Neighborhood< T > cursor, RandomAccessibleInterval< T > rai )
	{
		int numDimensions = rai.numDimensions();
		final long[] position = new long[ numDimensions ];
		cursor.localize( position );

		for ( int d = 0; d < numDimensions - 1; ++d )
		{
			if ( position[ d ] == rai.min( d ) ) return true;
			if ( position[ d ] == rai.max( d ) ) return true;
		}

		return false;

	}


	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > invertedView( RandomAccessibleInterval< T > input )
	{
		final double maximum = Algorithms.getMaximumValue( input );

		final RandomAccessibleInterval< T > inverted = Converters.convert( input, ( i, o ) -> {
			o.setReal( ( int ) ( maximum - i.getRealDouble() ) );
		},  Views.iterable( input ).firstElement() );

		return inverted;
	}


	public static long[] asLongs( double[] doubles )
	{
		final long[] longs = new long[ doubles.length ];

		for ( int i = 0; i < doubles.length; ++i )
		{
			longs[ i ] = (long) doubles[ i ];
		}

		return longs;
	}

	public static long[] divideAndReturnAsLongs( double[] doubles, double factor )
	{
		final long[] longs = new long[ doubles.length ];

		for ( int i = 0; i < doubles.length; ++i )
		{
			longs[ i ] = (long) ( doubles[ i ] / factor );
		}

		return longs;

	}

	public static void divide( double[] doubles, double factor )
	{
		for ( int i = 0; i < doubles.length; ++i )
		{
			doubles[ i ] /= factor;
		}
	}

	public static RandomAccessibleInterval< IntType >  asIntImg( ImgLabeling< Integer, IntType > labeling )
	{
		final RandomAccessibleInterval< IntType > intImg =
				Converters.convert( ( RandomAccessibleInterval< LabelingType< Integer > > ) labeling,
						( i, o ) -> {
							o.set( i.getIndex().getInteger() );
						}, new IntType() );

		return intImg;
	}

	public static double[] get2dDoubleArray( double value )
	{
		double[] array = new double[ 2 ];
		Arrays.fill( array, value );
		return array;
	}

	public static double[] get3dDoubleArray( double value )
	{
		double[] array = new double[ 3 ];
		Arrays.fill( array, value );
		return array;
	}

	public static FinalInterval getInterval( LabelRegion labelRegion )
	{
		final long[] min = Intervals.minAsLongArray( labelRegion );
		final long[] max = Intervals.maxAsLongArray( labelRegion );
		return new FinalInterval( min, max );
	}

	public static < T extends RealType< T > & NativeType< T > >
	void drawPoint( RandomAccessibleInterval< T > rai, double[] position, double radius, double calibration )
	{
		Shape shape = new HyperSphereShape( (int) ceil( radius / calibration ) );
		final RandomAccessible< Neighborhood< T > > nra = shape.neighborhoodsRandomAccessible( rai );
		final RandomAccess< Neighborhood< T > > neighborhoodRandomAccess = nra.randomAccess();

		neighborhoodRandomAccess.setPosition( asLongs( position )  );
		final Neighborhood< T > neighborhood = neighborhoodRandomAccess.get();

		final Cursor< T > cursor = neighborhood.cursor();
		while( cursor.hasNext() )
		{
			try
			{
				cursor.next().setReal( 200 );
			}
			catch ( ArrayIndexOutOfBoundsException e )
			{
				log( "[ERROR] Draw points out of bounds..." );
				break;
			}
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > getEnlargedRai( RandomAccessibleInterval< T > rai )
	{
		long[] min = new long[ 2 ];
		long[] max = new long[ 2 ];
		rai.max( max );
		for ( int d = 0; d < 2; ++d )
		{
			max[ d ] *= 1.2;
		}
		final FinalInterval interval = new FinalInterval( min, max );
		return Views.interval( Views.extendZero( rai ), interval );
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > getEnlargedRai2( RandomAccessibleInterval< T > rai, int border )
	{
		long[] min = new long[ 3 ];
		long[] max = new long[ 3 ];

		rai.min( min );
		rai.max( max );

		for ( int d = 0; d < 3; ++d )
		{
			min[ d ] -= border;
			max[ d ] += border;

		}

		final FinalInterval interval = new FinalInterval( min, max );
		return Views.interval( Views.extendZero( rai ), interval );
	}

	public static RandomAccessibleInterval< BitType > asMask( LabelRegion labelRegion )
	{
		RandomAccessibleInterval< BitType > rai = ArrayImgs.bits( Intervals.dimensionsAsLongArray( labelRegion ) );
		rai = Transforms.getWithAdjustedOrigin( labelRegion, rai  );
		final RandomAccess< BitType > randomAccess = rai.randomAccess();

		final LabelRegionCursor cursor = labelRegion.cursor();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			randomAccess.setPosition( cursor );
			randomAccess.get().set( true );
		}

		return rai;

	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > getMaskedAndCropped( RandomAccessibleInterval<T> image, LabelRegion labelRegion )
	{
		ImgFactory< T > imgFactory = new ArrayImgFactory( image.randomAccess().get()  );
		RandomAccessibleInterval< T > output = Views.translate( imgFactory.create( labelRegion ), Intervals.minAsLongArray( labelRegion )  ) ;

		final RandomAccess< T > imageRandomAccess = image.randomAccess();
		final RandomAccess< T > outputRandomAccess = output.randomAccess();
		final LabelRegionCursor cursor = labelRegion.cursor();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			imageRandomAccess.setPosition( cursor );
			outputRandomAccess.setPosition( cursor );
			outputRandomAccess.get().set( imageRandomAccess.get() );
		}

		return output;
	}

	public static < T extends RealType< T > & NativeType< T > >
	void setValues( RandomAccessibleInterval< T > rai, double value )
	{
		Cursor< T > cursor = Views.iterable( rai ).localizingCursor();

		double maxValue = Double.MIN_VALUE;

		while ( cursor.hasNext() )
		{
			cursor.next().setReal( value );
		}
	}

	public static < T extends RealType< T > & NativeType< T >  >
	ImgLabeling< Integer, IntType > asImgLabeling( RandomAccessibleInterval< T > rai )
	{
		RandomAccessibleInterval< IntType > labelImg = ArrayImgs.ints( Intervals.dimensionsAsLongArray( rai ) );
		labelImg = Transforms.getWithAdjustedOrigin( rai, labelImg );
		final ImgLabeling< Integer, IntType > labeling = new ImgLabeling<>( labelImg );

		final java.util.Iterator< Integer > labelCreator = new java.util.Iterator< Integer >()
		{
			int id = 0;

			@Override
			public boolean hasNext()
			{
				return true;
			}

			@Override
			public synchronized Integer next()
			{
				return id++;
			}
		};

		ConnectedComponents.labelAllConnectedComponents( ( RandomAccessible ) Views.extendBorder( rai ), labeling, labelCreator, ConnectedComponents.StructuringElement.EIGHT_CONNECTED );

		return labeling;
	}
}
