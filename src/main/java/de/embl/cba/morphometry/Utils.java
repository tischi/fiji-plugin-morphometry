package de.embl.cba.morphometry;

import de.embl.cba.morphometry.geometry.CentroidsParameters;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.process.LUT;
import net.imagej.Dataset;
import net.imagej.axis.LinearAxis;
import net.imglib2.*;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
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
import net.imglib2.img.basictypeaccess.array.IntArray;
import net.imglib2.img.basictypeaccess.array.LongArray;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.ops.parse.token.Int;
import net.imglib2.outofbounds.OutOfBounds;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.roi.labeling.*;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.AbstractIntegerType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.util.Intervals;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.util.Util;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import org.scijava.log.LogService;

import java.awt.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.List;

import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;
import static java.lang.Math.*;

public class Utils
{
	public static int imagePlusChannelDimension = 2;

	public static String logFilePath = null;


	public static void setNewLogFilePath( String aLogFilePath )
	{
		logFilePath = aLogFilePath;
		createLogFile();
	}

	public static void log( String message )
	{
		IJ.log( message );

		if ( logFilePath != null )
		{
			File logFile = new File( logFilePath );

			if ( ! logFile.exists() )
			{
				createLogFile();
			}
			else
			{
				writeToLogFile( message + "\n" );
			}
		}

	}

	public static void writeToLogFile( String message )
	{
		try {
			Files.write( Paths.get( logFilePath ), message.getBytes(), StandardOpenOption.APPEND);
		}
		catch (IOException e) {
			//exception handling left as an exercise for the reader
		}
	}

	public static void createLogFile()
	{
		PrintWriter writer = null;
		try
		{
			writer = new PrintWriter( logFilePath, "UTF-8" );
			writer.println( "Start logging..." );
		}
		catch ( FileNotFoundException e )
		{
			e.printStackTrace();
		}
		catch ( UnsupportedEncodingException e )
		{
			e.printStackTrace();
		}
		writer.close();
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
			RandomAccessibleInterval< BitType > mask,
			boolean showImages )
	{
		ArrayList< RandomAccessibleInterval< T > > maskedChannels = new ArrayList<>(  );

		long numChannels = channels.size();

		for ( int c = 0; c < numChannels; ++c )
		{
			final RandomAccessibleInterval< T > channel = Utils.copyAsArrayImg( channels.get( c ) );
			Utils.applyMask( channel, mask );
			maskedChannels.add( channel );

			if ( showImages ) show( channel, "masked channel " + c,  Transforms.origin(), 1.0 );

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


	public static double computeMinLoc( ArrayList< Double > coordinates, ArrayList< Double > values, double[] coordinateRangeMinMax )
	{
		double minValue = Double.MAX_VALUE;
		double minLoc = coordinates.get( 0 );

		for ( int i = 0; i < values.size(); ++i )
		{
			if ( coordinateRangeMinMax != null )
			{
				if ( coordinates.get( i ) < coordinateRangeMinMax[ 0 ] ) continue;
				if ( coordinates.get( i ) > coordinateRangeMinMax[ 1 ] ) continue;
			}

			if ( values.get( i ) < minValue )
			{
				minValue = values.get( i );
				minLoc = coordinates.get( i );
			}
		}

		return minLoc;
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


	public static double[] getCalibrationWithFixedZ( ImagePlus imp )
	{
		double[] calibration = new double[ 3 ];

		calibration[ X ] = imp.getCalibration().pixelWidth;
		calibration[ Y ] = imp.getCalibration().pixelHeight;
		calibration[ Z ] = imp.getCalibration().pixelDepth;

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
	RandomAccessibleInterval< T > createEmptyArrayImg( RandomAccessibleInterval< T > rai )
	{
		RandomAccessibleInterval< T > newImage = new ArrayImgFactory( rai.randomAccess().get() ).create( rai );
		newImage = Transforms.getWithAdjustedOrigin( rai, newImage );
		return newImage;
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

	public static double[] as2dDoubleArray( double value )
	{
		double[] array = new double[ 2 ];
		Arrays.fill( array, value );
		return array;
	}

	public static double[] as3dDoubleArray( double value )
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

	public static RandomAccessibleInterval< BitType > labelRegionAsMask( LabelRegion labelRegion )
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
		final ImgLabeling< Integer, IntType > imgLabeling = new ImgLabeling<>( labelImg );

		final java.util.Iterator< Integer > labelCreator = new java.util.Iterator< Integer >()
		{
			int id = 1;

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

		ConnectedComponents.labelAllConnectedComponents( ( RandomAccessible ) Views.extendBorder( rai ), imgLabeling, labelCreator, ConnectedComponents.StructuringElement.FOUR_CONNECTED );

		return imgLabeling;
	}

	public static ImgLabeling< Integer, IntType > labelMapAsImgLabeling( RandomAccessibleInterval< IntType > labelMap )
	{
		final ImgLabeling< Integer, IntType > imgLabeling = new ImgLabeling<>( labelMap );

		final double maximumLabel = Algorithms.getMaximumValue( labelMap );

		final ArrayList< Set< Integer > > labelSets = new ArrayList< >();

		labelSets.add( new HashSet<>() ); // empty 0 label
		for ( int label = 1; label <= maximumLabel; ++label )
		{
			final HashSet< Integer > set = new HashSet< >();
			set.add( label );
			labelSets.add( set );
		}

		new LabelingMapping.SerialisationAccess< Integer >( imgLabeling.getMapping() )
		{
			{
				super.setLabelSets( labelSets );
			}
		};

		return imgLabeling;
	}

	public static <T extends IntegerType<T> > ImgLabeling< Integer, IntType > labelMapAsImgLabelingRobert( RandomAccessibleInterval< T > labelMap )
	{
		final RandomAccessibleInterval< IntType > indexImg = ArrayImgs.ints( Intervals.dimensionsAsLongArray( labelMap ) );
		final ImgLabeling< Integer, IntType > imgLabeling = new ImgLabeling<>( indexImg );

		final Cursor< LabelingType< Integer > > labelCursor = Views.flatIterable( imgLabeling ).cursor();

		for ( final IntegerType input : Views.flatIterable( labelMap ) ) {

			final LabelingType< Integer > element = labelCursor.next();

			if ( input.getRealFloat() != 0 )
			{
				element.add( (int) input.getRealFloat() );
			}
		}

		return imgLabeling;
	}



	private static Set< Integer > getLabelSet( RandomAccessibleInterval< UnsignedShortType > labelMap )
	{
		final Cursor< UnsignedShortType > cursor = Views.iterable( labelMap ).cursor();

		final Set< Integer > labelSet = new HashSet<>();

		while ( cursor.hasNext() )
		{
			labelSet.add(  cursor.next().getInteger() );
		}

		return labelSet;
	}

	public static void drawObject( RandomAccessibleInterval< IntType > img,
								   LabelRegion labelRegion,
								   int value )
	{
		final Cursor< Void > regionCursor = labelRegion.cursor();
		final RandomAccess< IntType > access = img.randomAccess();
		BitType bitTypeTrue = new BitType( true );
		while ( regionCursor.hasNext() )
		{
			regionCursor.fwd();
			access.setPosition( regionCursor );
			access.get().set( value );
		}
	}

	public static RandomAccessibleInterval<BitType> asMask( ImgLabeling<Integer, IntType> imgLabeling )
	{
		final RandomAccessibleInterval< IntType > labeling = imgLabeling.getSource();

		RandomAccessibleInterval< BitType > mask = asMask( labeling );

		return mask;

	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval<BitType> asMask( RandomAccessibleInterval< T > rai )
	{
		RandomAccessibleInterval< BitType > mask = ArrayImgs.bits( Intervals.dimensionsAsLongArray( rai ) );
		mask = Transforms.getWithAdjustedOrigin( rai, mask  );
		final RandomAccess< BitType > maskAccess = mask.randomAccess();

		final Cursor< T > cursor = Views.iterable( rai ).cursor();

		while ( cursor.hasNext() )
		{
			cursor.fwd();

			if ( cursor.get().getRealDouble() > 0 )
			{
				maskAccess.setPosition( cursor );
				maskAccess.get().set( true );
			}
		}
		return mask;
	}

	public static int getNumObjects( RandomAccessibleInterval< BitType > mask )
	{
		final LabelRegions labelRegions = new LabelRegions( asImgLabeling( mask )  );
		return labelRegions.getExistingLabels().size() - 1;
	}

	public static < T extends RealType< T > & NativeType< T > >
	ImagePlus createIJ1Movie( ArrayList< RandomAccessibleInterval< T > > labelings, String title )
	{
		RandomAccessibleInterval movie = Views.stack( labelings );
		movie = Views.addDimension( movie, 0, 0);
		movie = Views.addDimension( movie, 0, 0);
		movie = Views.permute( movie, 2,4 );
		final ImagePlus imp = new Duplicator().run( ImageJFunctions.wrap( movie, title ) );
		imp.setTitle( title );
		return imp;
	}

	public static boolean acceptFile( String fileNameEndsWith, String file )
	{
		final String[] fileNameEndsWithList = fileNameEndsWith.split( "," );

		for ( String endsWith : fileNameEndsWithList )
		{
			if ( file.endsWith( endsWith.trim() ) )
			{
				return true;
			}
		}

		return false;
	}

	public static void error( String s )
	{
		IJ.showMessage( s );
	}

	public static ImagePlus asLabelImagePlus( RandomAccessibleInterval< IntType > indexImg )
	{
		final Duplicator duplicator = new Duplicator();
		final ImagePlus labelImagePlus = duplicator.run( ImageJFunctions.wrap( indexImg, "mask" ) );
		labelImagePlus.setLut( getGoldenAngleLUT() );
		return labelImagePlus;
	}


	public static LUT getGoldenAngleLUT()
	{
		byte[][] bytes = createGoldenAngleLut( 256 );
		final byte[][] rgb = new byte[ 3 ][ 256 ];

		for ( int c = 0; c < 3; ++c )
		{
			rgb[ c ][ 0 ] = 0; // Black background
		}

		for ( int c = 0; c < 3; ++c )
		{
			for ( int i = 1; i < 256; ++i )
			{
				rgb[ c ][ i ] = bytes[ i ][ c ];
			}
		}

		LUT lut = new LUT( rgb[ 0 ], rgb[ 1 ], rgb[ 2 ] );
		return lut;
	}


	/**
	 * Make lookup table with esthetically pleasing colors based on the golden
	 * angle
	 *
	 * From: MorphoLibJ
	 * // TODO: properly cite!
	 *
	 * @param nColors number of colors to generate
	 * @return lookup table with golden-angled-based colors
	 */
	public final static byte[][] createGoldenAngleLut( int nColors )
	{
		// hue for assigning new color ([0.0-1.0])
		float hue = 0.5f;
		// saturation for assigning new color ([0.5-1.0])
		float saturation = 0.75f;

		// create colors recursively by adding golden angle ratio to hue and
		// saturation of previous color
		Color[] colors = new Color[nColors];
		for (int i = 0; i < nColors; i++)
		{
			// create current color
			colors[i] = Color.getHSBColor(hue, saturation, 1);

			// update hue and saturation for next color
			hue += 0.38197f; // golden angle
			if (hue > 1)
				hue -= 1;
			saturation += 0.38197f; // golden angle
			if (saturation > 1)
				saturation -= 1;
			saturation = 0.5f * saturation + 0.5f;
		}

		// create map
		byte[][] map = new byte[nColors][3];

		// fill up the color map by converting color array
		for (int i = 0; i < nColors; i++)
		{
			Color color = colors[i];
			map[i][0] = (byte) color.getRed();
			map[i][1] = (byte) color.getGreen();
			map[i][2] = (byte) color.getBlue();
		}

		return map;
	}

	public static < T extends AbstractIntegerType< T > >
	Set< Long > computeUniqueValues( RandomAccessibleInterval< T > rai )
	{
		final Set< Long > unique = new HashSet<>(  );

		final Cursor< T > cursor = Views.iterable( rai ).cursor();

		while ( cursor.hasNext() )
		{
			unique.add( cursor.next().getIntegerLong() );
		}

		return unique;
	}

	public static ImagePlus asImagePlus( RandomAccessibleInterval alignedDapiMask )
	{
		return ImageJFunctions.wrap(
				Views.permute(
						Views.addDimension( alignedDapiMask, 0, 0 ),
						2, 3), "aligned_dapi_mask" );
	}
}
