package de.embl.cba.morphometry;

import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.morphology.Closing;
import net.imglib2.algorithm.morphology.distance.DistanceTransform;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.realtransform.RealViews;
import net.imglib2.realtransform.Scale;
import net.imglib2.roi.labeling.*;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import java.util.*;

import static de.embl.cba.morphometry.Constants.XYZ;
import static de.embl.cba.morphometry.Transforms.createTransformedInterval;
import static java.lang.Math.abs;
import static java.lang.Math.acos;

public class Algorithms
{

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createIsotropicArrayImg( RandomAccessibleInterval< T > input, double[] scalingFactors )
	{
		assert scalingFactors.length == input.numDimensions();

		/*
		 * Blur image
		 */

		RandomAccessibleInterval< T > blurred = createOptimallyBlurredArrayImg( input, scalingFactors );

		/*
		 * Sample values from blurred image
		 */

		final RandomAccessibleInterval< T > resampled = createResampledArrayImg( blurred, scalingFactors );

		return resampled;
	}

	private static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createResampledArrayImg( RandomAccessibleInterval< T > input, double[] scalingFactors )
	{
		// TODO: is there a simple way to do this?
		Scale scale = new Scale( scalingFactors );
		RealRandomAccessible< T > rra = Views.interpolate( Views.extendBorder( input ), new NLinearInterpolatorFactory<>() );
		rra = RealViews.transform( rra, scale );
		final RandomAccessible< T > raster = Views.raster( rra );
		final RandomAccessibleInterval< T > output = Views.interval( raster, createTransformedInterval( input, scale ) );

		return Utils.copyAsArrayImg( output  );
	}

	private static < T extends RealType< T > & NativeType< T > > RandomAccessibleInterval< T > createOptimallyBlurredArrayImg( RandomAccessibleInterval< T > input, double[] scalingFactors )
	{
		final long[] inputDimensions = Intervals.dimensionsAsLongArray( input );
		final double[] sigmas = new double[ inputDimensions.length ];

		for ( int d = 0; d < inputDimensions.length; ++d )
		{
			sigmas[ d ] = 0.5 / scalingFactors[ d ]; // From Saalfeld
		}

		ImgFactory< T > imgFactory = new ArrayImgFactory( input.randomAccess().get()  );
		RandomAccessibleInterval< T > blurred = Views.translate( imgFactory.create( inputDimensions ), Intervals.minAsLongArray( input )  ) ;

		Gauss3.gauss( sigmas, Views.extendBorder( input ), blurred ) ;

		return blurred;
	}


	public static < T extends RealType< T > & NativeType< T > >
	Point findMaximumLocation( RandomAccessibleInterval< T > rai, double[] calibration )
	{
		Cursor< T > cursor = Views.iterable( rai ).localizingCursor();

		double maxValue = Double.MIN_VALUE;

		long[] maxLoc = new long[ cursor.numDimensions() ];
		cursor.localize( maxLoc );

		while ( cursor.hasNext() )
		{
			final double value = cursor.next().getRealDouble();
			if ( value > maxValue )
			{
				maxValue = value;
				cursor.localize( maxLoc );
			}
		}


		for ( int d = 0; d < rai.numDimensions(); ++d )
		{
			maxLoc[ d ] *= calibration[ d ];
		}
		
		Point point = new Point( maxLoc );

		return point;
	}

	public static < T extends RealType< T > & NativeType< T > >
	double getMaximumValue( RandomAccessibleInterval< T > rai )
	{
		Cursor< T > cursor = Views.iterable( rai ).localizingCursor();

		double maxValue = Double.MIN_VALUE;

		while ( cursor.hasNext() )
		{
			final double value = cursor.next().getRealDouble();
			if ( value > maxValue )
			{
				maxValue = value;
			}
		}

		return maxValue;
	}

	public static < T extends RealType< T > & NativeType< T > >
	boolean isCenterLargest( T center, Neighborhood< T > neighborhood )
	{
		boolean centerIsLargest = true;

		for( T neighbor : neighborhood ) {
			if( neighbor.compareTo( center ) > 0 )
			{
				centerIsLargest = false;
				break;
			}
		}

		return centerIsLargest;
	}

	public static < T extends RealType< T > & NativeType< T > >
	List< RealPoint > findLocalMaximumValues( RandomAccessibleInterval< T > rai, Shape shape )
	{
		List< RealPoint > points = new ArrayList<>();

		RandomAccessible<Neighborhood<T>> neighborhoods = shape.neighborhoodsRandomAccessible( Views.extendBorder( rai ) );
		RandomAccessibleInterval<Neighborhood<T>> neighborhoodsInterval = Views.interval( neighborhoods, rai );

		LoopBuilder.setImages( neighborhoodsInterval, rai ).forEachPixel(
				(neighborhood, center) -> {
					if( isCenterLargest( center, neighborhood ) )
					{
						points.add( new RealPoint( neighborhood ) );
					}
				}
		);

		return points;
	}


	public static int getCentralLabelIndex( ImgLabeling< Integer, IntType > labeling )
	{
		final RandomAccess< LabelingType< Integer > > labelingRandomAccess = labeling.randomAccess();
		for ( int d : XYZ ) labelingRandomAccess.setPosition( labeling.dimension( d ) / 2, d );
		int centralIndex = labelingRandomAccess.get().getIndex().getInteger();
		return labeling.getMapping().labelsAtIndex( centralIndex ).iterator().next();
	}

	public static LabelRegion< Integer > getCentralObjectLabelRegion( ImgLabeling< Integer, IntType > labeling )
	{
		int centralLabel = getCentralLabelIndex( labeling );

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( labeling );

		return labelRegions.getLabelRegion( centralLabel );
	}

	public static LabelRegion< Integer > getLargestObject( ImgLabeling< Integer, IntType > labeling )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( labeling );

		long maxSize = Long.MIN_VALUE;
		LabelRegion largestRegion = null;

		for ( LabelRegion labelRegion : labelRegions )
		{
			if ( labelRegion.size() > maxSize)
			{
				largestRegion = labelRegion;
				maxSize = labelRegion.size();
			}
		}

		return largestRegion;
	}


	public static Img< UnsignedByteType > createUnsignedByteTypeMaskFromLabelRegion( LabelRegion< Integer > centralObjectRegion, long[] dimensions )
	{
		final Img< UnsignedByteType > centralObjectImg = ArrayImgs.unsignedBytes( dimensions );

		final Cursor< Void > regionCursor = centralObjectRegion.cursor();
		final RandomAccess< UnsignedByteType > access = centralObjectImg.randomAccess();
		while ( regionCursor.hasNext() )
		{
			regionCursor.fwd();
			access.setPosition( regionCursor );
			access.get().set( 255 );
		}
		return centralObjectImg;
	}

	public static ImgLabeling< Integer, IntType > removeSmallObjectsAndReturnImgLabeling( ImgLabeling< Integer, IntType > labeling, double size, double calibration )
	{
		RandomAccessibleInterval< BitType > sizeFilteredObjects = removeSmallObjectsAndReturnMask( labeling, size, calibration );

		ImgLabeling< Integer, IntType > labelImg = Utils.asImgLabeling( sizeFilteredObjects );
		return labelImg;
	}


	private static RandomAccessibleInterval< BitType > removeSmallObjectsAndReturnMask( ImgLabeling< Integer, IntType > labeling, double size, double calibration )
	{
		RandomAccessibleInterval< BitType > sizeFilteredObjectsMask = ArrayImgs.bits( Intervals.dimensionsAsLongArray( labeling ) );
		sizeFilteredObjectsMask = Transforms.getWithAdjustedOrigin( labeling.getSource(), sizeFilteredObjectsMask  );

		long minimalObjectPixelSize = ( long ) ( size / Math.pow( calibration, labeling.numDimensions() ) );

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( labeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			if ( labelRegion.size() > minimalObjectPixelSize )
			{
				drawObject( sizeFilteredObjectsMask, labelRegion );
			}
		}

		return sizeFilteredObjectsMask;
	}

	public static RandomAccessibleInterval< BitType >
	removeSmallObjectsAndReturnMask( RandomAccessibleInterval< BitType > img, double size, double calibration )
	{
		return removeSmallObjectsAndReturnMask( Utils.asImgLabeling( img ), size, calibration );
	}

	private static void drawObject( RandomAccessibleInterval< BitType > img, LabelRegion labelRegion )
	{
		final Cursor< Void > regionCursor = labelRegion.cursor();
		final RandomAccess< BitType > access = img.randomAccess();
		BitType bitTypeTrue = new BitType( true );
		while ( regionCursor.hasNext() )
		{
			regionCursor.fwd();
			access.setPosition( regionCursor );
			access.get().set( bitTypeTrue );
		}
	}


	public static Img< BitType > createMaskFromLabelRegion( LabelRegion< Integer > centralObjectRegion, long[] dimensions )
	{
		final Img< BitType > centralObjectImg = ArrayImgs.bits( dimensions );

		final Cursor< Void > regionCursor = centralObjectRegion.cursor();
		final RandomAccess< BitType > access = centralObjectImg.randomAccess();
		while ( regionCursor.hasNext() )
		{
			regionCursor.fwd();
			access.setPosition( regionCursor );
			access.get().set( true );
		}
		return centralObjectImg;
	}


	public static ArrayList< RealPoint > origin()
	{
		final ArrayList< RealPoint > origin = new ArrayList<>();
		origin.add( new RealPoint( new double[]{ 0, 0, 0 } ) );
		return origin;
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

	public static double angleInRadians( final double[] v1, final double[] v2 )
	{
		final double dot = LinAlgHelpers.dot( v1, v2 );

		final double angleInRadians = Math.signum( v1[ 0 ] ) * acos( dot / ( LinAlgHelpers.length( v1 ) * LinAlgHelpers.length( v2 ) ) );
		final double angleInDegrees = angleInRadians * 180.0 / Math.PI;

		return angleInRadians;
	}

	public static < T extends RealType< T > & NativeType< T > >
	void copy( RandomAccessibleInterval< T > source, RandomAccessibleInterval< T > target )
	{
		final Cursor< T > cursor = Views.iterable( source ).cursor();
		final RandomAccess< T > randomAccess = target.randomAccess();

		while( cursor.hasNext() )
		{
			cursor.fwd();
			randomAccess.setPosition( cursor );
			randomAccess.get().set( cursor.get() );
		}

	}


	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createMask( RandomAccessibleInterval< T > downscaled, double threshold )
	{
		Utils.log( "Creating mask...");

		RandomAccessibleInterval< BitType > mask =
				Converters.convert( downscaled, ( i, o )
						-> o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

		return mask;
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createSeeds( RandomAccessibleInterval< T > distance, Shape shape, double globalThreshold, double localThreshold )
	{

		RandomAccessibleInterval< BitType > seeds = ArrayImgs.bits( Intervals.dimensionsAsLongArray( distance ) );
		seeds = Transforms.getWithAdjustedOrigin( distance, seeds );

		RandomAccessible< Neighborhood< T > > neighborhoods = shape.neighborhoodsRandomAccessible( Views.extendPeriodic( distance ) );
		RandomAccessibleInterval< Neighborhood< T > > neighborhoodsInterval = Views.interval( neighborhoods, distance );

		final Cursor< Neighborhood< T > > neighborhoodCursor = Views.iterable( neighborhoodsInterval ).cursor();
		final RandomAccess< T > distanceRandomAccess = distance.randomAccess();
		final RandomAccess< BitType > seedsRandomAccess = seeds.randomAccess();

		double[] centerPosition = new double[ distance.numDimensions() ];
		double[] currentPosition = new double[ distance.numDimensions() ];

		for ( int d = 0; d < centerPosition.length; ++d )
		{
			centerPosition[ d ] = (long) (distance.dimension( d ) / 2);
		}

		while ( neighborhoodCursor.hasNext() )
		{
			final Neighborhood< T > neighborhood = neighborhoodCursor.next();
			neighborhood.localize( currentPosition );
			seedsRandomAccess.setPosition( neighborhood );
			distanceRandomAccess.setPosition( neighborhood );

			T centerValue = distanceRandomAccess.get();

			if ( centerValue.getRealDouble() > globalThreshold )
			{
				seedsRandomAccess.get().set( true );
			}
			else if ( LinAlgHelpers.distance( currentPosition, centerPosition  ) < 3 )
			{
				// seedsRandomAccess.get().set( true );
			}
			else if ( Utils.isLateralBoundaryPixel( neighborhood, distance ) && distanceRandomAccess.get().getRealDouble() >  0 )
			{
				seedsRandomAccess.get().set( true );
			}
			else if ( isCenterLargestOrEqual( centerValue, neighborhood ) )
			{
				if ( centerValue.getRealDouble() > localThreshold )
				{
					// local maximum and larger than local Threshold
					seedsRandomAccess.get().set( true );
				}
			}

		}

		return seeds;
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createCenterAndBoundarySeeds( RandomAccessibleInterval< T > distance, Shape shape, double globalThreshold, double localThreshold )
	{

		RandomAccessibleInterval< BitType > seeds = ArrayImgs.bits( Intervals.dimensionsAsLongArray( distance ) );
		seeds = Transforms.getWithAdjustedOrigin( distance, seeds );

		RandomAccessible< Neighborhood< T > > neighborhoods = shape.neighborhoodsRandomAccessible( Views.extendPeriodic( distance ) );
		RandomAccessibleInterval< Neighborhood< T > > neighborhoodsInterval = Views.interval( neighborhoods, distance );

		final Cursor< Neighborhood< T > > neighborhoodCursor = Views.iterable( neighborhoodsInterval ).cursor();
		final RandomAccess< T > distanceRandomAccess = distance.randomAccess();
		final RandomAccess< BitType > seedsRandomAccess = seeds.randomAccess();

		double[] centerPosition = new double[ distance.numDimensions() ];
		double[] currentPosition = new double[ distance.numDimensions() ];

		for ( int d = 0; d < centerPosition.length; ++d )
		{
			centerPosition[ d ] = (long) (distance.dimension( d ) / 2);
		}

		while ( neighborhoodCursor.hasNext() )
		{
			final Neighborhood< T > neighborhood = neighborhoodCursor.next();
			neighborhood.localize( currentPosition );
			seedsRandomAccess.setPosition( neighborhood );
			distanceRandomAccess.setPosition( neighborhood );

			T centerValue = distanceRandomAccess.get();

			if ( centerValue.getRealDouble() > globalThreshold )
			{
				// maximaRandomAccess.get().set( true );
			}
			else if ( LinAlgHelpers.distance( currentPosition, centerPosition  ) < 3 )
			{
				seedsRandomAccess.get().set( true );
			}
			else if ( Utils.isLateralBoundaryPixel( neighborhood, distance ) && distanceRandomAccess.get().getRealDouble() >  0 )
			{
				seedsRandomAccess.get().set( true );
			}
//			else if ( isCenterLargestOrEqual( centerValue, neighborhood ) )
//			{
//				if ( centerValue.getRealDouble() > localThreshold )
//				{
//					// local maximum and larger than local Threshold
//					// maximaRandomAccess.get().set( true );
//				}
//			}

		}

		return seeds;
	}


	public static < T extends RealType< T > & NativeType< T > >
	ArrayList< PositionAndValue > getLocalMaxima( RandomAccessibleInterval< T > rai, Shape shape, double threshold )
	{
		final ArrayList< PositionAndValue > maxima = new ArrayList<>();

		RandomAccessible< Neighborhood< T > > neighborhoods = shape.neighborhoodsRandomAccessible( Views.extendPeriodic( rai ) );
		final Cursor< Neighborhood< T > > neighborhoodCursor = Views.iterable( Views.interval( neighborhoods, rai ) ).cursor();
		final RandomAccess< T > randomAccess = rai.randomAccess();

		while ( neighborhoodCursor.hasNext() )
		{
			final Neighborhood< T > neighborhood = neighborhoodCursor.next();
			randomAccess.setPosition( neighborhood );

			T centerValue = randomAccess.get();

			if ( isCenterLargestOrEqual( centerValue, neighborhood ) )
			{
				if ( centerValue.getRealDouble() > threshold )
				{
					final PositionAndValue positionAndValue = new PositionAndValue();
					positionAndValue.position = new double[ rai.numDimensions() ];
					neighborhood.localize( positionAndValue.position );
					positionAndValue.value = centerValue.getRealDouble();
					maxima.add( positionAndValue );
				}
			}
		}

		maxima.sort( Comparator.comparing( PositionAndValue::getValue ).reversed() );

		return maxima;
	}

	private static < T extends RealType< T > & NativeType< T > >
	boolean isCenterLargestOrEqual( T center, Neighborhood< T > neighborhood )
	{
		for( T neighbor : neighborhood )
		{
			if( neighbor.compareTo( center ) > 0 )
			{
				return false;
			}
		}
		return true;
	}


	public static < T extends RealType< T > & NativeType< T > >
	void splitTouchingObjects(
			RandomAccessibleInterval< BitType > mask,
			ImgLabeling< Integer, IntType > imgLabeling,
			HashMap< Integer, Integer > numObjectsPerRegion,
			RandomAccessibleInterval< T > image,
			long minimalObjectPixelWidth,
			long minimalObjectSize,
			long maximalWatershedLength,
			OpService opService )
	{

		final LabelRegions labelRegions = new LabelRegions( imgLabeling );

		for ( int label : numObjectsPerRegion.keySet() )
		{
			if ( numObjectsPerRegion.get( label ) > 1 )
			{

				final RandomAccessibleInterval< T > maskedAndCropped = Views.zeroMin( Utils.getMaskedAndCropped( image, labelRegions.getLabelRegion( label ) ) );
				final RandomAccessibleInterval< BitType > labelRegionMask = Views.zeroMin( Utils.asMask( labelRegions.getLabelRegion( label ) ) );

				final ArrayList< PositionAndValue > localMaxima = Algorithms.getLocalMaxima( maskedAndCropped, new HyperSphereShape( minimalObjectPixelWidth ), 0.0 );

				final RandomAccessibleInterval< T > seeds = Utils.copyAsEmptyArrayImg( maskedAndCropped );
				final RandomAccess< T > randomAccess = seeds.randomAccess();
				for ( int i = 0; i < numObjectsPerRegion.get( label ); ++i )
				{
					randomAccess.setPosition( Utils.asLongs( localMaxima.get( i ).position ) );
					randomAccess.get().setOne();
				}

				final ImgLabeling< Integer, IntType > watershedImgLabeling = getEmptyImgLabeling( labelRegionMask );
				final ImgLabeling< Integer, IntType > seedsImgLabeling = Utils.asImgLabeling( seeds );

				opService.image().watershed(
						watershedImgLabeling,
						Utils.invertedView( maskedAndCropped ),
						seedsImgLabeling,
						true,
						true,
						labelRegionMask );


				LabelRegions< Integer > splitObjects = new LabelRegions( watershedImgLabeling );

				boolean isValidSplit = checkSplittingValidity( splitObjects, minimalObjectSize, maximalWatershedLength );

				if ( isValidSplit )
				{
					drawWatershedIntoMask( mask, labelRegions, label, splitObjects );
				}

				// Utils.applyMask( watershedImgLabeling.getSource(), mask );

//				ImageJFunctions.show( seedsImgLabeling.getSource() );
				ImageJFunctions.show( watershedImgLabeling.getSource() );

			}
		}

	}

	public static boolean checkSplittingValidity(  LabelRegions< Integer > splitObjects, long minimumObjectSize, long maximalWatershedLength )
	{
		boolean isValidSplit = true;

		for( LabelRegion region : splitObjects )
		{
			int splitObjectLabel = ( int ) region.getLabel();
			if ( splitObjectLabel == -1 )
			{
				if ( region.size() > maximalWatershedLength )
				{
					isValidSplit = false;
					break;
				}
			}
			else
			{
				if ( region.size() < minimumObjectSize )
				{
					isValidSplit = false;
					break;
				}
			}
		}

		return isValidSplit;
	}

	public static void drawWatershedIntoMask( RandomAccessibleInterval< BitType > mask, LabelRegions labelRegions, int label, LabelRegions< Integer > splitObjects )
	{
		final long[] regionOffset = Intervals.minAsLongArray( labelRegions.getLabelRegion( label ) );
		LabelRegion watershed = splitObjects.getLabelRegion( -1 );
		final LabelRegionCursor cursor = watershed.cursor();
		final RandomAccess< BitType > maskRandomAccess = mask.randomAccess();
		long[] position = new long[ watershed.numDimensions() ];
		while( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( position );
			for ( int d = 0; d < position.length; ++d )
			{
				position[ d ] += regionOffset[ d ];
			}
			maskRandomAccess.setPosition( position );
			maskRandomAccess.get().set( false );
		}
	}

	public static RandomAccessibleInterval< BitType > close( RandomAccessibleInterval< BitType > mask, int closingRadius )
	{
		RandomAccessibleInterval< BitType > closed = Utils.copyAsArrayImg( mask );
		Shape closingShape = new HyperSphereShape( closingRadius );
		Closing.close( Views.extendBorder( mask ), Views.iterable( closed ), closingShape, 1 );
		return closed;
	}


	public static ImgLabeling< Integer, IntType > getEmptyImgLabeling( RandomAccessibleInterval< BitType > mask )
	{
		RandomAccessibleInterval< IntType > watershedLabelImg = ArrayImgs.ints( Intervals.dimensionsAsLongArray( mask ) );
		watershedLabelImg = Transforms.getWithAdjustedOrigin( mask, watershedLabelImg );
		return new ImgLabeling<>( watershedLabelImg );
	}
}
