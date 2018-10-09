package de.embl.cba.morphometry;

import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.morphology.Closing;
import net.imglib2.algorithm.morphology.Dilation;
import net.imglib2.algorithm.morphology.Opening;
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

	public static final int WATERSHED = -1;

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createRescaledArrayImg(
			RandomAccessibleInterval< T > input,
			double[] scalingFactors )
	{
		assert scalingFactors.length == input.numDimensions();

		/**
		 * - In principle, writing a function that computes weighted averages
		 *   of an appropriate number of neighboring (not only nearest) pixels
		 *   around each requested (real) position in the new image appears to me
		 *   the most straight-forward way of rescaling.
		 * - However, in practice, blurring and subsequent re-sampling seems to be more commonly done,
		 *   maybe for implementation efficiency?
		 * - http://imagej.1557.x6.nabble.com/downsampling-methods-td3690444.html
		 * - https://github.com/axtimwalde/mpicbg/blob/050bc9110a186394ea15190fd326b3e32829e018/mpicbg/src/main/java/mpicbg/ij/util/Filter.java#L424
		 * - https://imagej.net/Downsample
		 */

		/*
		 * Blur image
		 */

		final RandomAccessibleInterval< T > blurred = createOptimallyBlurredArrayImg( input, scalingFactors );

		/*
		 * Sample values from blurred image
		 */

		final RandomAccessibleInterval< T > resampled = createResampledArrayImg( blurred, scalingFactors );

		return resampled;
	}

	private static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createResampledArrayImg( RandomAccessibleInterval< T > input, double[] scalingFactors )
	{
		// Convert to RealRandomAccessible such that we can obtain values at (infinite) non-integer coordinates
		RealRandomAccessible< T > rra = Views.interpolate( Views.extendBorder( input ), new NLinearInterpolatorFactory<>() );

		// Change scale such that we can sample from integer coordinates (for raster function below)
		Scale scale = new Scale( scalingFactors );
		RealRandomAccessible< T > rescaledRRA  = RealViews.transform( rra, scale );

		// Create view sampled at integer coordinates
		final RandomAccessible< T > rastered = Views.raster( rescaledRRA );

		// Put an interval to make it a finite "normal" image again
		final RandomAccessibleInterval< T > finiteRastered = Views.interval( rastered, createTransformedInterval( input, scale ) );

		// Convert from View to a "conventional" image in RAM
		// - Above code would also run on, e.g. 8 TB image, within ms
		// - Now, we really force it to create the image (we actually might now have to, depends...)
		final RandomAccessibleInterval< T > output = Utils.copyAsArrayImg( finiteRastered );

		return output;
	}

	private static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > createOptimallyBlurredArrayImg(
			RandomAccessibleInterval< T > input,
			double[] scalingFactors )
	{
		/**
		 * - https://en.wikipedia.org/wiki/Decimation_(signal_processing)
		 * - Optimal blurring is 0.5 / M, where M is the downsampling factor
		 */

		final double[] sigmas = new double[input.numDimensions() ];

		for ( int d = 0; d < input.numDimensions(); ++d )
		{
			sigmas[ d ] = 0.5 / scalingFactors[ d ];
		}

		// allocate output image
		RandomAccessibleInterval< T > output = Utils.createEmptyArrayImg( input );

		// blur input image and write into output image
		Gauss3.gauss( sigmas, Views.extendBorder( input ), output ) ;

		return output;
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
	RandomAccessibleInterval< BitType > fillHoles3Din2D(
			RandomAccessibleInterval< BitType > mask,
			int axis,
			OpService opService )
	{
		final ArrayList< RandomAccessibleInterval< BitType > > holesFilled = new ArrayList<>();

		final IntervalView< BitType > rotated = Views.rotate( mask, axis, 2 );

		for ( long coordinate = rotated.min( 2 ); coordinate <= rotated.max( 2 ); ++coordinate )
		{
			RandomAccessibleInterval< BitType > maskSlice = Views.hyperSlice( rotated, 2, coordinate );
			holesFilled.add( opService.morphology().fillHoles( maskSlice ) );
		}

		RandomAccessibleInterval< BitType > stack =  Views.stack( holesFilled );

		stack = Views.zeroMin( Views.rotate( stack, 2, axis ) );

		stack = Transforms.getWithAdjustedOrigin( mask, stack );

		return stack;

	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createWatershedSeeds(
			RandomAccessibleInterval< T > distance,
			Shape shape,
			double globalThreshold,
			double localThreshold )
	{

		RandomAccessibleInterval< BitType > seeds = ArrayImgs.bits( Intervals.dimensionsAsLongArray( distance ) );
		seeds = Transforms.getWithAdjustedOrigin( distance, seeds );

		RandomAccessible< Neighborhood< T > > neighborhoods = shape.neighborhoodsRandomAccessible( Views.extendPeriodic( distance ) );
		RandomAccessibleInterval< Neighborhood< T > > neighborhoodsInterval = Views.interval( neighborhoods, distance );

		final Cursor< Neighborhood< T > > neighborhoodCursor = Views.iterable( neighborhoodsInterval ).cursor();
		final RandomAccess< T > distanceRandomAccess = distance.randomAccess();
		final RandomAccess< BitType > seedsRandomAccess = seeds.randomAccess();


		double[] currentPosition = new double[ distance.numDimensions() ];

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
			else if ( Utils.isLateralBoundaryPixel( neighborhood, distance ) && distanceRandomAccess.get().getRealDouble() >  0 )
			{
				seedsRandomAccess.get().set( true );
			}
			else if ( isCenterLargestOrEqual( centerValue, neighborhood ) )
			{
				if ( centerValue.getRealDouble() > localThreshold )
				{
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
	ArrayList< PositionAndValue > getLocalMaxima(
			RandomAccessibleInterval< T > rai,
			double minimalDistanceBetweenMaxima,
			double threshold )
	{

		Shape shape = new HyperSphereShape( (long) minimalDistanceBetweenMaxima );
		RandomAccessible< Neighborhood< T > > neighborhoods = shape.neighborhoodsRandomAccessible( Views.extendPeriodic( rai ) );
		final Cursor< Neighborhood< T > > neighborhoodCursor = Views.iterable( Views.interval( neighborhoods, rai ) ).cursor();
		final RandomAccess< T > randomAccess = rai.randomAccess();

		final ArrayList< PositionAndValue > maxima = new ArrayList<>();

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

					boolean isMaximumFarEnoughAwayFromOtherMaxima = true;

					for ( int i = 0; i < maxima.size(); ++i )
					{
						double distance = LinAlgHelpers.distance( maxima.get( i ).position, positionAndValue.position );
						if ( distance < minimalDistanceBetweenMaxima )
						{
							isMaximumFarEnoughAwayFromOtherMaxima = false;
							break;
						}
					}

					if ( isMaximumFarEnoughAwayFromOtherMaxima )
					{
						maxima.add( positionAndValue );
					}

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
			ImgLabeling< Integer, IntType > imgLabeling,
			RandomAccessibleInterval< T > intensity,
			RandomAccessibleInterval< BitType > mask, // This will be changed, i.e. the split(s) will be drawn into it
			HashMap< Integer, Integer > numObjectsPerRegion,
			long minimalObjectWidth,
			long minimalObjectSize,
			long maximalWatershedBoundaryLength,
			OpService opService,
			boolean forceSplit,
			boolean showSplittingAttempts )
	{

		final LabelRegions labelRegions = new LabelRegions( imgLabeling );

		for ( int label : numObjectsPerRegion.keySet() )
		{
			if ( numObjectsPerRegion.get( label ) > 1 )
			{

				final RandomAccessibleInterval< T > maskedAndCropped = Views.zeroMin( Utils.getMaskedAndCropped( intensity, labelRegions.getLabelRegion( label ) ) );
				final RandomAccessibleInterval< BitType > labelRegionMask = Views.zeroMin( Utils.labelRegionAsMask( labelRegions.getLabelRegion( label ) ) );

				final ArrayList< PositionAndValue > localMaxima =
						computeLocalIntensityMaxima( 2 * minimalObjectWidth, maskedAndCropped, showSplittingAttempts );

				if ( localMaxima.size() < 2 )
				{
					Utils.log( "Not enough local maxima found for object: " + label );
					continue; // TODO: check these cases
				}

				final RandomAccessibleInterval< T > seeds = getSeeds( numObjectsPerRegion.get( label ), maskedAndCropped, localMaxima );

				//ImageJFunctions.show( seeds, "seeds" );

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

				if ( ! splitObjects.getExistingLabels().contains( -1 ) )
				{
					continue; // TODO: examine these cases
				}


				boolean isValidSplit;

				if ( forceSplit )
				{
					isValidSplit = true;
				}
				else
				{
					// TODO: add integrated intensity along watershed as criterium
					isValidSplit = checkSplittingValidity(
							splitObjects,
							minimalObjectSize,
							maximalWatershedBoundaryLength );
				}

				Utils.log( "Valid split found: " + isValidSplit );

				if ( showSplittingAttempts )
				{
					ImageJFunctions.show( watershedImgLabeling.getSource(), "" + label + "-" + isValidSplit );

				}

				if ( isValidSplit )
				{
					drawWatershedIntoMask( mask, labelRegions, label, splitObjects );
				}

			}
		}

	}

	public static < T extends RealType< T > & NativeType< T > >
	ArrayList< PositionAndValue > computeLocalIntensityMaxima(
			long minimalObjectWidth,
			RandomAccessibleInterval< T > maskedAndCropped,
			boolean showSplittingAttempts )
	{
		double blurSimga = minimalObjectWidth / 2.0;
		final RandomAccessibleInterval< T > blurred = Utils.createBlurredRai( maskedAndCropped, blurSimga, 1.0 );
		//if ( showSplittingAttempts ) ImageJFunctions.show( blurred, "blurred" );
		final ArrayList< PositionAndValue > sortedLocalMaxima =
				Algorithms.getLocalMaxima(
						blurred,
						 minimalObjectWidth,
						0.0 );
		return sortedLocalMaxima;
	}

	private static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > getSeeds( int numSeeds,
											RandomAccessibleInterval< T > maskedAndCropped,
											ArrayList< PositionAndValue > localMaxima )
	{
		final RandomAccessibleInterval< T > seeds = Utils.createEmptyArrayImg( maskedAndCropped );
		final RandomAccess< T > randomAccess = seeds.randomAccess();
		for ( int i = 0; i < numSeeds; ++i )
		{
			randomAccess.setPosition( Utils.asLongs( localMaxima.get( i ).position ) );
			randomAccess.get().setOne();
		}
		return seeds;
	}

	public static boolean checkSplittingValidity(
			LabelRegions< Integer > splitObjects,
			long minimumObjectSize,
			long maximalWatershedLength )
	{

		if ( ! isWatershedValid( splitObjects, maximalWatershedLength ) ) return false;

		ArrayList< Long > regionSizes = new ArrayList<>(  );

		for( LabelRegion region : splitObjects )
		{
			regionSizes.add( region.size() );
		}


		if ( regionSizes.size() >=2 )
		{
			Collections.sort( regionSizes );
			Collections.reverse( regionSizes );

			if ( regionSizes.get( 1 ) < minimumObjectSize )
			{
				// 2nd largest object too small
				return false;
			}
		}

		return true;
	}

	public static boolean isWatershedValid( LabelRegions< Integer > splitObjects, long maximalWatershedLength )
	{
		boolean isValidSplit = true;

		for( LabelRegion region : splitObjects )
		{
			int splitObjectLabel = ( int ) region.getLabel();

			if ( splitObjectLabel == WATERSHED )
			{
				final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( Utils.labelRegionAsMask( region ) );
				final LabelRegions< Integer > splitRegions = new LabelRegions( imgLabeling );

				long maximalLength = 0;
				for ( LabelRegion splitRegion : splitRegions )
				{
					if ( splitRegion.size() > maximalLength )
					{
						maximalLength = splitRegion.size();
					}
				}

				if ( maximalLength > maximalWatershedLength )
				{
					isValidSplit = false;
					break;
				}
			}
		}
		return isValidSplit;
	}

	public static void drawWatershedIntoMask( RandomAccessibleInterval< BitType > mask,
											  LabelRegions labelRegions,
											  int label,
											  LabelRegions< Integer > splitObjects )
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

	public static < T extends RealType< T > & NativeType< T > > ArrayList< RandomAccessibleInterval< T > >
	createMaximumProjectedIntensitiesAssumingImagePlusDimensionOrder(
			RandomAccessibleInterval< T > inputImages,
			long c,
			long tMin, long tMax )
	{
		ArrayList<  RandomAccessibleInterval< T > > intensities = new ArrayList<>();

		for ( long t = tMin; t <= tMax; ++t )
		{
			final IntervalView< T > channelView = Views.hyperSlice( inputImages, 2, c );
			final IntervalView< T > timePointView = Views.hyperSlice( channelView, 3, t );
			final RandomAccessibleInterval maximum = new Projection( timePointView, 2 ).maximum();
			intensities.add( maximum );
		}

		return intensities;
	}

	public static RandomAccessibleInterval< DoubleType > computeDistanceTransform( RandomAccessibleInterval< BitType > mask )
	{
		final RandomAccessibleInterval< DoubleType > doubleBinary = Converters.convert( mask, ( i, o ) -> o.set( i.get() ? Double.MAX_VALUE : 0 ), new DoubleType() );

		final RandomAccessibleInterval< DoubleType > distance = ArrayImgs.doubles( Intervals.dimensionsAsLongArray( doubleBinary ) );

		DistanceTransform.transform( doubleBinary, distance, DistanceTransform.DISTANCE_TYPE.EUCLIDIAN, 1.0D );
		return distance;
	}

	public static RandomAccessibleInterval< BitType > open(
			RandomAccessibleInterval< BitType > mask,
			int radius )
	{
		// TODO: Bug(?!) in imglib2 Closing.close makes this necessary
		RandomAccessibleInterval< BitType > morphed = ArrayImgs.bits( Intervals.dimensionsAsLongArray( mask ) );
		final RandomAccessibleInterval< BitType > enlargedMask = Utils.getEnlargedRai2( mask, 2 * radius );
		final RandomAccessibleInterval< BitType > enlargedMorphed = Utils.getEnlargedRai2( morphed, 2 * radius );

		if ( radius > 0 )
		{
			Utils.log( "Morphological opening...");
			Shape shape = new HyperSphereShape( radius );
			Opening.open( Views.extendZero( enlargedMask ), Views.iterable( enlargedMorphed ), shape, 1 );
		}

		return Views.interval( enlargedMorphed, mask );
	}

	public static RandomAccessibleInterval< BitType > dilate(
			RandomAccessibleInterval< BitType > mask,
			int radius )
	{
		RandomAccessibleInterval< BitType > morphed = ArrayImgs.bits( Intervals.dimensionsAsLongArray( mask ) );

		if ( radius > 0 )
		{
			Shape shape = new HyperSphereShape( radius );
			Dilation.dilate( Views.extendZero( mask ), Views.iterable( morphed ), shape, 1 );
		}

		return morphed;
	}
}
