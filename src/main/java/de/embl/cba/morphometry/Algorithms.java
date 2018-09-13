package de.embl.cba.morphometry;

import de.embl.cba.morphometry.objects.Measurements;
import net.imglib2.*;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
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
import net.imglib2.util.Intervals;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.Views;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
	void setValues( RandomAccessibleInterval< T > rai, double value )
	{
		Cursor< T > cursor = Views.iterable( rai ).localizingCursor();

		double maxValue = Double.MIN_VALUE;

		while ( cursor.hasNext() )
		{
			cursor.next().setReal( value );
		}
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

		ImgLabeling< Integer, IntType > labelImg = createImgLabeling( sizeFilteredObjects );
		return labelImg;
	}


	private static RandomAccessibleInterval< BitType > removeSmallObjectsAndReturnMask( ImgLabeling< Integer, IntType > labeling, double size, double calibration )
	{
		RandomAccessibleInterval< BitType > sizeFilteredObjectsMask = ArrayImgs.bits( Intervals.dimensionsAsLongArray( labeling ) );
		sizeFilteredObjectsMask = Transforms.getWithAdjustedOrigin( labeling.getSource(), sizeFilteredObjectsMask  );

		long pixelSize = ( long ) ( size / Math.pow( calibration, labeling.numDimensions() ) );

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( labeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			if ( labelRegion.size() > pixelSize )
			{
				drawObject( sizeFilteredObjectsMask, labelRegion );
			}
		}

		return sizeFilteredObjectsMask;
	}

	public static RandomAccessibleInterval< BitType >
	removeSmallObjectsAndReturnMask( RandomAccessibleInterval< BitType > img, double size, double calibration )
	{
		return removeSmallObjectsAndReturnMask( createImgLabeling( img ), size, calibration );
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

	public static < T extends RealType< T > & NativeType< T >  >
	ImgLabeling< Integer, IntType > createImgLabeling( RandomAccessibleInterval< T > rai )
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

}
