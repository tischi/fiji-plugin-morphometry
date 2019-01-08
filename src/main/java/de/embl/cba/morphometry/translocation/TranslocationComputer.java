package de.embl.cba.morphometry.translocation;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.CoordinateAndValue;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.morphometry.segmentation.SignalOverBackgroundSegmenter;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.converter.Converters;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.LongArray;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

import java.util.ArrayList;

public class TranslocationComputer< T extends RealType< T > & NativeType< T > >
{
	final ArrayList< RandomAccessibleInterval< T > > movie;
	final ArrayList<FinalInterval> intervals;
	private double signalToNoise;
	private int closingRadius;
	private int maxDistanceOfMembraneToOutside;
	private int resolutionBlurWidthInPixel;
	private int minimumGradient;

	public ArrayList< TranslocationResult > getResults()
	{
		return results;
	}

	final ArrayList<TranslocationResult> results;

	final OpService opService;
	private int minimalRegionSize;

	public TranslocationComputer(
			ArrayList< RandomAccessibleInterval< T > > movie,
			ArrayList< FinalInterval > intervals,
			OpService opService )
	{
		this.movie = movie;
		this.intervals = intervals;
		this.opService = opService;

		minimalRegionSize = 50;
		signalToNoise = 25.0;
		closingRadius = 5;
		minimumGradient = 20;

		resolutionBlurWidthInPixel = 2;
		maxDistanceOfMembraneToOutside = Integer.MAX_VALUE; // TODO: how to cope with the ruffling? maybe find cell boundary with a gradient

		results = new ArrayList<TranslocationResult>();

		computeTranslocations();
	}

	private void computeTranslocations()
	{
		for ( FinalInterval interval : intervals )
		{
			computeTranslocationsInRoi( interval );
		}
	}

	private void computeTranslocationsInRoi( FinalInterval interval )
	{
		final TranslocationResult result = new TranslocationResult();

		result.regionCenterX = Utils.getCenterLocation( interval )[ 0 ];
		result.regionCenterY = Utils.getCenterLocation( interval )[ 1 ];

		for ( int t = 0; t < movie.size(); t++ )
		{

			createMembraneMask( interval, result, t );

//			Algorithms.createLocalMaximaMask(
//						(RandomAccessibleInterval< T >) result.gradients.get( t ),
//						shape,
//						minimumGradient ) );

			computeIntensities( Views.interval( movie.get( t ), interval ), result );

			double membrane = ( double ) result.membraneIntensities.get( t );
			double outside = ( double ) result.outsideIntensities.get( t );
			double inside = ( double ) result.insideIntensities.get( t );


			final double translocation = ( membrane - outside ) / ( inside - outside );

			if ( translocation < 0 )
			{
				int b = 1;
			}

			result.translocations.add( translocation );

		}

		results.add( result );
	}

	private void createMembraneMask( FinalInterval interval, TranslocationResult result, int t )
	{
		final RandomAccessibleInterval< T > gauss = (RandomAccessibleInterval) opService.filter().gauss( Views.interval( movie.get( t ), interval ), 1 );

		final HyperSphereShape shape = new HyperSphereShape( resolutionBlurWidthInPixel );

		result.gradients.add( Algorithms.computeGradient( gauss, shape ) );

		RandomAccessibleInterval< BitType > mask = ArrayImgs.bits( Intervals.dimensionsAsLongArray( interval ) );
		mask = Views.translate( mask, Intervals.minAsLongArray( interval ) );

		opService.threshold().isoData(
				Views.iterable( mask ),
				Views.iterable( ( RandomAccessibleInterval< T > ) result.gradients.get( t ) ) );

		Regions.removeSmallRegionsInMask( mask, minimalRegionSize );

		result.membraneMasks.add( mask );
	}

	private void computeIntensities( RandomAccessibleInterval< T > image, TranslocationResult result )
	{
		final RandomAccessibleInterval< BitType > mask = createCellMask( image );
		result.cellMasks.add( mask );

		final RandomAccessibleInterval< IntType > distances = computeDistanceMap( mask );

		result.outsideIntensities.add( computeMeanIntensityWithinDistanceRange(
				image,
				distances,
				0,
				0 ) );

		final CoordinateAndValue membraneDistanceAndIntensity =
				computeMembraneIntensity( image, distances );

		result.membraneIntensities.add( membraneDistanceAndIntensity.value );

		result.insideIntensities.add( computeMeanIntensityWithinDistanceRange(
				image,
				distances,
				(int) membraneDistanceAndIntensity.coordinate + resolutionBlurWidthInPixel,
				Integer.MAX_VALUE ) );
	}

	private CoordinateAndValue computeMembraneIntensity(
			RandomAccessibleInterval< T > image,
			RandomAccessibleInterval< IntType > distances )
	{
		double maximalGradient = 0;
		double membraneIntensity = 0;
		double distanceAtMaximalGradient = 0;

		ArrayList< Double > intensities = new ArrayList<>(  );

		final double dMax = Algorithms.getMaximumValue( distances );
		for ( int d = 0; d <= dMax; d++ )
		{
			double intensity = computeMeanIntensityWithinDistanceRange( image, distances, d, d );

			intensities.add( intensity );

			int previousDistance = d > resolutionBlurWidthInPixel ? d - resolutionBlurWidthInPixel : 0;

			double gradient = intensity - intensities.get( previousDistance );

			if ( gradient > maximalGradient )
			{
				maximalGradient = gradient;
				distanceAtMaximalGradient = d;
				membraneIntensity = intensity;
			}

		}

		return new CoordinateAndValue( distanceAtMaximalGradient, membraneIntensity );
	}

	private RandomAccessibleInterval< IntType > computeDistanceMap( RandomAccessibleInterval< BitType > mask )
	{
		final RandomAccessibleInterval< DoubleType > squaredDistances = Algorithms.computeSquaredDistances( mask );

		return Converters.convert(
				squaredDistances,
				( i, o ) -> o.set( (int) Math.sqrt( i.get() ) ),
				new IntType() );
	}

	private double computeMeanIntensityWithinDistanceRange(
			final RandomAccessibleInterval< T > image,
			final RandomAccessibleInterval< IntType > distances,
			final int dMin,
			final int dMax )
	{
		final Cursor< IntType > cursor = Views.iterable( distances ).cursor();
		final RandomAccess< T > access = image.randomAccess();

		double sum = 0.0;
		int n = 0;

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			if ( cursor.get().getInteger() >= dMin && cursor.get().getInteger() <= dMax )
			{
				access.setPosition( cursor );
				sum += access.get().getRealDouble();
				n++;
			}
		}

		if ( n == 0 )
		{
			final double maximumValue = Algorithms.getMaximumValue( distances );
			int a = 1;
		}

		return sum / n;

	}

	private RandomAccessibleInterval< BitType > createCellMask( RandomAccessibleInterval< T > image )
	{
		final SignalOverBackgroundSegmenter< T > segmenter =
				new SignalOverBackgroundSegmenter<>(
						image,
						signalToNoise,
						minimalRegionSize );

		RandomAccessibleInterval< BitType > mask = segmenter.createMask();

		mask = Algorithms.close( mask, closingRadius );

		return mask;
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > copyAsArrayImg( RandomAccessibleInterval< T > orig )
	{
		final RandomAccessibleInterval< T > copy =
				Views.translate(
						new ArrayImgFactory( Util.getTypeFromInterval( orig ) ).create( orig ),
						Intervals.minAsLongArray( orig ) );

		LoopBuilder.setImages( copy, orig ).forEachPixel( Type::set );

		return copy;
	}


}
