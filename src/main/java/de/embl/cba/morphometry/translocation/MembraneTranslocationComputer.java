package de.embl.cba.morphometry.translocation;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.CoordinateAndValue;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.regions.RegionAndSize;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.morphometry.segmentation.SignalOverBackgroundSegmenter;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.converter.Converters;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

import java.util.ArrayList;
import java.util.Collections;

public class MembraneTranslocationComputer< R extends RealType< R > & NativeType< R > >
{
	final ArrayList< RandomAccessibleInterval< R > > movie;
	final ArrayList<FinalInterval> intervals;
	private double signalToNoise;
	private int closingRadius;
	private int resolutionBlurWidthInPixel;
	private int blurRadius;
	private int erosionRadius;
	private int region;

	public ArrayList< TranslocationResult > getResults()
	{
		return results;
	}

	final ArrayList<TranslocationResult> results;

	final OpService opService;
	private int minimalRegionSize;

	public MembraneTranslocationComputer(
			ArrayList< RandomAccessibleInterval< R > > movie,
			ArrayList< FinalInterval > intervals,
			OpService opService )
	{
		this.movie = movie;
		this.intervals = intervals;
		this.opService = opService;

		minimalRegionSize = 50;
		signalToNoise = 25.0;
		closingRadius = 5;
		resolutionBlurWidthInPixel = 3;
		blurRadius = 4;
		erosionRadius = 4;

		results = new ArrayList<TranslocationResult>();

		computeTranslocations();
	}

	private void computeTranslocations()
	{
		for ( region = 0; region < intervals.size(); region++ )
		{
			computeTranslocationsInRoi( intervals.get( region ) );
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

			measureMembraneIntensity( result, t );

			createInsideOutsideMasksAndMeasureIntensities( result, t );

			computeTranslocation( result, t );
		}

		results.add( result );
	}

	private void computeTranslocation( TranslocationResult result, int t )
	{
		double membrane = ( double ) result.membraneIntensities.get( t );
		double outside = ( double ) result.outsideIntensities.get( t );
		double inside = ( double ) result.insideIntensities.get( t );

		if ( inside > 0 || outside > 0 )
		{
			final double translocation = ( membrane - outside ) / ( inside - outside );
			result.translocations.add( translocation );
		}
		else
		{
			// something went wrong
			result.translocations.add( -1.0 );
		}
	}

	private void createMembraneMask( FinalInterval interval, TranslocationResult result, int t )
	{
		final RandomAccessibleInterval< R > intensity = Views.interval( movie.get( t ) , interval );

//		final RandomAccessibleInterval< R > smoothed =
//				opService.filter().gauss(
//						intensity,
//						1.0 );

		RandomAccessibleInterval< R > smoothed = Algorithms.median(
				intensity, 4, opService );

		final RandomAccessibleInterval< R > erode = Algorithms.erode(
				smoothed, erosionRadius );

		RandomAccessibleInterval< R > gradient = Utils.createEmptyCopy( erode );

		LoopBuilder.setImages( gradient, smoothed, erode ).forEachPixel( ( g, i, e ) ->
				{
					g.setReal( i.getRealDouble() - e.getRealDouble() );
				}
		);

		RandomAccessibleInterval< BitType > binaryGradient = ArrayImgs.bits( Intervals.dimensionsAsLongArray( interval ) );
		binaryGradient = Views.translate( binaryGradient, Intervals.minAsLongArray( interval ) );

		opService.threshold().isoData(
				Views.iterable( binaryGradient ),
				Views.iterable( gradient ) );

		Regions.onlyKeepLargestRegion( binaryGradient );

		// smooth for helping the thinning
		final RandomAccessibleInterval< BitType > openedBinaryGradient
				= Algorithms.dilate( Algorithms.erode( binaryGradient, 1 ), 1 );

		final RandomAccessibleInterval< BitType > thin = Algorithms.thin( openedBinaryGradient, opService );

		result.binaryGradients.add( openedBinaryGradient );
		result.intensities.add( smoothed );
		result.gradients.add( gradient );
		result.membranes.add( thin );
	}

	private void measureMembraneIntensity( TranslocationResult result, int t )
	{
		final Cursor< BitType > cursor = Views.iterable(
				(RandomAccessibleInterval) result.membranes.get( t ) ).cursor();
		final RandomAccess< R > access = movie.get( t ).randomAccess();

		double sum = 0.0;
		int n = 0;
		while ( cursor.hasNext() )
		{
			if ( cursor.next().get() )
			{
				access.setPosition( cursor );
				sum += access.get().getRealDouble();
				n++;
			}
		}

		result.membraneIntensities.add( sum / n );
	}

	private void createInsideOutsideMasksAndMeasureIntensities( TranslocationResult result, int t )
	{

		final RandomAccessibleInterval< BitType > membraneMask = ( RandomAccessibleInterval< BitType > ) result.membranes.get( t );

		final RandomAccessibleInterval< BitType > invertedMembraneMask = Converters.convert(
				Algorithms.dilate( membraneMask, 2 * resolutionBlurWidthInPixel ), // dilate to exclude membrane blur from inside or outside
				( i, o ) -> o.set( !i.get() ),
				new BitType() );

		final ArrayList< RegionAndSize > sortedRegions = Regions.getSizeSortedRegions( invertedMembraneMask );

		final RandomAccessibleInterval< BitType > insideOutsideMask = Utils.createEmptyCopy( membraneMask );

		final RandomAccess< R > intensityAccess = movie.get( t ).randomAccess();

		final ArrayList< Double > intensities = new ArrayList<>();

		if ( sortedRegions.size() < 2 )
		{
			// Something went wrong
			Logger.log( "Could not detect membrane in region " + ( results.size() + 1 ) + ", frame " + ( t + 1) );
//			ImageJFunctions.show( ( RandomAccessibleInterval ) result.gradients.get( t ), "Failure: Gradient Timepoint " + (t+1) );
			intensities.add( -1.0 );
			intensities.add( -1.0 );
		}
		else
		{
			for ( int region = 0; region < 2; region++ )
			{
				intensities.add( getMean( sortedRegions, intensityAccess, region ) );
				Regions.drawRegionInMask( sortedRegions.get( region ).getRegion(), insideOutsideMask );
			}

			Collections.sort( intensities );
		}

		result.outsideIntensities.add( intensities.get( 0 ) ); // the smaller one
		result.insideIntensities.add( intensities.get( 1 ) ); // the larger one
		result.insideOutsideMasks.add( insideOutsideMask );

	}

	private double getMean( ArrayList< RegionAndSize > regionsAndSizes, RandomAccess< R > access, int i )
	{
		final LabelRegionCursor cursor = regionsAndSizes.get( i ).getRegion().cursor();

		double sum = 0.0;
		int n = 0;

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			access.setPosition( cursor );
			sum += access.get().getRealDouble();
			n++;
		}

		return sum / n;
	}

	private void computeIntensities( RandomAccessibleInterval< R > image, TranslocationResult result )
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
			RandomAccessibleInterval< R > image,
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
			final RandomAccessibleInterval< R > image,
			final RandomAccessibleInterval< IntType > distances,
			final int dMin,
			final int dMax )
	{
		final Cursor< IntType > cursor = Views.iterable( distances ).cursor();
		final RandomAccess< R > access = image.randomAccess();

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

	private RandomAccessibleInterval< BitType > createCellMask( RandomAccessibleInterval< R > image )
	{
		final SignalOverBackgroundSegmenter< R > segmenter =
				new SignalOverBackgroundSegmenter<>(
						image,
						signalToNoise,
						minimalRegionSize );

		RandomAccessibleInterval< BitType > mask = segmenter.createMask();

		mask = Algorithms.close( mask, closingRadius );

		return mask;
	}

}
