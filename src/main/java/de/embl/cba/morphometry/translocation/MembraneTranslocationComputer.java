package de.embl.cba.morphometry.translocation;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.CoordinateAndValue;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.morphometry.segmentation.SignalOverBackgroundSegmenter;
import net.imagej.ops.OpService;
import net.imagej.ops.Ops;
import net.imglib2.*;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.converter.Converters;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
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
import java.util.Collection;
import java.util.Collections;

public class MembraneTranslocationComputer< T extends RealType< T > & NativeType< T > >
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

	public MembraneTranslocationComputer(
			ArrayList< RandomAccessibleInterval< T > > movie,
			ArrayList< FinalInterval > intervals,
			OpService opService )
	{
		this.movie = movie;
		this.intervals = intervals;
		this.opService = opService;

		minimalRegionSize = 10;
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
			createMembraneMaskAndMeasureIntensity( interval, result, t );

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

		final double translocation = ( membrane - outside ) / ( inside - outside );

		if ( translocation < 0 )
		{
			int b = 1;
		}

		result.translocations.add( translocation );
	}

	private void createMembraneMaskAndMeasureIntensity( FinalInterval interval, TranslocationResult result, int t )
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

		final RandomAccessibleInterval< BitType > thin = Utils.createEmptyArrayImg( mask );
		opService.morphology().thinGuoHall( thin, mask );

		measureMembraneIntensity( result, t, thin );

		result.membraneMasks.add( thin );
	}

	private void measureMembraneIntensity( TranslocationResult result, int t, RandomAccessibleInterval< BitType > thin )
	{
		final Cursor< BitType > cursor = Views.iterable( thin ).cursor();
		final RandomAccess< T > access = movie.get( t ).randomAccess();

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

	class RegionAndSize implements Comparable< RegionAndSize >
	{
		LabelRegion region;
		Long size;

		public RegionAndSize( LabelRegion region, Long size )
		{
			this.region = region;
			this.size = size;
		}

		public int compareTo( RegionAndSize r )
		{
			return ( (int) ( r.size - this.size ) );
		}
	}

	private void createInsideOutsideMasksAndMeasureIntensities( TranslocationResult result, int t )
	{

		final RandomAccessibleInterval< BitType > membraneMask = ( RandomAccessibleInterval< BitType > ) result.membraneMasks.get( t );

		final RandomAccessibleInterval< BitType > invertedMembraneMask = Converters.convert(
				Algorithms.dilate( membraneMask, 2 ),
				( i, o ) -> o.set( !i.get() ),
				new BitType() );

		final ArrayList< RegionAndSize > sortedRegions = getSizeSortedRegions( invertedMembraneMask );

		final RandomAccessibleInterval< BitType > insideOutsideMask = Utils.createEmptyArrayImg( membraneMask );

		final RandomAccess< T > intensityAccess = movie.get( t ).randomAccess();

		final ArrayList< Double > intensities = new ArrayList<>();

		for ( int region = 0; region < 2; region++ )
		{
			intensities.add( getMean( sortedRegions, intensityAccess, region ) );
			Regions.drawRegionInMask( sortedRegions.get( region ).region, insideOutsideMask );
		}

		Collections.sort( intensities );

		result.outsideIntensities.add( intensities.get( 0 ) ); // the smaller one
		result.insideIntensities.add( intensities.get( 1 ) ); // the larger one

		result.insideOutsideMasks.add( insideOutsideMask );

	}

	private ArrayList< RegionAndSize > getSizeSortedRegions( RandomAccessibleInterval< BitType > invertedMembraneMask )
	{
		final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( invertedMembraneMask );

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );

		final ArrayList< RegionAndSize > regionsAndSizes = new ArrayList<>();

		for ( LabelRegion labelRegion : labelRegions )
		{
			regionsAndSizes.add( new RegionAndSize( labelRegion, labelRegion.size() ) );
		}

		Collections.sort( regionsAndSizes );

		return regionsAndSizes;
	}

	private double getMean( ArrayList< RegionAndSize > regionsAndSizes, RandomAccess< T > access, int i )
	{
		final LabelRegionCursor cursor = regionsAndSizes.get( i ).region.cursor();

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
