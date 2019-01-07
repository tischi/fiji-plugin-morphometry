package de.embl.cba.morphometry.translocation;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.segmentation.SignalOverBackgroundSegmenter;
import net.imagej.ops.OpService;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.converter.Converters;
import net.imglib2.img.array.ArrayImgFactory;
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
import java.util.Collections;
import java.util.HashMap;

public class TranslocationComputer< T extends RealType< T > & NativeType< T > >
{
	final ArrayList< RandomAccessibleInterval< T > > movie;
	final ArrayList<FinalInterval> intervals;
	private double signalToNoise;
	private int closingRadius;
	private int maxDistanceOfMembraneToOutside;

	public ArrayList< TranslocationResult > getResults()
	{
		return results;
	}

	final ArrayList<TranslocationResult> results;

	final OpService opService;
	private int minimalObjectSize;

	public TranslocationComputer(
			ArrayList< RandomAccessibleInterval< T > > movie,
			ArrayList< FinalInterval > intervals,
			OpService opService )
	{
		this.movie = movie;
		this.intervals = intervals;
		this.opService = opService;

		minimalObjectSize = 50;
		signalToNoise = 25.0;
		closingRadius = 5;
		maxDistanceOfMembraneToOutside = 6;

		results = new ArrayList<TranslocationResult>();

		computeTranslocations();
	}

	public void computeTranslocations()
	{
		for ( FinalInterval interval : intervals )
		{
			computeTranslocationsInRoi( interval );
		}
	}

	private void computeTranslocationsInRoi( FinalInterval interval )
	{
		final TranslocationResult result = new TranslocationResult();

		for ( int t = 0; t < movie.size(); t++ )
		{
			computeIntensities( Views.interval( movie.get( t ), interval ), result );

			double membrane = ( double ) result.membraneIntensities.get( t );
			double outside = ( double ) result.outsideIntensities.get( t );
			double inside = ( double ) result.insideIntensities.get( t );

			result.translocation.add( ( membrane - outside ) / ( inside - outside ) );
		}

		results.add( result );
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

		result.membraneIntensities.add( computeMembraneIntensity( image, distances ) );

		result.insideIntensities.add( computeMeanIntensityWithinDistanceRange(
				image,
				distances,
				maxDistanceOfMembraneToOutside,
				Integer.MAX_VALUE ) );
	}

	private Double computeMembraneIntensity( RandomAccessibleInterval< T > image, RandomAccessibleInterval< IntType > distances )
	{
		final HashMap< Integer, Double > distanceIntensityMap = new HashMap<>();

		for ( int d = 1; d < maxDistanceOfMembraneToOutside; d++ )
		{
			double mean = computeMeanIntensityWithinDistanceRange( image, distances, d, d );

			distanceIntensityMap.put( d, mean );
		}

		return Collections.max( distanceIntensityMap.values() );
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
			RandomAccessibleInterval< T > image,
			RandomAccessibleInterval< IntType > distances,
			int dMin,
			int dMax )
	{
		final Cursor< IntType > cursor = Views.iterable( distances ).cursor();
		final RandomAccess< T > access = image.randomAccess();

		double sum = 0.0;
		int n = 0;

		while ( cursor.hasNext() )
		{
			if ( cursor.next().getInteger() >= dMin && cursor.next().getInteger() <= dMax )
			{
				access.setPosition( cursor );
				sum += access.get().getRealDouble();
				n++;
			}
		}

		if ( n > 0 ) return sum / n;
		else return -1;
	}

	private RandomAccessibleInterval< BitType > createCellMask( RandomAccessibleInterval< T > image )
	{
		final SignalOverBackgroundSegmenter< T > segmenter =
				new SignalOverBackgroundSegmenter<>(
						image,
						signalToNoise,
						minimalObjectSize );

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
