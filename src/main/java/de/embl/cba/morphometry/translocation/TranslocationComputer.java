package de.embl.cba.morphometry.translocation;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.regions.Regions;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import java.util.ArrayList;

public class TranslocationComputer< T extends RealType< T > & NativeType< T > >
{
	final ArrayList< RandomAccessibleInterval< T > > movie;
	final ArrayList<FinalInterval> intervals;

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

		minimalObjectSize = 5;

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
			computeTranslocation( getCroppedImageAtTimePoint( interval, t ), result );
		}

		results.add( result );
	}

	private IntervalView< T > getCroppedImageAtTimePoint( FinalInterval interval, int t )
	{
		return Views.interval(
			Views.dropSingletonDimensions( movie.get( t ) ),
				interval );
	}

	private void computeTranslocation( RandomAccessibleInterval< T > image, TranslocationResult result )
	{
		final RandomAccessibleInterval< BitType > mask = createCellMask( image );
		result.cellMasks.add( mask );
	}

	private RandomAccessibleInterval< BitType > createCellMask( RandomAccessibleInterval< T > image )
	{
		final Histogram1d< T > histogram = opService.image().histogram( Views.iterable( image ) );
		double huang = opService.threshold().huang( histogram ).getRealDouble();
		//double yen = opService.threshold().yen( histogram ).getRealDouble();
		final RandomAccessibleInterval< BitType > mask = Algorithms.createMask( image, huang );
		Regions.removeSmallRegionsInMask( mask, minimalObjectSize );
		opService.morphology().fillHoles( mask );
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
