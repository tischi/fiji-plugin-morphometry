package de.embl.cba.morphometry.regions;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.transforms.utils.Transforms;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public abstract class Regions
{

	public static LabelRegion< Integer > getCentralRegion( ImgLabeling< Integer, IntType > labeling, long radius )
	{
		int centralLabel = Algorithms.getCentralLabelIndex( labeling, radius );

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( labeling );

		return labelRegions.getLabelRegion( centralLabel );
	}

	public static LabelRegion< Integer > getLargestRegion( ImgLabeling< Integer, IntType > labeling )
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

	public static long size( LabelRegion labelRegion )
	{

		final LabelRegionCursor cursor = labelRegion.cursor();

		long size = 0;
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			size++;
		}

		return size;
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
}
