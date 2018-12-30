package de.embl.cba.morphometry.regions;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.transforms.utils.Transforms;
import net.imglib2.Cursor;
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

import java.util.ArrayList;

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

	public static < T extends RealType< T > & NativeType< T > >
	void removeSmallRegionsInMasks(
			ArrayList< RandomAccessibleInterval< T > > masks,
			double sizeInCalibratedUnits,
			double calibration )
	{

		for ( RandomAccessibleInterval mask : masks )
		{
			removeSmallRegionsInMask( mask, sizeInCalibratedUnits, calibration );
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	void removeSmallRegionsInMask(
			RandomAccessibleInterval< T > mask,
			double sizeInCalibratedUnits,
			double calibration )
	{
		final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( mask );

		long minimalObjectSize = ( long ) ( sizeInCalibratedUnits / Math.pow( calibration, imgLabeling.numDimensions() ) );

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			if ( labelRegion.size() < minimalObjectSize )
			{
				removeRegion( mask, labelRegion );
			}
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	void removeSmallRegionsInMasks(
			ArrayList< RandomAccessibleInterval< T > > masks,
			long minimalObjectSize )
	{

		for ( RandomAccessibleInterval mask : masks )
		{
			removeSmallRegionsInMask( mask, minimalObjectSize );
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	void removeSmallRegionsInMask(
			RandomAccessibleInterval< T > mask,
			long minimalObjectSize )
	{
		final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( mask );

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );
		for ( LabelRegion labelRegion : labelRegions )
		{
			if ( labelRegion.size() < minimalObjectSize )
			{
				removeRegion( mask, labelRegion );
			}
		}
	}

	private static void drawRegion( RandomAccessibleInterval< BitType > img, LabelRegion labelRegion )
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

	private static < T extends RealType< T > & NativeType< T > >
	void removeRegion( RandomAccessibleInterval< T > img, LabelRegion labelRegion )
	{
		final Cursor regionCursor = labelRegion.cursor();
		final RandomAccess< T > access = img.randomAccess();
		while ( regionCursor.hasNext() )
		{
			regionCursor.fwd();
			access.setPosition( regionCursor );
			access.get().setReal( 0 );
		}
	}
}
