package de.embl.cba.morphometry.tracking;

import de.embl.cba.morphometry.Utils;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.HashMap;

public class TrackingUtils
{
	public static void addOverlap( HashMap< Integer, Long > overlaps, int integer )
	{
		if ( overlaps.keySet().contains( integer ) )
		{
			overlaps.put( integer,  overlaps.get( integer ).longValue() + 1 );
		}
		else
		{
			overlaps.put( integer,  1L );
		}
	}

	public static int computeObjectId( HashMap< Integer, Long > overlaps, int nextId )
	{
		int objectId;

		if( overlaps.size() == 0 )
		{
			objectId = nextId++;
		}
		else
		{
			objectId = getMaxOverlapLabel( overlaps );
		}

		return objectId;
	}

	public static int getMaxOverlapLabel( HashMap< Integer, Long > overlaps )
	{
		int maxOverlapLabel = 0;
		long maxOverlapCount = Long.MIN_VALUE;

		for ( Integer label : overlaps.keySet() )
		{
			final long overlapCount = overlaps.get( label ).longValue();

			if ( overlapCount > maxOverlapCount )
			{
				maxOverlapCount = overlapCount;
				maxOverlapLabel = label;
			}
		}
		return maxOverlapLabel;
	}

	public static HashMap< Integer, Long > computeOverlaps(
			RandomAccessibleInterval< IntType > previousLabeling,
			LabelRegion region )
	{
		final HashMap< Integer, Long > overlaps = new HashMap<>();

		final RandomAccess< IntType > previousLabelsAccess = previousLabeling.randomAccess();
		final LabelRegionCursor cursor = region.cursor();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			previousLabelsAccess.setPosition( cursor );

			final int label = previousLabelsAccess.get().getInteger();

			if ( label != 0 )
			{
				addOverlap( overlaps, label );
			}
		}
		return overlaps;
	}

	public static int getNumObjects( RandomAccessibleInterval< BitType > mask )
	{
		final LabelRegions labelRegions = new LabelRegions( Utils.asImgLabeling( mask )  );
		return labelRegions.getExistingLabels().size() - 1;
	}
}
